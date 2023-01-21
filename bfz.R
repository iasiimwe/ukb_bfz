# Libraries and relevant functions
# --------------------------------
library(data.table) # reading tables faster
library(tidyverse) # data processing 
library(qqman) # for generating Manhattan and QQ plots
library(gridExtra) # for arranging multiple ggplots on the same page
library(grid)
library(rms) # assessing non-linearity using restricted cubic splines
library(pwr) # for sample size calculations
source("relevant_functions.R")


# Minimum sample size
# --------------------
# Based on one of the primary outcomes (Blood glucose). Assuming:
  # 80% power
  # significance threshold of 5e-08
  # Effect size of 1.1 mmol/L (10% of 11.1 mmol/L, the hyperglycaemia threshold, PMID: 27631769)
  # Standard deviation of 1.2 mmol/L (for 431 764 UK Biobank participants)
# Different MAFs would require the following minimum sample sizes:
ukb_min_sample(1.1, 1.2, maf = 0.2, n_times = 2) # 20% MAF - 310 participants
ukb_min_sample(1.1, 1.2, maf = 0.1, n_times = 3) # 10% MAF - 539 participants
ukb_min_sample(1.1, 1.2, maf = 0.05, n_times = 5) # 5% MAF - 1008 participants
ukb_min_sample(1.1, 1.2, maf = 0.01, n_times = 24) # 5% MAF - 4776 participants


# Hypertensive participants
# -------------------------
# Import converted UK Biobank dataset
bd <- as_tibble(fread("bd_enc.csv.gz")) 
length(unique(bd$f.eid)) # 502413 UK Biobank participants in total

withdrawn_consent <- as_tibble(fread("withdrawn_consent.csv")) # this is emailed to
# approved users or can be found accessed through the UK Biobank Access Management System.
bd <- bd %>% filter(!f.eid %in% withdrawn_consent$V1)
length(unique(bd$f.eid)) # 502406 UK Biobank participants with consent

# Turn into long-format for diseases to get HTN at baseline (requires fields 53 and 20002)
# In the codings, 1065 = hypertension and 1072 = essential hypertension. 1073 (gestational hypertension/pre-eclampsia) not included
dict <- as_tibble(fread("Data_Dictionary_Showcase.csv")) # import dictionary
bd_HTN <- select(bd, contains(c("f.eid", "f.53.", "f.20002."))) %>%
  ukb_reshape_long(dict) %>%
  filter(`Non-cancer illness code, self-reported` == 1065 | `Non-cancer illness code, self-reported` == 1072) %>%
  filter(I == ".0.") %>%
  select(f.eid) %>%
  unique()
nrow(bd_HTN) # 133271 participants with self-reported hypertension at baseline


# Exclude patients taking specific drugs (only done for sensitivity analysis)
# ---------------------------------------------------------------------------
# Did not include conditions such as gout or diabetes since not all patients will be on medications

# Get medications (required field is 20003)
bd_meds <- select(bd, contains(c("f.eid", "f.53.", "f.20003."))) %>%
  filter(f.eid %in% bd_HTN$f.eid) %>%
  ukb_reshape_long(dict) %>%
  rename(Date = `Date of attending assessment centre`,
         Coding = `Treatment/medication code`) %>%
  filter(I == ".0.")
dict_meds <- as_tibble(fread("medication_GWAS_dict.csv")) # Derived from PMID: 31015401
bd_meds <- left_join(bd_meds, dict_meds) %>%
  select(-Date, -Coding) # Combine bd_meds with dict_meds and removing columns that are no longer needed

# Antidiabetics - relevant ATC code is (A10) 'Antidiabetes preparations'
ATC_code <- "(^|\\|)A10"
bd_meds_diabetes <- bd_meds %>%
  filter(str_detect(.$Medication_ATC_code, ATC_code))
length(unique(bd_meds_diabetes$f.eid)) # 12340 participants
  
# Anti-gout - relevant ATC code is (M04A) 'Antigout preparations' but colchicine has no effect on uric acid metabolism.
# Consider only allopurinol or sulfinpyrazone (PMID: 35128470)
# Allopurinol inhibits uric acid production, sulfinpyrazone increases uric acid excretion
bd_meds <- select(bd, contains(c("f.eid", "f.53.", "f.20003."))) %>%
  filter(f.eid %in% bd_HTN$f.eid) %>%
  ukb_reshape_long(dict)
dict_meds <- as_tibble(fread("Codings.csv")) # codings for medications == 4
dict_meds <- dict_meds %>%
  subset(Coding == 4) %>%
  select(Value, Meaning)
colnames(dict_meds) <- c("Treatment/medication code","Treatment/medication")
dict_meds$`Treatment/medication code` <-
  parse_integer(dict_meds$`Treatment/medication code`)
bd_meds <- bd_meds %>%
  left_join(dict_meds) %>%
  select(!`Treatment/medication code`) %>%
  filter(I == ".0.") 
drug <- "allopurinol|sulfinpyrazone"
index <- str_detect(tolower(bd_meds$`Treatment/medication`), drug)
index[is.na(index)] <- FALSE
bd_meds_gout <- bd_meds[index, ] %>%
  select(f.eid) %>%
  unique()
nrow(bd_meds_gout) # 3537 patients taking these (3517 allopurinol, 20 sulfinpyrazone)

# Potassium supplements
drug <- "\\+potassium|potassium product|potassium citrate 3g|potassium bicarbonate\\+citric acid|potassium aminobenzoate|movicol oral powder|sando-k effervescent tablet|cymalon cranberry 1.5g/|slow-k 600mg m|effercitrate soluble tablet"
index <- str_detect(tolower(bd_meds$`Treatment/medication`), drug)
index[is.na(index)] <- FALSE
bd_meds_potassium <- bd_meds[index, ] %>%
  select(f.eid) %>%
  unique()
nrow(bd_meds_potassium) # 914 patients taking potassium supplements

# Patients not taking any of the above drugs
with_excluded_drugs <- unique(c(bd_meds_diabetes$f.eid, bd_meds_gout$f.eid, bd_meds_potassium$f.eid))
bd_HTN <- bd_HTN %>%
  filter(!f.eid %in% with_excluded_drugs) 
nrow(bd_HTN) # 117133 patients not taking drugs affecting any of the outcomes


# Participant quality control checks (informed consent, with genetic data etc)
# -----------------------------------------------------------------------------
# Requires fields 22027 (outliers for heterozygosity or missing rate), 53 (date), 31 (sex), 
# 22001 (genetic sex), 22006 (genetic ethnic grouping) and 22019 (sex chromosome aneuploidy)

bd_check <- select(bd, contains(c("f.eid", "f.53.", "f.31.", "f.22001.", "f.22006."))) %>%
  filter(f.eid %in% bd_HTN$f.eid) %>%
  ukb_reshape_long(dict) %>%
  filter(I == ".0.")
nrow(bd_check) # 133271 participants

# Exclude those who did not cluster with Caucasians
bd_caucasian <- bd_check %>%
  filter(`Genetic ethnic grouping` == "Caucasian") %>%
  select(-`Genetic ethnic grouping`)
nrow(bd_caucasian) # 108850 participants

# Exclude those with sex-discrepancies 
bd_sex_checked <- bd_caucasian %>%
  filter(Sex == `Genetic sex`) %>%
  select(f.eid) 
nrow(bd_sex_checked) # 108753 participants

# Sex chromosome aneuploidy (was downloaded separately)
x_aneuploidy <- as_tibble(fread("x_aneuploidy.csv.gz")) %>%
  filter(is.na(`22019-0.0`) == FALSE) %>%
  select(eid)
bd_sex_checked_ane <- bd_sex_checked %>%
  filter(!f.eid %in% x_aneuploidy$eid) 
nrow(bd_sex_checked_ane) # 108632 participants

# Exclude heterozygosity/missing rate outliers (was downloaded separately)
outliers <- as_tibble(fread("outliers.csv")) %>% 
  filter(is.na(f.22027.0.0) == FALSE) 
bd_not_outliers <- bd_sex_checked_ane %>%
  filter(!f.eid %in% outliers$f.eid) 
nrow(bd_not_outliers) # 108453 participants 

# Relatedness 
missingness <- as_tibble(fread("missingness.csv"))
related <- as_tibble(fread("UKBB_relatedness.csv")) %>%
  filter(Kinship > 0.0884) %>%
  select(ID1, ID2)

related_ids <- unique(c(related$ID1, related$ID2))
related_ids <- related_ids[related_ids %in% bd_not_outliers$f.eid] #to save time, limit
# to only those ids present in the bd_not_outliers dataframe, n = 15,920

# First deal with the ids that appear once (2 relatives)
related.df <- tibble(ID1 = vector("integer", length(related_ids)),
                     include = vector("character", length(related_ids))
                     )
for (i in seq_along(related_ids)) {
  id <- related_ids[[i]]
  id_df <- related %>%
    filter(ID1 == id | ID2 == id)
  
  all_ids <- id_df %>%
    as.matrix() %>%
    as.integer()
  
  id_df_3_way <- related %>%
    filter(ID1 %in% all_ids | ID2 %in% all_ids) # this should capture 2+-way
  # relations
  
  all_IDS <- id_df_3_way %>%
    as.matrix() %>%
    as.integer()
  
  id_df_3_way <- related %>%
    filter(ID1 %in% all_IDS | ID2 %in% all_IDS) # additional iteration in case
  # of multiple relations
  
  id_count <- id_df_3_way %>%
    as.matrix() %>%
    as.integer() %>%
    as_tibble() %>%
    count(value)
  
  related.df[i, 1] <- id
  
  if(nrow(id_df_3_way) == 1) {  # one pair of related individuals
    
    related.df[i, 2] <- ukb_related(id_df_3_way, missingness, bd_not_outliers)
    
  } else if(sum(id_count$n) == 4){ # for this 3-way, sum of counts is 4 (1 + 1 + 2)
    
    id_df_n_way <- tibble(ID1 = id_count$value[id_count$n == 1][1],
                          ID2 = id_count$value[id_count$n == 1][2]) # to maximize
    # sample size, remove the relation that is common to both
    in_bd_not_outliers <- bd_not_outliers %>%
      filter(f.eid %in% as.integer(as.matrix(id_df_n_way))) %>%
      nrow() # to check if the two values are present
    if (in_bd_not_outliers < 2) {
      ids_BFZ <- id_df_3_way %>%
        as.matrix() %>%
        as.integer() %>%
        unique()
      ids_BFZ <- ids_BFZ[ids_BFZ %in% bd_not_outliers$f.eid]
      id_df_n_way <- tibble(ID1 = ids_BFZ[1],
                            ID2 = ids_BFZ[length(ids_BFZ)])
    }
    
    related.df[i, 2] <- ukb_related(id_df_n_way, missingness, bd_not_outliers)
    
  } else if(nrow(id_count) == 3 && sum(id_count$n) == 6){ # for this 3-way,
    # sum of counts is 6 (2 + 2 + 2)
    related.df[i, 2] <- ukb_related(id_df_3_way, missingness, bd_not_outliers)
    # just chose the least missing id of the 3
    
  } else { # any other relations
    in_bd_not_outliers <- bd_not_outliers %>%
      filter(f.eid %in% as.integer(as.matrix(id_df_3_way))) %>%
      select(f.eid)
    id_df_n_way <- tibble(ID1 = in_bd_not_outliers$f.eid[1],
                          ID2 = in_bd_not_outliers$f.eid[nrow(in_bd_not_outliers)])
    related.df[i, 2] <- ukb_related(id_df_n_way, missingness, bd_not_outliers)
  }
  
  message(paste(round(i / length(related_ids) * 100, 3), "% complete", sep = ""))
}

for_inclusion <- related.df %>%
  filter(is.na(include) == FALSE) %>%
  unique()
length(unique(for_inclusion$include)) # 12855 participants

for_exclusion <- related_ids[!related_ids %in% for_inclusion$include] 
length(for_exclusion) # 3065 participants

bd_unrelated <- bd_not_outliers %>%
  filter(!f.eid %in% for_exclusion) 
nrow(bd_unrelated) # 105388 participants


# Bendroflumethiazide (BFZ) and other thiazide BFZ users 
# ------------------------------------------------------
# Get medications for the HTN patients (required field is 20003)
bd_drugs <- select(bd, contains(c("f.eid", "f.53.", "f.20003."))) %>%
  filter(f.eid %in% bd_unrelated$f.eid) %>%
  ukb_reshape_long(dict) 

dict_meds <- as_tibble(fread("Codings.csv")) # codings for medications == 4
dict_meds <- dict_meds %>%
  subset(Coding == 4) %>%
  select(Value, Meaning)
colnames(dict_meds) <- c("Treatment/medication code","Treatment/medication")
dict_meds$`Treatment/medication code` <-
  parse_integer(dict_meds$`Treatment/medication code`) # Change to numeric to ensure compatibility with bd_drugs

bd_drugs <- bd_drugs %>%
  left_join(dict_meds) %>%
  select(!`Treatment/medication code`) %>%
  filter(I == ".0.") #Join data set with the dict and retain only first instance

# Taking BFZ
drug <- "bendrofluazide|bendroflumethiazide|neo-naclex|prestim|neo-bendromax"
index <- str_detect(tolower(bd_drugs$`Treatment/medication`), drug)
index[is.na(index)] <- FALSE
bd_BFZ <- bd_drugs[index, ] %>%
  select(!`Treatment/medication`) %>%
  unique()
nrow(bd_BFZ) # 24298 patients taking BFZ at recruitment

# Taking other thiazides
drug <- "accuretic tablet|amiloride hcl+cyclopenthiazide 2.5mg/250micrograms tablet|bisoprolol fumarate+hydrochlorothiazide 10mg/6.25mg tablet|capozide 50mg tablets x28|capozide tablet|captopril+hydrochlorothiazide 25mg/12.5mg tablet|carace 10 plus tablet|coaprovel 150mg/12.5mg tablet|co-diovan 80mg/12.5mg tablet|co-prenozide|cozaar-comp 50mg/12.5mg tablet|co-zidocapt 25mg/12.5mg tablet|cyclopenthiazide|diltiazem hcl+hydrochlorothiazide 150mg/12.5mg m/r capsule|dytide capsule|enalapril maleate+hydrochlorothiazide 20mg/12.5mg tablet|hydrochlorothiazide|hydroflumethiazide|innozide tablet|irbesartan+hydrochlorothiazide 150mg/12.5mg tablet|lisicostad hct 10/12.5mg tablet|lisinopril+hydrochlorothiazide 10mg/12.5mg tablet|losartan potassium+hydrochlorothiazide 50mg/12.5mg tablet|methyclothiazide|metoprolol tartrate+hydrochlorothiazide 100mg/12.5mg tablet|micardisplus 40mg/12.5mg tablet|moducren tablet|moduret 25 tablet|moduretic tablet|monozide 10 tablet|navidrex 500mcg tablet|navispare tablet|sotalol hydrochloride+hydrochlorothiazide 80mg/12.5mg tablet|telmisartan+hydrochlorothiazide 40mg/12.5mg tablet|timolol maleate+co-amilozide 10mg/2.5mg/25mg tablet|tolerzide tablet|trasidrex tablet|valsartan+hydrochlorothiazide 80mg/12.5mg tablet|zestoretic 10 tablet|amil-co tablet|atenixco 50mg/12.5mg tablet|atenolol+chlortalidone|atenolol+chlorthalidone|atenolol+co-amilozide|brinaldix k tablet|chlortalidone|chlorthalidone|co-amilozide|co-tenidone|coversyl plus 4mg/1.25mg tablet|diurexan 20mg tablet|hygroton 50mg tablet|indapamide|kalspare tablet|kalten capsule|metenix-5 tablet|metolazone|metoprolol tartrate+chlorthalidone 100mg/12.5mg tablet|natrilix sr 1.5mg m/r tablet|nindaxa 2.5mg tablet|perindopril+indapamide|tenoret 50 tablet|tenoretic tablet|viskaldix tablet|xipamide"
index <- str_detect(tolower(bd_drugs$`Treatment/medication`), drug)
index[is.na(index)] <- FALSE
bd_other_thiazides <- bd_drugs[index, ] %>%
  select(!`Treatment/medication`) %>%
  unique()
nrow(bd_other_thiazides)  # 3693 patients

# Make HTN data frame
bd_HTN <- bd_drugs %>%
  select(!`Treatment/medication`) %>%
  mutate(thiazide = "No thiazide") %>%
  unique()
bd_HTN$thiazide[bd_HTN$f.eid %in% bd_other_thiazides$f.eid] <- "Other thiazides"
bd_HTN$thiazide[bd_HTN$f.eid %in% bd_BFZ$f.eid] <- "BFZ"
bd_HTN$thiazide <- factor(bd_HTN$thiazide)

nrow(bd_HTN) # 105388 analyzed_participants
sum(bd_HTN$thiazide == "BFZ") # 24298 BFZ_participants
sum(bd_HTN$thiazide == "Other thiazides") # 3471 other_thiazides_participants
sum(bd_HTN$thiazide == "No thiazide") # 77619 no_thiazides_participants


# Covariates data frame
# ---------------------
# Obtain Sex (field 31), age (21003), potassium (30520), sodium (30530), glucose (30740) and urate (30880)
bd_long <- bd %>%
  filter(f.eid %in% bd_HTN$f.eid) %>%
  select(contains(c("f.eid", "f.31.", "f.53.", "f.21003.", "f.30520.", 
                    "f.30530.", "f.30740.", "f.30880.")
                  )
         ) %>%
  ukb_reshape_long(dict) %>%
  filter(I == ".0.")

bd_cov <- bd_long %>%
  select(f.eid, I, `Date of attending assessment centre`,
         `Age when attended assessment centre`, Sex
  ) %>%
  unique() %>%
  rename(date = `Date of attending assessment centre`,
         age = `Age when attended assessment centre`,
         sex = Sex
  )

# Combine with drug status
bd_HTN <- bd_HTN %>%
  select(f.eid, thiazide)
bd_cov <- left_join(bd_cov, bd_HTN)

# Add the principal components of genetic ancestry
pcs <- as_tibble(fread("UKBB_principal_components.csv.gz")) #load data
colnames(pcs) <- c("f.eid", paste0("C", 1:40))
pcs <- pcs[1:11] %>% filter(f.eid %in% bd_cov$f.eid)
bd_cov <- left_join(bd_cov, pcs)

# Deal with factors
factors <- c("sex", "thiazide")
levels <- list(c("Female", "Male"), c("No thiazide", "BFZ", "Other thiazides"))
bd_cov[factors] <- bd_cov %>%
  select(all_of(factors)) %>%
  map2(levels, factor) 


# Outcomes (log transform at the end)
# -----------------------------------
bd_outcomes <- bd_long %>%
  select(f.eid, I, `Date of attending assessment centre`,
         Glucose,
         Urate,
         `Potassium in urine`,
         `Sodium in urine`) %>%
  unique() %>%
  rename(date = `Date of attending assessment centre`,
         glucose = Glucose,
         urate = Urate,
         potassium = `Potassium in urine`,
         sodium = `Sodium in urine`
  )  %>%
  group_by(f.eid) %>%
  mutate_at(vars(glucose:sodium), median, na.rm = TRUE) %>%
  unique() %>%
  mutate_at(vars(glucose:sodium), log)

bd_all_df <- left_join(bd_cov, bd_outcomes) 
write_rds(bd_all_df, "bd_all_df.rds") # save for downstream work


# Obtain discovery and validation cohorts
# ---------------------------------------
bd_HTN <- read_rds("bd_all_df.rds")

# discovery and validation cohorts
bd_HTN %>% filter(thiazide == "BFZ") %>% nrow()  # 24298
bd_HTN %>% filter(thiazide == "Other thiazides") %>% nrow()  # 3471
# i.e. a ratio of 3471:24298 or 1:7 so split the control group that way 12.5:87.5
# discovery cohort control group to be 87.5% of the 'No thiazides' group.
control <- bd_HTN %>% filter(thiazide == "No thiazide")  # 77619

set.seed(7) # for reproducibility
index <- sample(1:nrow(control), nrow(control) * 0.875, replace = FALSE) # randomly sample without replacement
control_discovery <- control[index, ]  # 67916
control_validation <- control[-index, ]  # 9703

bd_HTN_dis <- bd_HTN %>% 
  filter(thiazide == "BFZ") %>%
  bind_rows(control_discovery) %>%
  mutate(thiazide = factor(thiazide))  %>%
  mutate(thiazide = fct_recode(thiazide, "No" = "No thiazide", "Yes" = "BFZ")) # 92214

bd_HTN_val <- bd_HTN %>% 
  filter(thiazide == "Other thiazides") %>%
  bind_rows(control_validation)  %>%
  mutate(thiazide = factor(thiazide)) %>%
  mutate(thiazide = fct_recode(thiazide, "No" = "No thiazide", "Yes" = "Other thiazides")) # 13174


# Descriptive analysis 
# --------------------

# Plots for variables (by cohorts and thiazides)
cohorts <- list(bd_HTN_dis, bd_HTN_val)
file_names <- c("dis", "val")
for (i in seq_along(cohorts)){
  bd_cohort <- cohorts[[i]] %>%
    select(age:thiazide, glucose:sodium) %>%
    mutate(glucose = exp(glucose),
           urate = exp(urate),
           potassium = exp(potassium),
           sodium = exp(sodium),
           thiazide = fct_recode(thiazide, "No thiazides" = "No", "Thiazides" = "Yes")
    ) 
  
  # continuous
  new_covariates <- bd_cohort %>% keep(is.numeric) %>% colnames()
  output <- list("list", length(new_covariates))
  for(j in seq_along(new_covariates)){
    output[[j]] <- ukb_box_plots(bd_cohort, new_covariates[[j]])
  }
  png(paste0(file_names[[i]],"_continuous_variables.png"), width = 1500, height = 1500, res = 120)
  n_columns <- 3
  grid.arrange(grobs = output, 
               ncol = n_columns)
  dev.off()
}
# categorical
png(paste0(file_names[[1]],"_categorical_variables.png"), width = 750, height = 750, res = 120)
ukb_bar_plots(bd_HTN_dis, "sex") # only sex
dev.off()
png(paste0(file_names[[2]],"_categorical_variables.png"), width = 750, height = 750, res = 120)
ukb_bar_plots(bd_HTN_val, "sex") # only sex
dev.off()

# Plots by cohorts
cohort_dis <- bd_HTN_dis %>%
  select(age:thiazide, glucose:sodium) %>%
  mutate (cohort = "Discovery cohort")
bd_cohort <- bd_HTN_val %>%
  select(age:thiazide, glucose:sodium) %>%
  mutate (cohort = "Validation cohort") %>%
  rbind(cohort_dis)
bd_cohort$cohort <- factor(bd_cohort$cohort)
# continuous
new_covariates <- bd_cohort %>% keep(is.numeric) %>% colnames()
output <- list("list", length(new_covariates))
for(j in seq_along(new_covariates)){
  output[[j]] <- ukb_box_plots(bd_cohort, new_covariates[[j]], "cohort")
}
png("continuous_variables.png", width = 1500, height = 1500, res = 120)
n_columns <- 3
grid.arrange(grobs = output, 
             ncol = n_columns)
dev.off()
# categorical
new_covariates <- bd_cohort %>% keep(is.factor) %>% select(-cohort) %>% colnames()
output <- list("list", length(new_covariates))
for(j in seq_along(new_covariates)){
  output[[j]] <- ukb_bar_plots(bd_cohort, new_covariates[[j]], "cohort")
}
png("categorical_variables.png", width = 1500, height = 750, res = 120)
n_columns <- 2
grid.arrange(grobs = output, 
             ncol = n_columns)
dev.off()

# Descriptive summary by cohort
cohorts <- list(bd_HTN_dis, bd_HTN_val)
file_names <- c("Discovery cohort", "Validation cohort")
for (i in seq_along(cohorts)){
  bd_cohort <- cohorts[[i]] %>%
    select(age:thiazide, glucose:sodium) %>%
    mutate(glucose = exp(glucose),
           urate = exp(urate),
           potassium = exp(potassium),
           sodium = exp(sodium),
           thiazide = fct_recode(thiazide, "No thiazides" = "No", "Thiazides" = "Yes")
    ) 
  cohort_summary <- ukb_summary(bd_cohort)
  colnames(cohort_summary) <- c("Independent variables", "statistic", file_names[[i]])
  if (i == 1) cohort_summary_all <- cohort_summary 
  else cohort_summary_all <- full_join(cohort_summary_all, cohort_summary)
}
write.csv(cohort_summary_all, "cohort_summary_all.csv", row.names = FALSE)
  
# Included vs excluded
outcomes <- c("glucose", "urate", "potassium", "sodium")
cohort_summary_outcomes <- vector("list", length(outcomes))
names(cohort_summary_outcomes) <- outcomes
nrow_included_outcomes <- vector("list", length(outcomes))
names(nrow_included_outcomes) <- outcomes
nrow_excluded_outcomes <- vector("list", length(outcomes)) 
names(nrow_excluded_outcomes) <- outcomes

for(j in seq_along(outcomes)){
  outcome <- outcomes[[j]]
  nrow_included <- vector("double", length(cohorts))
  nrow_excluded <- vector("double", length(cohorts))
  for(i in seq_along(cohorts)){
    bd_cohort <- cohorts[[i]] %>%
      select(age:thiazide, all_of(outcome))
    included <- bd_cohort[!is.na(bd_cohort[ncol(bd_cohort)]),] %>%
      select(-all_of(outcome))
    excluded <- bd_cohort[is.na(bd_cohort[ncol(bd_cohort)]),] %>%
      select(-all_of(outcome))
    nrow_included[[i]] <- nrow(included)
    nrow_excluded[[i]] <- nrow(excluded)
    included_summary <- ukb_summary(included) %>%
      rename(included = value)
    excluded_summary <- ukb_summary(excluded) %>%
      select(value) %>%
      rename(excluded = value)
    all_summary <- bind_cols(included_summary, excluded_summary) 
    colnames(all_summary) <- c("Independent variables", "statistic", paste0(file_names[[i]], "_included"),
                               paste0(file_names[[i]], "_excluded"))
    if (i == 1) cohort_summary_all <- all_summary
    else cohort_summary_all <- full_join(cohort_summary_all, all_summary)
  }
  nrow_included_outcomes[[j]] <- nrow_included
  nrow_excluded_outcomes[[j]] <- nrow_excluded
  cohort_summary_outcomes[[j]] <- cohort_summary_all
}
write.csv(cohort_summary_outcomes, "cohort_summary_outcomes.csv", row.names = FALSE)

# Quantile-quantile (QQ)  plots for normality for the outcomes 
# Non-transformed
png("QQ_outcomes.png", width = 2000, height = 1000, res = 120)
par(mfrow = c(2, 4))
for(i in seq_along(cohorts)){
  bd_cohort <- cohorts[[i]] %>%
    select(glucose:sodium) %>%
    mutate(glucose = exp(glucose),
           urate = exp(urate),
           potassium = exp(potassium),
           sodium = exp(sodium)
           ) 
  ukb_qq_plots(bd_cohort, cohort = file_names[[i]])
}
dev.off()

# Transformed
png("QQ_log_outcomes.png", width = 2000, height = 1000, res = 120)
par(mfrow = c(2, 4))
for(i in seq_along(cohorts)){
  bd_cohort <- cohorts[[i]] %>%
    select(glucose:sodium) 
  ukb_qq_plots(bd_cohort, cohort = file_names[[i]])
}
dev.off()

# Linearity assumptions 
file_names <- c("dis", "val")
continuous_covariates <- c("age")
for(i in seq_along(cohorts)){
  png(paste("Non_linearity_", file_names[[i]], ".png", sep = ""),
      width = 1500, height = 1200, res = 120)
  par(mfrow = c(2, 2))
  for(j in seq_along(outcomes)){
    outcome <- outcomes[[j]]
    bd_cohort <- cohorts[[i]]
    bd_cohort <- bd_cohort[is.na(bd_cohort[[outcome]]) == FALSE, ]
    for (k in seq_along(continuous_covariates)){
      continuous_covariate <- continuous_covariates[[k]]
      # ordinary least squares
      model_ols <- ols(bd_cohort[[outcome]] ~ bd_cohort[[continuous_covariate]])
      # restricted cubic splines, 5 knots
      model_rcs <- ols(bd_cohort[[outcome]] ~ rcs(bd_cohort[[continuous_covariate]],
                                                  quantile(bd_cohort[[continuous_covariate]],
                                                           c(0.05, 0.275, 0.5, 0.725, 0.95)
                                                           )
                                                  )
                       )
      # make plot
      plot(x = bd_cohort[[continuous_covariate]],
           y = bd_cohort[[outcome]],
           xlab = ukb_label_names(continuous_covariate),
           ylab = paste(str_to_title(outcome), " (", ukb_label_names_2(outcome),
                        "); ", "log transformed ", sep = ""),
           main = paste(ukb_title_names(outcome), ": ",
                        ukb_title_names(continuous_covariate), sep = "")
      )
      lines(bd_cohort[[continuous_covariate]],
            model_ols$fitted.values,
            col = "red", lwd = 2)
      points(bd_cohort[[continuous_covariate]],
             model_rcs$fitted.values,
             col = "blue", pch = 15, cex = 1)
    }
  }
  dev.off()
}


# Prepare files for GWAS that includes testing for interactions (to be run using PLINK)
# -------------------------------------------------------------------------------------
cohorts <- list(bd_HTN_dis, bd_HTN_val)
file_names <- c("dis", "val")
for (i in seq_along(cohorts)){
  bd_pheno <- cohorts[[i]] %>%
    select(f.eid, glucose:sodium) %>%
    mutate(FID = f.eid,
           IID = f.eid
           ) %>%
    select(-f.eid) %>%
    relocate(glucose:sodium, .after = last_col())
  write.table(bd_pheno, file = paste("pheno_BFZ_", file_names[[i]], ".txt", sep = ""), 
              sep = " ", row.names = FALSE, quote = FALSE, na = "NA")
  
  bd_keep <- bd_pheno %>% select(FID, IID)
  names(bd_keep) <- NULL
  write.table(bd_keep, file = paste("keep_BFZ_", file_names[[i]], ".txt", sep = ""), 
              sep = " ", row.names = FALSE, quote = FALSE)
  
  bd_covar <- cohorts[[i]] %>%
    select(age:C10) %>%
    mutate(FID = cohorts[[i]]$f.eid,
           IID = cohorts[[i]]$f.eid) %>%
    relocate(age:C10, .after = last_col()) %>%
    filter(FID %in% bd_pheno$FID)
  write.table(bd_covar, file = paste("covar_BFZ_", file_names[[i]], ".txt", sep = ""), 
              sep = " ", row.names = FALSE, quote = FALSE, na = "NA")
}


# GWAS including interactions (glucose as outcome, but applies to other outcomes)
# -------------------------------------------------------------------------------
# Note: Linux code is commented out.

# # On Linux platform (go to directory with all relevant files)
# cohort="dis"
# outcome="glucose"
# for i in {1..22}; do
#  plink2 \
#   --bgen ukb_imp_chr${i}_v3.bgen ref-first \
#   --sample ukbb.sample \
#   --pheno pheno_BFZ_${cohort}.txt \
#   --pheno-name ${outcome} \
#   --glm interaction \
#   --keep keep_BFZ_${cohort}.txt \
#   --covar covar_BFZ_${cohort}.txt \
#   --covar-name age,sex,thiazide,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10 
#   --parameters 1-14,17 \
#   --out ${cohort}_chr_${i} \
#   --geno 0.05 \
#   --maf 0.01 \
#   --hwe 0.000001 \
#   --minimac3-r2-filter 0.4
# done

## After running the above script
## SNP effects
# cohort="dis"
# outcome="glucose"
# echo "CHR BP SNP alleleA alleleB A1 TEST OBS_CT BETA SE T_STAT P" > UKBB_plink_${cohort}_${outcome}_Manhattan_info.txt
# for i in {1..22}; do
#   grep -E -Ew ADD ${cohort}_chr_${i}.${outcome}.glm.linear >> UKBB_plink_${cohort}_${outcome}_Manhattan_info.txt
# done
# sed '1d' UKBB_plink_${cohort}_${outcome}_Manhattan_info.txt > UKBB_plink_${cohort}_${outcome}_Manhattan.txt 
# bgzip --threads 32 UKBB_plink_${cohort}_${outcome}_Manhattan.txt

# # Interaction effects
# cohort="dis"
# outcome="glucose"
# echo "CHR BP SNP alleleA alleleB A1 TEST OBS_CT BETA SE T_STAT P" > UKBB_plink_${cohort}_${outcome}_Manhattan_int_info.txt
# for i in {1..22}; do
#   grep -E -Ew ADDxthiazideYes ${cohort}_chr_${i}.${outcome}.glm.linear >> UKBB_plink_${cohort}_${outcome}_Manhattan_int_info.txt
# done
# sed '1d' UKBB_plink_${cohort}_${outcome}_Manhattan_int_info.txt > UKBB_plink_${cohort}_${outcome}_Manhattan_int.txt 
# bgzip --threads 32 UKBB_plink_${cohort}_${outcome}_Manhattan_int.txt

# In an R environment
cohort <- "dis"
outcome <- "glucose"
results <- as_tibble(fread(paste("UKBB_plink_", cohort, "_", outcome, 
                                 "_Manhattan.txt.gz", sep = "")))  # use "_Manhattan_int.txt.gz" for interaction effects
colnames(results) <- c("CHR", "BP", "SNP", "alleleA", "alleleB", "A1", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P")
results2 <- results %>%
  arrange(P) 
write.csv(results2[results2$P < 1e-5,],  
          file = paste("UKBB_plink_",cohort, "_", outcome, ".csv", sep = ""), 
          row.names = FALSE)
nrow(results) # analyzed SNPs
nrow(results2[results2$P < 5e-8,])  # genome significant SNPs
nrow(results2[results2$P < 1e-5,])  # nominally significant SNPs
summary(results$OBS_CT) # number of analyzed participants (differs per SNP)

# Manhattan plot
ukb_manhattan_plots(results, cohort, outcome)
# Quantile-quantile (QQ) plot
ukb_qq_gwas_plots(results, cohort, outcome)
# Genomic inflation factor (lambda)
ukb_lambda(results) 

# For FUMA submission
cohort <- "dis"
outcome <- "glucose"
results <- results %>% select(rsID, A1, Beta, SE, P)
write.table(results, file = paste0("FUMA_", cohort, "_", outcome, ".txt"), row.names = FALSE, quote = FALSE)