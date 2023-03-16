## Genome-wide gene-environment interaction study for thiazide-associated electrolyte disorders in the UK Biobank.  
## ---------------------------------------------------------------------------------------------------------------

# Libraries and relevant functions
# --------------------------------
library(data.table) # reading tables faster
library(tidyverse) # data processing 
library(qqman) # for generating Manhattan and QQ plots
library(gridExtra) # for arranging multiple ggplots on the same page
library(grid)
library(rms) # assessing non-linearity using restricted cubic splines
library(genpwr) # for sample size calculations
library(mice) # for imputation
library(magick) # image annotation - not needed for github
library(rvest) # web scraping
library(bigsnpr) # analysis of genotype data
source("relevant_functions.R")


# Minimum sample size
# --------------------
# Based on one of the primary outcomes (blood glucose). Assuming:
  # 80% power
  # significance threshold of 2.5e-09 (genome-wide vQTL/QTL analysis) and 0.05 (one GEI test)
  # Standard deviation of 1.2 mmol/L (for 431 764 UK Biobank participants, https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=30740)
  # R2 of 1% in outcome and for genetic, environmental and genetic/environment interaction values
  # MAFs ranging from 5 to 20%
  # Additive mode of inheritance
  # 25% of hypertensive UK Biobank participants being prescribed thiazides

# Genome-wide vQTL and QTL analysis
ss <- ss.calc.linear(
  power = 0.8,
  MAF = c(0.05, 0.1, 0.2),
  ES = NULL,
  R2 = 0.01,
  sd_y = 1.2,
  Alpha = 2.5e-9,
  True.Model = "Additive",
  Test.Model = "Additive"
)
ceiling(range(ss["N_total_at_Alpha_2.5e-09"])) # 4 606 participants

# GEI analysis
ss <- ss_envir.calc.linear_outcome(
  pow = 0.8,
  MAF = c(0.05, 0.10, 0.20),
  ES_G = NULL,
  ES_E = NULL,
  ES_GE = NULL,
  P_e = 0.25,
  R2_G = 0.01,
  R2_E = 0.01,
  R2_GE = 0.01,
  sd_y = 1.2,
  Alpha = c(0.05, 0.005, 0.0005, 0.00005),
  True.Model = "Additive",
  Test.Model = "Additive"
)
ceiling(range(ss$N_at_Alpha_0.05)) # 766 (1 GEI test)
ceiling(range(ss$N_at_Alpha_0.005)) # 1 298 (10 GEI tests)
ceiling(range(ss["N_at_Alpha_5e-04"])) # 1 822 (100 GEI tests)
ceiling(range(ss["N_at_Alpha_5e-05"])) # 2 339 (1 000 GEI tests)


# Hypertensive participants
# -------------------------
# Import converted UK Biobank dataset
bd <- as_tibble(fread("bd_enc.csv.gz")) 
length(unique(bd$f.eid)) # 502 413 UK Biobank participants in total

withdrawn_consent <- as_tibble(fread("withdrawn_consent.csv")) # this is emailed to
# approved users or can be found accessed through the UK Biobank Access Management System.
bd <- bd %>% filter(!f.eid %in% withdrawn_consent$V1)
length(unique(bd$f.eid)) # 502 406 UK Biobank participants with consent

# Turn into long-format for diseases to get HTN at baseline (requires fields 53 and 20002) self-reported, verbal interview
# In the codings, 1065 = hypertension and 1072 = essential hypertension. 1073 (gestational hypertension/pre-eclampsia) not included
dict <- as_tibble(fread("Data_Dictionary_Showcase.csv")) # import dictionary
bd_HTN <- select(bd, contains(c("f.eid", "f.53.", "f.20002."))) %>%
  ukb_reshape_long(dict) %>%
  filter(`Non-cancer illness code, self-reported` == 1065 | `Non-cancer illness code, self-reported` == 1072) %>%
  filter(I == ".0.") %>%
  select(f.eid) %>%
  unique()
nrow(bd_HTN) # 133 271 participants with self-reported hypertension at baseline

# Self-reported, touch screen Data-Field 6150: Description:	Vascular/heart problems diagnosed by doctor
bd_HTN_ts <- select(bd, contains(c("f.eid", "f.53.", "f.6150."))) %>%
  ukb_reshape_long(dict) %>%
  filter(I == ".0." & `Vascular/heart problems diagnosed by doctor` == "High blood pressure") %>%
  select(f.eid) %>%
  unique() # 135 735 participants

# Add both interview and touch screen
bd_HTN <- rbind(bd_HTN, bd_HTN_ts) %>%
  unique()
nrow(bd_HTN) # 138 183 participants


# Participant quality control checks 
# ----------------------------------
# Requires fields 22027 (outliers for heterozygosity or missing rate), 53 (date), 31 (sex), 
# 22001 (genetic sex), 22006 (genetic ethnic grouping) and 22019 (sex chromosome aneuploidy)

bd_check <- select(bd, contains(c("f.eid", "f.53.", "f.31.", "f.22001.", "f.22006."))) %>%
  filter(f.eid %in% bd_HTN$f.eid) %>%
  ukb_reshape_long(dict) %>%
  filter(I == ".0.")
nrow(bd_check) # 138 183 participants

# Exclude those who did not cluster with Caucasians
bd_caucasian <- bd_check %>%
  filter(`Genetic ethnic grouping` == "Caucasian") %>%
  select(-`Genetic ethnic grouping`)
nrow(bd_caucasian) # 112 794 participants

# Exclude those with sex-discrepancies 
bd_sex_checked <- bd_caucasian %>%
  mutate(Sex = as.numeric(as.character(fct_recode(factor(Sex), "1" = "Male", "0" = "Female")))) %>%
  filter(Sex == `Genetic sex`) %>%
  select(f.eid) 
nrow(bd_sex_checked) # 112 697 participants

# Sex chromosome aneuploidy (was downloaded separately)
x_aneuploidy <- as_tibble(fread("x_aneuploidy.csv.gz")) %>%
  filter(is.na(`22019-0.0`) == FALSE) %>%
  select(eid)
bd_sex_checked_ane <- bd_sex_checked %>%
  filter(!f.eid %in% x_aneuploidy$eid) 
nrow(bd_sex_checked_ane) # 112 573 participants

# Exclude heterozygosity/missing rate outliers (was downloaded separately)
outliers <- as_tibble(fread("outliers.csv")) %>% 
  filter(is.na(f.22027.0.0) == FALSE) 
bd_not_outliers <- bd_sex_checked_ane %>%
  filter(!f.eid %in% outliers$f.eid) 
nrow(bd_not_outliers) # 112 388 participants 

# Relatedness 
missingness <- as_tibble(fread("missingness.csv"))
related <- as_tibble(fread("UKBB_relatedness.csv")) %>%
  filter(Kinship > 0.0884) %>%
  select(ID1, ID2)

related_ids <- unique(c(related$ID1, related$ID2))
related_ids <- related_ids[related_ids %in% bd_not_outliers$f.eid] #to save time, limit
# to only those ids present in the bd_not_outliers dataframe, n = 16 497

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
length(unique(for_inclusion$include)) # 13 249 participants

for_exclusion <- related_ids[!related_ids %in% for_inclusion$include] 
length(for_exclusion) # 3 248 participants

bd_unrelated <- bd_not_outliers %>%
  filter(!f.eid %in% for_exclusion) 
nrow(bd_unrelated) # 109 140 participants


# Bendroflumethiazide (BFZ) and other thiazide users 
# --------------------------------------------------
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
nrow(bd_BFZ) # 24 313 patients taking BFZ at recruitment

# Taking other thiazides
drug <- "accuretic tablet|amiloride hcl+cyclopenthiazide 2.5mg/250micrograms tablet|bisoprolol fumarate+hydrochlorothiazide 10mg/6.25mg tablet|capozide 50mg tablets x28|capozide tablet|captopril+hydrochlorothiazide 25mg/12.5mg tablet|carace 10 plus tablet|coaprovel 150mg/12.5mg tablet|co-diovan 80mg/12.5mg tablet|co-prenozide|cozaar-comp 50mg/12.5mg tablet|co-zidocapt 25mg/12.5mg tablet|cyclopenthiazide|diltiazem hcl+hydrochlorothiazide 150mg/12.5mg m/r capsule|dytide capsule|enalapril maleate+hydrochlorothiazide 20mg/12.5mg tablet|hydrochlorothiazide|hydroflumethiazide|innozide tablet|irbesartan+hydrochlorothiazide 150mg/12.5mg tablet|lisicostad hct 10/12.5mg tablet|lisinopril+hydrochlorothiazide 10mg/12.5mg tablet|losartan potassium+hydrochlorothiazide 50mg/12.5mg tablet|methyclothiazide|metoprolol tartrate+hydrochlorothiazide 100mg/12.5mg tablet|micardisplus 40mg/12.5mg tablet|moducren tablet|moduret 25 tablet|moduretic tablet|monozide 10 tablet|navidrex 500mcg tablet|navispare tablet|sotalol hydrochloride+hydrochlorothiazide 80mg/12.5mg tablet|telmisartan+hydrochlorothiazide 40mg/12.5mg tablet|timolol maleate+co-amilozide 10mg/2.5mg/25mg tablet|tolerzide tablet|trasidrex tablet|valsartan+hydrochlorothiazide 80mg/12.5mg tablet|zestoretic 10 tablet|amil-co tablet|atenixco 50mg/12.5mg tablet|atenolol+chlortalidone|atenolol+chlorthalidone|atenolol+co-amilozide|brinaldix k tablet|chlortalidone|chlorthalidone|co-amilozide|co-tenidone|coversyl plus 4mg/1.25mg tablet|diurexan 20mg tablet|hygroton 50mg tablet|indapamide|kalspare tablet|kalten capsule|metenix-5 tablet|metolazone|metoprolol tartrate+chlorthalidone 100mg/12.5mg tablet|natrilix sr 1.5mg m/r tablet|nindaxa 2.5mg tablet|perindopril+indapamide|tenoret 50 tablet|tenoretic tablet|viskaldix tablet|xipamide"
index <- str_detect(tolower(bd_drugs$`Treatment/medication`), drug)
index[is.na(index)] <- FALSE
bd_other_thiazides <- bd_drugs[index, ] %>%
  select(!`Treatment/medication`) %>%
  unique()
nrow(bd_other_thiazides)  # 3 700 patients

# Make HTN data frame
bd_HTN <- bd_drugs %>%
  select(!`Treatment/medication`) %>%
  mutate(thiazide = "No thiazide") %>%
  unique()
bd_HTN$thiazide[bd_HTN$f.eid %in% bd_other_thiazides$f.eid] <- "Other thiazides"
bd_HTN$thiazide[bd_HTN$f.eid %in% bd_BFZ$f.eid] <- "BFZ"
bd_HTN$thiazide <- factor(bd_HTN$thiazide)

nrow(bd_HTN) # 109 140 analyzed participants
sum(bd_HTN$thiazide == "BFZ") # 24 313 BFZ participants
sum(bd_HTN$thiazide == "Other thiazides") # 3 478 other thiazides participants
sum(bd_HTN$thiazide == "No thiazide") # 81 349 no thiazides participants


# Covariates data frame
# ---------------------
# Obtain Sex (field 31), Townsend index (189), Types of physical activity in last 4 weeks (6164), 
# bmi (21001), age (21003), smoking status (20116), alcohol status (20117), genotyping array (22000, BiLEVE vs Axiom)
# medications 
# potassium (30520), sodium (30530), glucose (30740) and urate (30880)
# primary analysis: age, gender, 10 PCs 
# Sensitivity: genotype array (BiLEVE and Axiom array), BMI, PCs 1-40, smoking, alcohol, pa, Townsend index, medications

bd_long <- bd %>%
  filter(f.eid %in% bd_HTN$f.eid) %>%
  select(contains(c("f.eid", "f.31.", "f.53.", "f.189.", "f.6164.", "f.21001.", "f.21003.", "f.20116.", 
                    "f.20117.", "f.30520.", "f.30530.", "f.30740.", "f.30880.", "f.22000.")
                  )
         ) %>%
  ukb_reshape_long(dict) %>%
  filter(I == ".0.")

# First sort out physical exercise (the most extreme of exercises chosen as an individual could engage in many)
pa_levels <- c("Strenuous sports",
               "Other exercises (eg: swimming, cycling, keep fit, bowling)",
               "Heavy DIY (eg: weeding, lawn mowing, carpentry, digging)",
               "Light DIY (eg: pruning, watering the lawn)",
               "Walking for pleasure (not as a means of transport)",
               "None of the above", "Prefer not to answer")
bd_pa <- bd_long %>%
  mutate (pa = factor(`Types of physical activity in last 4 weeks`,
                      levels = pa_levels)) %>%
  select(f.eid, I, `Date of attending assessment centre`, pa) %>%
  group_by(f.eid) %>%
  summarise(pa = pa_levels[min(as.numeric(pa))])

# Then genotyping batch
dict_array <- as_tibble(fread("Codings.csv")) %>%
  subset(Coding == 22000) %>%
  select(Value, Meaning) %>%
  rename(array = Meaning) %>%
  mutate(Value = as.integer(Value))

bd_array <- bd_long %>%
  select(f.eid, I, `Date of attending assessment centre`,
         `Genotype measurement batch`) %>%
  unique() %>%
  rename(date = `Date of attending assessment centre`,
         Value = `Genotype measurement batch`
  ) %>% 
  left_join(dict_array) %>%
  select(f.eid, array) %>%
  mutate(array = case_when(startsWith(array, "UKBiLEVE") ~ "bileve", startsWith(array, "Batch") ~ "axiom"))

# Sort out other covariates
bd_cov <- bd_long %>%
  select(f.eid, I, `Date of attending assessment centre`,
         `Age when attended assessment centre`,
         Sex, `Body mass index (BMI)`, `Smoking status`,
         `Alcohol drinker status`,
         `Townsend deprivation index at recruitment`
         ) %>%
  unique() %>%
  rename(date = `Date of attending assessment centre`,
         age = `Age when attended assessment centre`,
         sex = Sex, 
         bmi = `Body mass index (BMI)`,
         smoking = `Smoking status`,
         alcohol = `Alcohol drinker status`,
         townsend = `Townsend deprivation index at recruitment`
         )

# Combine with PA and array
bd_cov <- left_join(bd_cov, bd_pa) %>%
  left_join(bd_array)

# Add the principal components of genetic ancestry
pcs <- as_tibble(fread("UKBB_principal_components.csv.gz")) #load data
colnames(pcs) <- c("f.eid", paste0("C", 1:40))
bd_cov <- left_join(bd_cov, pcs)

## Check medication statuses
# Separate by class https://pubmed.ncbi.nlm.nih.gov/22240117/, 
# Losartan (effect size 0.81) considered separately
# ACEIs combined with non-losartan ARBs due to similar effect sizes (1.24 vs 1.29)
# Other classes are beta blockers and CCB

bd_HTN <- bd_HTN %>%
  select(f.eid, thiazide)

# Get medications again (required field is 20003)
bd_meds <- select(bd, contains(c("f.eid", "f.53.", "f.20003."))) %>%
  filter(f.eid %in% bd_HTN$f.eid) %>%
  ukb_reshape_long(dict) %>%
  rename(Date = `Date of attending assessment centre`,
         Coding = `Treatment/medication code`) %>%
  filter(I == ".0.")
dict_meds <- as_tibble(fread("medication_GWAS_dict.csv")) # Derived from PMID: 31015401
bd_meds <- left_join(bd_meds, dict_meds) %>%
  select(-Date, -Coding) # Combine bd_meds with dict_meds and removing columns that are no longer needed

# Losartan cohort ("C09CA01" for losartan only, "C09DA01" for losartan with HCT)
ATC_code <- "C09CA01|C09DA01" 
bd_meds_losartan <- bd_meds %>% 
  filter(str_detect(.$Medication_ATC_code, ATC_code))
length(unique(bd_meds_losartan$f.eid)) # 4 140

# non-losartan RAAS
ATC_code <- "(^|\\|)C09" 
bd_meds_raas <- bd_meds %>% 
  filter(str_detect(.$Medication_ATC_code, ATC_code))
length(unique(bd_meds_raas$f.eid)) # 50 855
bd_meds_raas <- bd_meds_raas %>%
  filter(!f.eid %in% bd_meds_losartan$f.eid)
length(unique(bd_meds_raas$f.eid)) # 46 715, excludes those taking losartan

# Beta blockers
ATC_code <- "(^|\\|)C07" 
bd_meds_bb <- bd_meds %>% 
  filter(str_detect(.$Medication_ATC_code, ATC_code))
length(unique(bd_meds_bb$f.eid)) # 20 630

# CCB
ATC_code <- "(^|\\|)C08" 
bd_meds_ccb <- bd_meds %>% 
  filter(str_detect(.$Medication_ATC_code, ATC_code))
length(unique(bd_meds_ccb$f.eid)) # 26 263

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

drug <- "allopurinol"
index <- str_detect(tolower(bd_meds$`Treatment/medication`), drug)
index[is.na(index)] <- FALSE
bd_meds_allo <- bd_meds[index, ] %>%
  select(f.eid) %>%
  unique()
nrow(bd_meds_allo) # 2 827 

drug <- "sulfinpyrazone"
index <- str_detect(tolower(bd_meds$`Treatment/medication`), drug)
index[is.na(index)] <- FALSE
bd_meds_sul <- bd_meds[index, ] %>%
  select(f.eid) %>%
  unique()
nrow(bd_meds_sul) # 16 - not included due to very low numbers


# Add the columns
bd_HTN <- bd_HTN %>%
  mutate(losartan = "no",
         raas = "no",
         bb = "no",
         ccb = "no",
         allo = "no")
bd_HTN$losartan[bd_HTN$f.eid %in% bd_meds_losartan$f.eid] <- "yes"
bd_HTN$raas[bd_HTN$f.eid %in% bd_meds_raas$f.eid] <- "yes"
bd_HTN$bb[bd_HTN$f.eid %in% bd_meds_bb$f.eid] <- "yes"
bd_HTN$ccb[bd_HTN$f.eid %in% bd_meds_ccb$f.eid] <- "yes"
bd_HTN$allo[bd_HTN$f.eid %in% bd_meds_allo$f.eid] <- "yes"


# Combine with drug status
bd_cov <- left_join(bd_cov, bd_HTN)


# Deal with factors
bd_cov <- bd_cov %>%
  mutate(pa = fct_collapse(pa,
                           strenuous = c("Strenuous sports"),
                           moderate = c("Other exercises (eg: swimming, cycling, keep fit, bowling)", 
                                        "Heavy DIY (eg: weeding, lawn mowing, carpentry, digging)"),
                           mild = c("Light DIY (eg: pruning, watering the lawn)",
                                    "Walking for pleasure (not as a means of transport)"),
                           none = c("None of the above")
                           )
         )
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8682505/ coded the above as:
# NA: "None of the above", "Prefer not to answer" ("Prefer not to answer" left as NA in our analysis)
# 1: "Walking for pleasure (not as a means of transport)", "Light DIY (eg: pruning, watering the lawn)"
# 2: "Other exercises (eg: swimming, cycling, keep fit, bowling)", "Heavy DIY (eg: weeding, lawn mowing, carpentry, digging)"
# 3: "Strenuous sports"

factors <- c("sex", "smoking", "alcohol", "pa", "thiazide", "array", "losartan", "raas", "bb", "ccb", "allo")
levels <- list(c("Female", "Male"),
               c("Never", "Previous", "Current", "Prefer not to answer"),
               c("Never", "Previous", "Current", "Prefer not to answer"),
               c("none", "mild", "moderate", "strenuous", "Prefer not to answer"),
               c("No thiazide", "BFZ", "Other thiazides"),
               c("axiom", "bileve"),
               c("no", "yes"), c("no", "yes"), c("no", "yes"), c("no", "yes"), c("no", "yes")
               )
bd_cov[factors] <- bd_cov %>%
  select(all_of(factors)) %>%
  map2(levels, factor) %>%
  ukb_recode_factor_na


# Outcomes 
# --------
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
  unique() 

bd_all_df <- left_join(bd_cov, bd_outcomes) 
# BMI (484 NA's), Smoking (488 NA's), Alcohol (124 NA's), Townsend (137 NA's), PA (824 NA's) imputed.

tempData <- bd_all_df  %>%
  select(-f.eid, -I, -date) %>% # Not used for imputation - outcomes left in
  mice(m = 1, maxit = 10, seed = 7) # Impute missing data

bd_imp <- tempData %>%
  complete(1, include = FALSE) %>%
  as_tibble() %>%
  select(-glucose, -urate, -potassium, -sodium) %>%
  mutate(f.eid = bd_all_df$f.eid,
         glucose = bd_all_df$glucose,
         urate = bd_all_df$urate,
         potassium = bd_all_df$potassium,
         sodium = bd_all_df$sodium
         )
write_rds(bd_imp, "bd_all_df_imp.rds") # Save for downstream work


# Obtain discovery and validation cohorts
# ---------------------------------------
bd_HTN <- read_rds("bd_all_df_imp.rds")

# Discovery (or first) and validation (or second) cohorts
bd_HTN %>% filter(thiazide == "BFZ") %>% nrow()  # 24 313
bd_HTN %>% filter(thiazide == "Other thiazides") %>% nrow()  # 3 478
# i.e. a ratio of 3478:24313 or 1:7 so split the control group that way 12.5:87.5
# discovery cohort control group to be 87.5% of the 'No thiazides' group.
control <- bd_HTN %>% filter(thiazide == "No thiazide")  # 81 349

set.seed(7) # for reproducibility
index <- sample(1:nrow(control), nrow(control) * 0.875, replace = FALSE) # randomly sample without replacement
control_discovery <- control[index, ]  # 71 180
control_validation <- control[-index, ]  # 10 169

bd_HTN_dis <- bd_HTN %>% 
  filter(thiazide == "BFZ") %>%
  bind_rows(control_discovery) %>%
  mutate(thiazide = factor(thiazide))  %>%
  mutate(thiazide = fct_recode(thiazide, "No" = "No thiazide", "Yes" = "BFZ")) # 95 493

bd_HTN_val <- bd_HTN %>% 
  filter(thiazide == "Other thiazides") %>%
  bind_rows(control_validation)  %>%
  mutate(thiazide = factor(thiazide)) %>%
  mutate(thiazide = fct_recode(thiazide, "No" = "No thiazide", "Yes" = "Other thiazides")) # 13 674


# Descriptive analysis 
# --------------------
# Figures S1 and S2 (participant characteristics)
cohorts <- list(bd_HTN_dis, bd_HTN_val)
file_names <- c("dis", "val")
for (i in seq_along(cohorts)){
  bd_cohort <- cohorts[[i]] %>%
    select(age:array, thiazide:allo, glucose:sodium) %>%
    mutate(thiazide = fct_recode(thiazide, "No Thiazide" = "No", "Thiazide" = "Yes")) 
  
  # continuous
  new_covariates <- bd_cohort %>% keep(is.numeric) %>% colnames()
  output <- list("list", length(new_covariates))
  for(j in seq_along(new_covariates)){
    output[[j]] <- ukb_box_plots(bd_cohort, new_covariates[[j]])
  }
  png(paste0(file_names[[i]],"_continuous_variables.png"), width = 2000, height = 1000, res = 120)
  n_columns <- 4
  grid.arrange(grobs = output, 
               ncol = n_columns)
  dev.off()
  
  # categorical
  new_covariates <- bd_cohort %>% keep(is.factor) %>% colnames()
  new_covariates <- new_covariates[new_covariates != "thiazide"]
  output <- list("list", length(new_covariates))
  for(j in seq_along(new_covariates)){
    output[[j]] <- ukb_bar_plots(bd_cohort, new_covariates[[j]])
  }
  png(paste0(file_names[[i]],"_categorical_variables.png"), width = 2000, height = 1000, res = 120)
  n_columns <- 4
  grid.arrange(grobs = output, 
               ncol = n_columns)
  dev.off()
}

image_1 <- image_read("dis_continuous_variables.png") %>% 
  image_trim() %>%
  image_border(color = "white", "10x50") %>%
  ukb_image_annotate("A. First cohort (N = 95,493)", 0, 0, style = "normal", bold = TRUE, size = 30) %>%
  image_trim() %>%
  image_border(color = "white") 
image_2 <- image_read("val_continuous_variables.png") %>% 
  image_trim() %>%
  image_border(color = "white", "10x50") %>%
  ukb_image_annotate("B. Second cohort (N = 13,647)", 0, 0, style = "normal", bold = TRUE, size = 30) %>%
  image_trim() %>%
  image_border(color = "white") 
image_append(c(image_1, image_2), stack = TRUE) %>%
  image_write(path = "Fig_S1.png", format = "png")

image_1 <- image_read("dis_categorical_variables.png") %>% 
  image_trim() %>%
  image_border(color = "white", "10x50") %>%
  ukb_image_annotate("A. First cohort (N = 95,493)", 0, 0, style = "normal", bold = TRUE, size = 30) %>%
  image_trim() %>%
  image_border(color = "white") 
image_2 <- image_read("val_categorical_variables.png") %>% 
  image_trim() %>%
  image_border(color = "white", "10x50") %>%
  ukb_image_annotate("B. Second cohort (N = 13,647)", 0, 0, style = "normal", bold = TRUE, size = 30) %>%
  image_trim() %>%
  image_border(color = "white") 
image_append(c(image_1, image_2), stack = TRUE) %>%
  image_write(path = "Fig_S2.png", format = "png")


# Figure S3 (linearity assumptions) 
file_names <- c("dis", "val")
continuous_covariates <- c("age")
for(i in seq_along(cohorts)){
  png(paste("Non_linearity_", file_names[[i]], ".png", sep = ""),
      width = 1500, height = 1500, res = 120)
  par(mfrow = c(4, 3))
  outcomes <- c("glucose", "urate", "potassium", "sodium")
  for(j in seq_along(outcomes)){
    outcome <- outcomes[[j]]
    bd_cohort <- cohorts[[i]] %>%
      select(age:array, thiazide:allo, glucose:sodium)
    bd_cohort <- bd_cohort[is.na(bd_cohort[[outcome]]) == FALSE, ]
    continuous_covariates <- bd_cohort %>% keep(is.numeric) %>% colnames()
    continuous_covariates <- continuous_covariates[!continuous_covariates %in% outcomes]
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
                        ") ", sep = ""),
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

image_read("Non_linearity_dis.png") %>% 
  image_trim() %>%
  image_border(color = "white", "10x50") %>%
  ukb_image_annotate("A. First cohort (N = 95,493)", 0, 0, style = "normal", bold = TRUE, size = 30) %>%
  image_trim() %>%
  image_border(color = "white") %>%
  image_write(path = "Fig_S3A.png", format = "png")
image_read("Non_linearity_val.png") %>% 
  image_trim() %>%
  image_border(color = "white", "10x50") %>%
  ukb_image_annotate("B. Second cohort (N = 13,647)", 0, 0, style = "normal", bold = TRUE, size = 30) %>%
  image_trim() %>%
  image_border(color = "white") %>%
  image_write(path = "Fig_S3B.png", format = "png")


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
    select(age:allo) %>%
    mutate(FID = cohorts[[i]]$f.eid,
           IID = cohorts[[i]]$f.eid) %>%
    relocate(age:allo, .after = last_col()) %>%
    filter(FID %in% bd_pheno$FID)
  write.table(bd_covar, file = paste("covar_BFZ_", file_names[[i]], ".txt", sep = ""), 
              sep = " ", row.names = FALSE, quote = FALSE, na = "NA")
  
  # pre-adjusted phenotypes (main analysis)
  bd_pheno <- cohorts[[i]] %>%
    select(f.eid, glucose:sodium) %>%
    mutate(FID = f.eid,
           IID = f.eid,
           glucose = ukb_model_outcome(cohorts[[i]], "glucose"),
           urate = ukb_model_outcome(cohorts[[i]], "urate"),
           potassium = ukb_model_outcome(cohorts[[i]], "potassium"),
           sodium = ukb_model_outcome(cohorts[[i]], "sodium")
    ) %>%
    select(-f.eid) %>%
    relocate(glucose:sodium, .after = last_col())
  
  # # pre-adjusted phenotypes (sensitivity analysis)
  # bd_pheno <- cohorts[[i]] %>%
  #   select(f.eid, glucose:sodium) %>%
  #   mutate(FID = f.eid,
  #          IID = f.eid,
  #          glucose = ukb_model_outcome2(cohorts[[i]], "glucose"),
  #          urate = ukb_model_outcome2(cohorts[[i]], "urate"),
  #          potassium = ukb_model_outcome2(cohorts[[i]], "potassium"),
  #          sodium = ukb_model_outcome2(cohorts[[i]], "sodium")
  #   ) %>%
  #   select(-f.eid) %>%
  #   relocate(glucose:sodium, .after = last_col())
  
  outcomes <- c("glucose", "urate", "potassium", "sodium")
  for (j in seq_along(outcomes)) {
    write.table(bd_pheno[, c("FID", "IID", outcomes[[j]])], file = paste(outcomes[[j]], "_BFZ_", file_names[[i]], "_osca.txt", sep = ""), 
                sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, na = "NA")
  }
}


# Table S6 (SNPs and PCs excluded)
# --------------------------------
model <- "urate~age+sex+bmi+array+smoking+alcohol+pa+townsend+thiazide+losartan+raas+bb+ccb+allo" %>%
  as.formula %>%
  lm(., data = bd_HTN_dis) %>%
  summary
model


# Genome-wide vQTL, QTL and GEIs (glucose as outcome, but applies to other outcomes)
# ---------------------------------------------------------------------------------
# Note: Linux code is commented out.

# # On Linux platform (go to directory with all relevant files)

## Obtain genotype data in PLINK format
# cohort="dis"
# outcome="glucose"
#
# for i in {1..22}; do
#  plink2 \
#   --bgen ukb_imp_chr${i}_v3.bgen ref-first \
#   --sample ukbb.sample \
#   --keep keep_BFZ_${cohort}.txt \
#   --make-bed \
#   --out ${cohort}_chr_${i} \
#   --geno 0.05 \
#   --maf 0.05 \
#   --hwe 0.00001 \
#   --minimac3-r2-filter 0.3 \
#   --hard-call-threshold 0.1 \
#   --rm-dup 'force-first'
#
## Conduct vQTL analysis (assumes osca-0.46.1 has been installed from https://yanglab.westlake.edu.cn/software/osca/#Download)
#  osca-0.46.1 \
#   --vqtl \
#   --bfile ${cohort}_chr_${i} \
#   --pheno ${outcome}_BFZ_${cohort}_osca.txt \
#   --vqtl-mtd 2 \
#   --out osca_${cohort}_chr_${i}.${outcome}
#
## Conduct QTL analysis
#
#  plink2 \
#   --bfile ${cohort}_chr_${i} \
#   --pheno ${outcome}_BFZ_${cohort}_osca.txt \
#   --glm \
#   --out plink_${cohort}_chr_${i}
# done
#
## After running the above script
## vQTL
# sed -n 1p osca_${cohort}_chr_1.${outcome}.vqtl > osca_${cohort}_${outcome}.txt
# for i in {1..22}; do
# sed "1d" osca_${cohort}_chr_${i}.${outcome}.vqtl >> osca_${cohort}_${outcome}.txt
# done
# bgzip --threads 32 osca_${cohort}_${outcome}.txt
#
## QTL
# awk 'NR==1' plink_${cohort}_chr_1.PHENO1.glm.linear > plink_${cohort}_${outcome}_Manhattan.txt
# for i in {1..22}; do
# grep -E -Ew ADD plink_${cohort}_chr_${i}.PHENO1.glm.linear >> plink_${cohort}_${outcome}_Manhattan.txt
# done
# bgzip --threads 32 plink_${cohort}_${outcome}_Manhattan.txt
#
## In an R environment
## vQTL
cohort <- "dis"
outcome <- "glucose"
p_threshold <- 1e-8/4
type <- "vqtl"
results <- as_tibble(fread(paste("osca_", cohort, "_", outcome, ".txt.gz", sep = "")))  %>%
  rename(CHR = Chr,
         MAF = freq,
         BP = bp)
write.csv(results[results$P < p_threshold,],  
          file = paste("osca_",cohort, "_", outcome, ".csv", sep = ""), 
          row.names = FALSE)
nrow(results) # Number of analyzed SNPs
nrow(results[results$P < p_threshold,]) # Number of significant vQTL SNPs
summary(results$NMISS) # Number of analyzed participants
# Manhattan plot
ukb_manhattan_plots(results, cohort, outcome, type, genomewideP = p_threshold, suggestive = FALSE)
# Quantile-quantile (QQ) plot
ukb_qq_gwas_plots(results, cohort, outcome, type)
# Genomic inflation factor (lambda)
ukb_lambda(results) 

## QTL
type <- "qtl"
results <- as_tibble(fread(paste("plink_", cohort, "_", outcome, 
                                 "_Manhattan.txt.gz", sep = "")))  
colnames(results) <- c("CHR", "BP", "SNP", "alleleA", "alleleB", "A1", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P")
write.csv(results[results$P < p_threshold,],  
          file = paste("plink_",cohort, "_", outcome, ".csv", sep = ""), 
          row.names = FALSE)
nrow(results)
nrow(results[results$P < p_threshold,])  
summary(results$OBS_CT) 
# Manhattan plot
ukb_manhattan_plots(results, cohort, outcome, type, genomewideP = p_threshold, suggestive = FALSE)
# Quantile-quantile (QQ) plot
ukb_qq_gwas_plots(results, cohort, outcome, type)
# Genomic inflation factor (lambda)
ukb_lambda(results)


## Downstream analysis below is for vQTLs but the same pipeline applies to QTLs
# Get top SNPs
top_snps <- fread(paste0("osca_",cohort, "_", outcome, ".csv")) %>% 
  select(SNP)
write.table(top_snps, paste0(cohort, "_top_snps_BFZ_osca.txt"), sep = " ", row.names = FALSE, 
            col.names = FALSE, quote = FALSE)

## Get PLINK files for the top SNPs
# cohort="dis"
# for i in {1..22}; do
# plink2 \
#  --bgen ukb_imp_chr${i}_v3.bgen ref-first \
#  --sample ukbb.sample \
#  --keep keep_BFZ_${cohort}.txt \
#  --extract ${cohort}_top_snps_BFZ_osca.txt \
#  --make-bed 
#  --out BFZ_${cohort}_chr_${i}
# done
#
## Merge the files (used plink1.9)
# /plink1.9/plink \
#  --bfile BFZ_${cohort}_chr_1 \
#  --merge-list allfiles_BFZ.txt \
#  --make-bed \
#  --out BFZ_${cohort} 
## Note: 1. allfiles_BFZ.txt contains the files to be merged.
##       2. If SNPs producing an error (missnps) exist, exclude them before merging again.

## In R, get top SNPs file for LD analysis
outcome <- "glucose"
p_threshold <- 1e-8/4
results <- as_tibble(fread(paste0("osca_", cohort, "_", outcome, ".txt.gz")))  %>%
  rename(CHR = Chr,
         MAF = freq,
         BP = bp)
write.table(results[results$P < p_threshold,], 
            file = paste0("osca_",cohort, "_", outcome, "_ld.txt"), 
            sep = " ", 
            row.names = FALSE, 
            quote = FALSE)

## Get near-independent SNPs
# /plink1.9/plink \
#  --bfile BFZ_${cohort} \
#  --clump-p1 2.5e-9 \
#  --clump-p2 2.5e-9 \
#  --clump-r2 0.01 \
#  --clump-kb 5000 \
#  --clump osca_${cohort}_${outcome}_ld.txt \
#  --clump-snp-field SNP \
#  --clump-field P \
#  --out osca${cohort}_${outcome}_vqtl
# 
# awk 'NR!=1{print $3}' osca${cohort}_${outcome}_vqtl.clumped >  osca${cohort}_${outcome}_vqtl.snps


## GEIs for the top near-independent vQTL SNPs (done in R)
# --------------------------------------------------------
outcome <- "glucose"
cohort <- "dis" # similar analysis done for the "val" or second cohort
phenotype <- fread(paste0(outcome, "_BFZ_", cohort, "_osca.txt"))
colnames(phenotype) <- c("FID", "IID", outcome)
cov <- fread(paste0("covar_BFZ_", cohort, ".txt"))
pheno <- left_join(phenotype, cov) 
# Preprocess the bed file (only need to do once for each data set)
snp_readBed(paste0("BFZ_", cohort, ".bed"))
# Attach the genotype object
obj.bigSNP <- snp_attach(paste0("BFZ_", cohort, ".rds"))
# Assign the genotype to a variable 
genotype <- obj.bigSNP$genotypes
# Assume the fam order is the same across different chromosomes
fam.order <- as.data.table(obj.bigSNP$fam)
# Rename fam order
setnames(fam.order,
         c("family.ID", "sample.ID"),
         c("FID", "IID"))
# Reformat the phenotype file such that y is of the same order as the sample ordering in the genotype file
y <- fam.order[, 1]
y <- left_join(y, pheno)
# Get list of snps
snps_list <- fread(paste0("osca_", cohort, "_", outcome, "_ld.txt"), header = T)
snps_list$int_p <- NA
snps_list$n_thiazide <- NA
snps_list$beta_thiazide <- NA
snps_list$se_thiazide <- NA
snps_list$n_not_thiazide <- NA
snps_list$beta_not_thiazide <- NA
snps_list$se_not_thiazide <- NA
snps_ind <- fread(paste0("osca", cohort, "_", outcome, "_vqtl.snps"), header = F)
snps_list <- snps_list[snps_list$SNP %in% snps_ind$V1]
snps_ind <- fread(paste0("osca", cohort, "_", outcome, "_vqtl.snps"), header = F)
snps_list <- snps_list[snps_list$SNP %in% snps_ind$V1]
# Get the model for the SNPs (uses pre-adjusted phenotype)
for (i in seq_along(snps_list$SNP)) {
  snp_position <- which((obj.bigSNP$map)$marker.ID == snps_list$SNP[i])
  y2 <- y %>% mutate(snp = genotype[, snp_position])
  model_snp <- as.formula(paste0(outcome, " ~ snp + thiazide + snp:thiazide")) %>%
    lm(., data = y2) %>%
    summary
  snps_list$int_p[i] <- model_snp$coefficients[16]
  # Effects for SNP in each thiazide group
  snps_list$n_thiazide[i] <- nrow(y2 %>% filter(thiazide == "Yes" & is.na(snp) == FALSE & is.na(glucose) == FALSE))
  snps_list$n_not_thiazide[i] <- nrow(y2 %>% filter(thiazide == "No" & is.na(snp) == FALSE & is.na(glucose) == FALSE))
  model_snp <- as.formula(paste0(outcome, " ~ snp")) %>%
   lm(., data = y2 %>% filter(thiazide == "Yes")) %>%
  summary
  snps_list$beta_thiazide[i] <- model_snp$coefficients[2]
  snps_list$se_thiazide[i] <- model_snp$coefficients[4]
  model_snp <- as.formula(paste0(outcome, " ~ snp")) %>%
   lm(., data = y2 %>% filter(thiazide == "No")) %>%
  summary
  snps_list$beta_not_thiazide[i] <- model_snp$coefficients[2]
  snps_list$se_not_thiazide[i] <- model_snp$coefficients[4]
  }
snps_list # Contains the GEI P value and effect sizes in thiazide- vs non-thiazide-treated participants


## Code for some Figures and Tables (done after conducting genome-wide vQTLs and QTLs for all outcomes)
# ----------------------------------------------------------------------------------------------------
# Figure 2 (annotations done in word, but can be done e.g. using ukb_image_annotate())
outcomes <- c("glucose", "urate")
outcome_names <- c("Blood glucose", "Serum urate")
plot <- "manhattan"
type <- "vqtl"
for (i in seq_along(outcomes)) {
  img <- image_read(paste0("UKBB_dis_", outcomes[[i]], "_", plot, "_", type, ".png")) %>% 
    image_trim() %>%
    image_border(color = "white", "10x50") %>%
    ukb_image_annotate(paste0(LETTERS[i], ". ", outcome_names[i]), 0, 0, style = "normal", bold = TRUE, size = 30) %>%
    image_trim() %>%
    image_border(color = "white") 
  if (i == 1) img_pooled <- img else img_pooled <- image_append(c(img_pooled, img), stack = TRUE)
}
img_pooled %>%
  image_write(path = "Fig_2.png", format = "png")


# Figure S4 (vQTLs for urine potassium and urine sodium, Manhattan plots)
outcomes <- c("potassium", "sodium")
outcome_names <- c("Urine potassium", "Urine sodium")
plot <- "manhattan"
type <- "vqtl"
for (i in seq_along(outcomes)) {
  img <- image_read(paste0("UKBB_dis_", outcomes[[i]], "_", plot, "_", type, ".png")) %>% 
    image_trim() %>%
    image_border(color = "white", "10x50") %>%
    ukb_image_annotate(paste0(LETTERS[i], ". ", outcome_names[i]), 0, 0, style = "normal", bold = TRUE, size = 30) %>%
    image_trim() %>%
    image_border(color = "white") 
  if (i == 1) img_pooled <- img else img_pooled <- image_append(c(img_pooled, img), stack = TRUE)
}
img_pooled %>%
  image_write(path = "Fig_S4.png", format = "png")

# Figure S5 (vQTLs for all outcomes, QQ plots)
outcomes <- c("glucose", "urate", "potassium", "sodium")
outcome_names <- c("Blood glucose", "Serum urate", "Urine potassium", "Urine sodium")
lambdas <- round(c(1.061777, 1.035757, 1.033079, 1.049142), 3)
plot <- "qqplot"
type <- "vqtl"
for (i in seq_along(outcomes)) {
  img <- image_read(paste0("UKBB_", cohort, "_", outcomes[[i]], "_", plot, "_", type, ".png")) %>% 
    image_trim() %>%
    image_border(color = "white", "10x50") %>%
    ukb_image_annotate(paste0(LETTERS[i], ". ", outcome_names[[i]], " (λ = ", lambdas[[i]], ")"), 0, 0, style = "normal", bold = TRUE, size = 40) %>%
    image_trim() %>%
    image_border(color = "white") 
  if (i == 1) img_pooled <- img else {
    if (i == 2) img_pooled <- image_append(c(img_pooled, img)) else {
      if (i == 3) {
        img_1_2 <- img_pooled
        img_pooled <- image_append(c(img_pooled, img), stack = TRUE)
        img_3 <- img
      } else { # only 4 outcomes for now, so no need for a more general solution
        img_3_4 <- image_append(c(img_3, img))
        img_pooled <- image_append(c(img_1_2, img_3_4), stack = TRUE)
      }
    }
}
}
img_pooled %>%
  image_write(path = "Fig_S5.png", format = "png")

# Figure S6 (GWAS for urine potassium and urine sodium, Manhattan plots)
outcomes <- c("glucose", "urate", "potassium", "sodium")
outcome_names <- c("Blood glucose", "Serum urate", "Urine potassium", "Urine sodium")
plot <- "manhattan"
type <- "qtl"
for (i in seq_along(outcomes)) {
  img <- image_read(paste0("UKBB_dis_", outcomes[[i]], "_", plot, "_", type, ".png")) %>% 
    image_trim() %>%
    image_border(color = "white", "10x50") %>%
    ukb_image_annotate(paste0(LETTERS[i], ". ", outcome_names[i]), 0, 0, style = "normal", bold = TRUE, size = 30) %>%
    image_trim() %>%
    image_border(color = "white") 
  if (i == 1) img_pooled <- img else img_pooled <- image_append(c(img_pooled, img), stack = TRUE)
}
img_pooled %>%
  image_write(path = "Fig_S6.png", format = "png")

# Figure S7 (GWAS for all outcomes, QQ plots)
outcomes <- c("glucose", "urate", "potassium", "sodium")
outcome_names <- c("Blood glucose", "Serum urate", "Urine potassium", "Urine sodium")
lambdas <- round(c(1.106262, 1.166802, 1.088766, 1.118625), 3)
plot <- "qqplot"
type <- "qtl"
for (i in seq_along(outcomes)) {
  img <- image_read(paste0("UKBB_dis_", outcomes[[i]], "_", plot, "_", type, ".png")) %>% 
    image_trim() %>%
    image_border(color = "white", "10x50") %>%
    ukb_image_annotate(paste0(LETTERS[i], ". ", outcome_names[[i]], " (λ = ", lambdas[[i]], ")"), 0, 0, style = "normal", bold = TRUE, size = 40) %>%
    image_trim() %>%
    image_border(color = "white") 
  if (i == 1) img_pooled <- img else {
    if (i == 2) img_pooled <- image_append(c(img_pooled, img)) else {
      if (i == 3) {
        img_1_2 <- img_pooled
        img_pooled <- image_append(c(img_pooled, img), stack = TRUE)
        img_3 <- img
      } else { # only 4 outcomes for now, so no need for a more general solution
        img_3_4 <- image_append(c(img_3, img))
        img_pooled <- image_append(c(img_1_2, img_3_4), stack = TRUE)
      }
    }
  }
}
img_pooled %>%
  image_write(path = "Fig_S7.png", format = "png")

# Figure S8 (Sensitivity analysis for serum urate, Manhattan plots)
image_1 <- image_read("UKBB_dis_urate_manhattan_vqtl_sen.png") %>% 
  image_trim() %>%
  image_border(color = "white", "10x50") %>%
  ukb_image_annotate("A. vQTL analysis", 0, 0, style = "normal", bold = TRUE, size = 30) %>%
  image_trim() %>%
  image_border(color = "white") 
image_2 <- image_read("UKBB_dis_urate_manhattan_gwas_sen.png") %>% 
  image_trim() %>%
  image_border(color = "white", "10x50") %>%
  ukb_image_annotate("B. GWAS analysis", 0, 0, style = "normal", bold = TRUE, size = 30) %>%
  image_trim() %>%
  image_border(color = "white") 
image_append(c(image_1, image_2), stack = TRUE) %>%
  image_write(path = "Fig_S8.png", format = "png")


# Table S4 (Top vQTLs)
cohort <- "dis"
outcomes <- c("glucose", "urate", "potassium", "sodium")
n_round <- 4
for (i in seq_along(outcomes)) {
  file <- read_csv(paste0("osca_", cohort, "_", outcomes[[i]], ".csv"), show_col_types = FALSE)
  if (nrow(file) > 0 ) {
    file_dbsnp <- ukb_dbsnp(file$SNP)
    file <- left_join(file, file_dbsnp)
    file <- file %>% 
      arrange(P) %>%
      mutate(Outcome = str_to_title(outcomes[[i]]),
             `#` = c(1:nrow(file)),
             `SNP (reference/ alternative alleles)` = paste0(SNP, " (", A1, "/", A2, ")"),
             `Chromosome (position)a` = paste0(CHR, " (", BP, ")"),
             `Gene (functional consequence)b` = gsub("NA \\(NA\\)", "NA", paste0(Gene, " (", `Functional consequence`, ")")),
             N = NMISS,
             MAF = str_pad(round(MAF, n_round), 
                           n_round + 1 + str_count(MAF, "[-.]"), 
                           "right", pad = "0"),
             `Effect sizec (standard error)` = paste0(str_pad(round(beta, n_round), 
                                                              n_round + 1 + str_count(beta, "[-.]"), 
                                                              "right", pad = "0"), " (",
                                                      str_pad(round(se, n_round), 
                                                              n_round + 1 + str_count(se, "[-.]"), 
                                                              "right", pad = "0"), ")"),
             `P-value` = format(P, scientific = TRUE)
      ) %>%
      select(Outcome:N, MAF, `Effect sizec (standard error)`, `P-value`)
  } else {
    file <- tibble()
  }
  if (i == 1) all_files <- file else all_files <- rbind(all_files, file)
}
write_csv(all_files, "Table_S4.csv")

# Table S5 (Top QTLs)
for (i in seq_along(outcomes)) {
  file <- read_csv(paste0("osca_plink_", cohort, "_", outcomes[[i]], ".csv"), show_col_types = FALSE)
  if (nrow(file) > 0 ) {
    file_dbsnp <- ukb_dbsnp(file$SNP)
    file <- left_join(file, file_dbsnp)
    file <- file %>% 
      arrange(P) %>%
      mutate(Outcome = str_to_title(outcomes[[i]]),
             `#` = c(1:nrow(file)),
             `SNP (reference/ alternative alleles)` = paste0(SNP, " (", A1, "/", ifelse(alleleA == A1, alleleB, alleleA), ")"),
             `Chromosome (position)a` = paste0(CHR, " (", BP, ")"),
             `Gene (functional consequence)b` = gsub("NA \\(NA\\)", "NA", paste0(Gene, " (", `Functional consequence`, ")")),
             N = OBS_CT,
             `Effect sizec (standard error)` = paste0(str_pad(round(BETA, n_round), 
                                                              n_round + 1 + str_count(BETA, "[-.]"), 
                                                              "right", pad = "0"), " (",
                                                      str_pad(round(SE, n_round), 
                                                              n_round + 1 + str_count(SE, "[-.]"), 
                                                              "right", pad = "0"), ")"),
             `P-value` = format(P, scientific = TRUE)
      ) %>%
      select(Outcome:`P-value`)
  } else {
    file <- tibble()
  }
  if (i == 1) all_files <- file else all_files <- rbind(all_files, file)
}
write_csv(all_files, "Table_S5.csv")


# Tables S7 and S8 - Sensitivity analysis
outcomes <- c("urate")
for (i in seq_along(outcomes)) {
  file <- read_csv(paste0("osca_", cohort, "_", outcomes[[i]], "_sen.csv"), show_col_types = FALSE)
  if (nrow(file) > 0 ) {
    file_dbsnp <- ukb_dbsnp(file$SNP)
    file <- left_join(file, file_dbsnp)
    file <- file %>% 
      arrange(P) %>%
      mutate(Outcome = str_to_title(outcomes[[i]]),
             `#` = c(1:nrow(file)),
             `SNP (reference/ alternative alleles)` = paste0(SNP, " (", A1, "/", A2, ")"),
             `Chromosome (position)a` = paste0(CHR, " (", BP, ")"),
             `Gene (functional consequence)b` = gsub("NA \\(NA\\)", "NA", paste0(Gene, " (", `Functional consequence`, ")")),
             N = NMISS,
             MAF = str_pad(round(MAF, n_round), 
                           n_round + 1 + str_count(MAF, "[-.]"), 
                           "right", pad = "0"),
             `Effect sizec (standard error)` = paste0(str_pad(round(beta, n_round), 
                                                              n_round + 1 + str_count(beta, "[-.]"), 
                                                              "right", pad = "0"), " (",
                                                      str_pad(round(se, n_round), 
                                                              n_round + 1 + str_count(se, "[-.]"), 
                                                              "right", pad = "0"), ")"),
             `P-value` = format(P, scientific = TRUE)
      ) %>%
      select(Outcome:N, MAF, `Effect sizec (standard error)`, `P-value`)
  } else {
    file <- tibble()
  }
  if (i == 1) all_files <- file else all_files <- rbind(all_files, file)
}
all_files <- all_files %>%
  mutate(MAF = paste0(MAF, "$$$")) # To maintain formatting of zeros during saving
write_csv(all_files, "Table_S7.csv")
for (i in seq_along(outcomes)) {
  file <- read_csv(paste0("osca_plink_", cohort, "_", outcomes[[i]], "_sen.csv"), show_col_types = FALSE)
  if (nrow(file) > 0 ) {
    file_dbsnp <- ukb_dbsnp(file$SNP)
    file <- left_join(file, file_dbsnp)
    file <- file %>% 
      arrange(P) %>%
      mutate(Outcome = str_to_title(outcomes[[i]]),
             `#` = c(1:nrow(file)),
             `SNP (reference/ alternative alleles)` = paste0(SNP, " (", A1, "/", ifelse(alleleA == A1, alleleB, alleleA), ")"),
             `Chromosome (position)a` = paste0(CHR, " (", BP, ")"),
             `Gene (functional consequence)b` = gsub("NA \\(NA\\)", "NA", paste0(Gene, " (", `Functional consequence`, ")")),
             N = OBS_CT,
             `Effect sizec (standard error)` = paste0(str_pad(round(BETA, n_round), 
                                                              n_round + 1 + str_count(BETA, "[-.]"), 
                                                              "right", pad = "0"), " (",
                                                      str_pad(round(SE, n_round), 
                                                              n_round + 1 + str_count(SE, "[-.]"), 
                                                              "right", pad = "0"), ")"),
             `P-value` = format(P, scientific = TRUE)
      ) %>%
      select(Outcome:`P-value`)
  } else {
    file <- tibble()
  }
  if (i == 1) all_files <- file else all_files <- rbind(all_files, file)
}
write_csv(all_files, "Table_S8.csv")