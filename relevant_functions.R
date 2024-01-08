# Below are project-specific functions used during analysis. 


# 1. Change to long format function
# ---------------------------------
# Aim: uses the pivot_longer() function {tidyr} to "lengthen" UKBiobank datasets.
# Input(s): a UKBiobank dataset (bd_enc.csv) converted (using ukbconv) - one participant per row,
#           the same variable can have multiple columns.  
#         : a dictionary to identify the Field IDs (Data_Dictionary_Showcase.csv).
# Output: a UKBiobank dataset with more rows (one participant has multiple rows that correspond to 
#           the different instances/assessments) and fewer columns (each variable is present in its
#           own column).
ukb_reshape_long <- function(data, dict) {
  fields <- data %>% 
    colnames() %>% 
    str_extract_all(., "f.[0-9]*") %>% 
    unlist() %>% 
    str_replace_all(., "f.", "") %>%
    unique() %>%
    parse_integer() %>%
    na.omit()
  fields <- append(fields[!fields %in% 53], 53, after = 0) # remove date field, so that it can be added  
          # to the start of the vector - this enables all assessments to be added using left_join().
  for (k in seq_along(fields)) {
    columns <- unlist(str_extract_all(colnames(data), 
                                      str_c("f.", fields[k], ".[0-9]*.[0-9]*", 
                                            sep=""))) # columns to pivot into longer format
    values_name <- dict %>% filter(FieldID == fields[k]) %>% .$Field
    dataQ <- data %>% 
      select(contains(c("f.eid", str_c("f.", fields[k], ".", sep="")))) %>%
      pivot_longer(cols = all_of(columns), 
                   names_to = "I", 
                   values_to = values_name, 
                   values_drop_na = TRUE
                   ) %>%
      arrange(f.eid) %>%
      mutate(I = str_replace(.$I, str_c("f.", fields[k], sep=""), "")) %>% 
      mutate(I = unlist(str_extract(.$I, ".[0-9].")))
    if (k == 1) 
      data_long <- dataQ 
    else 
      data_long <- left_join(data_long, dataQ)   #Join by "f.eid" and "I"
  }
  return(data_long)
}


# 2. Relatedness function
# -----------------------
# Aim: from highly related individuals (Kinship > 0.0884, greater than 3rd-degree relatedness), 
#      to exclude those with lower call rates/higher missingness levels.
# Input(s): relatedness dataset,containing the headers ID1 and ID2 (HetHet, IBS0, and Kinship removed
#           after the related individuals are identified).           
#         : missingness dataset.
#         : a dataset containing individuals with the desired phenotype e.g. bd_BFZ.df contains 
#           bendroflumethiazide-treated participants.
# Output: a dataset containing IDs to be included (of the related individuals, those with lower 
#         missingness levels).
ukb_related <- function(related.df, missingness.df, bd_BFZ.df) {
  ukb_missingness <- function(id, missingness.df) {
    colnames(missingness.df)[1:2] <- c("ID1", "missingness")
    id_missing <- missingness.df %>% 
      filter(ID1 == id) %>% 
      select(missingness)
    return(id_missing$missingness)
  }
  colnames(missingness.df)[1:2] <- c("ID1", "missingness")
  missingness.df$missingness[is.na(missingness.df$missingness) == TRUE] <- 1
  ids <- related.df %>% 
    as.matrix() %>% 
    as.integer() %>% 
    unique() 
  ids <- ids[ids %in% bd_BFZ.df$f.eid]
  if (length(ids) > 0) {
    output <- vector("double", length(ids))
    for (i in seq_along(ids)) output[[i]] <- ukb_missingness(ids[[i]], missingness.df)
    return(as.character(ids[[which.min(output)]]))
  } else {
    return("not_in_bd_BFZ")
  }
}


# 3. Label and title names functions
# ----------------------------------
# Aim: Provide labels and titles for plotting graphs, limited to the variables used during analysis, 
#      but can be accordingly expanded.
# Input(s): the variable being plotted.
# Output: label to use in a plot.

ukb_label_names <- function(covariate) {
  switch(covariate,
         age = "Age (years)",
         sex = "Sex",
         array = "Genotyping array",
         bmi = expression(Kg/m^2),
         smoking = "Smoking status",
         alcohol = "Alcohol status",
         townsend = "Townsend index score",
         pa = "Types of physical activity in last 4 weeks",
         glucose = "mmol/L",
         urate = expression(paste(mu, "mol/L", sep = "")),
         potassium = "mmol/L",
         sodium = "mmol/L",
         thiazide = "Thiazides",
         losartan = "Losartan",
         raas = "Non-losartan RAASi",
         bb = "Beta-blockers",
         ccb = "Calcium channel blockers",
         allo = "Allopurinol",
         stop("Unknown covariate!"))
}

ukb_label_names_2 <- function(covariate) {
  switch(covariate,
         age = "Age (years)",
         sex = "Sex",
         array = "Genotyping array",
         bmi = expression(Kg/m^2),
         smoking = "Smoking status",
         alcohol = "Alcohol status",
         townsend = "Townsend index score",
         pa = "Types of physical activity in last 4 weeks",
         glucose = "mmol/L",
         urate = "umol/L",
         potassium = "mmol/L",
         sodium = "mmol/L",
         thiazide = "Thiazides",
         losartan = "Losartan",
         raas = "Non-losartan RAASi",
         bb = "Beta-blockers",
         ccb = "Calcium channel blockers",
         allo = "Allopurinol",
         stop("Unknown covariate!"))
}

ukb_title_names <- function(covariate) {
  switch(covariate,
         age = "Age",
         sex = "Sex",
         array = "Genotyping array",
         bmi = "BMI",
         smoking = "Smoking status",
         alcohol = "Alcohol status",
         townsend = "Townsend index",
         pa = "Physical activity",
         glucose = "Blood glucose",
         urate = "Blood urate",
         potassium = "Urine potassium",
         sodium = "Urine sodium",
         thiazide = "Thiazides",
         losartan = "Losartan",
         raas = "Non-losartan RAASi",
         bb = "Beta-blockers",
         ccb = "Calcium channel blockers",
         allo = "Allopurinol",
         stop("Unknown covariate!"))
}


# 4. Plotting functions
# ---------------------
# Aim: Functions to make project-specific plotting easier and consistent.      
# Input(s): a dataset containing individuals with the desired phenotype e.g. bd_BFZ.df contains 
#           bendroflumethiazide-treated participants.
#         : the variable to be plotted (continuous for box/qq plots; categorical for bar plots).
# Output: variable-specific box, bar, manhattan or quantile-quantile (qq) plots.
ukb_box_plots <- function(bd_BFZ, covariate, grouping = "thiazide") {
  ggplot(data = bd_BFZ, 
         mapping = aes(x = .data[[grouping]], y = .data[[covariate]], fill = .data[[grouping]])
  ) +
    geom_boxplot(colour = "grey20", alpha = 0.7) +
    labs(title = ukb_title_names(covariate),
         x = NULL,
         y = ukb_label_names(covariate)
    ) +
    theme(axis.text.x = element_text(colour = "black"), 
          axis.text.y = element_text(colour = "black"), 
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold.italic"),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.position = "none"  # remove key/legend
    ) +
    scale_fill_brewer(palette = "Set1")
}

ukb_bar_plots <- function(bd_BFZ, covariate, grouping = "thiazide") {
  ggplot(data = bd_BFZ, mapping = aes(x = .data[[grouping]], 
                                      fill = fct_rev(.data[[covariate]]))) +
    geom_bar(position = "fill", colour = "grey20", alpha = 0.7) +
    labs(title = ukb_title_names(covariate),
         x = NULL,
         y = "Percentage",
         fill = ukb_title_names(covariate)
    ) +
    theme(axis.text.x = element_text(colour = "black"), 
          axis.text.y = element_text(colour = "black"), 
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold.italic"),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.position = "right",  # remove key/legend
          legend.title = element_text(face = "italic")
    ) +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00),
                       labels = c("0%", "25%", "50%", "75%", "100%")
    )
}

ukb_manhattan_plots <- function(data, cohort, outcome, type, y_max = NULL,
                                genomewideP = 5e-8, suggestive = -log10(1e-05), 
                                snps_to_highlight = NUL) {
  y_limit <- max(ceiling(-log10(5e-8)), ceiling(-log10(min(data$P))))
  if (!is.null(y_max)) y_limit <- y_max
  png(paste0("UKBB_", cohort, "_", outcome, "_manhattan_", type, ".png"), 
      width = 1500, height = 800, res = 120)
  qqman::manhattan(data,
                   genomewideline = -log10(genomewideP),
                   col = c("blue", "deepskyblue3"),
                   ylim = c(0, y_limit),
                   suggestiveline = suggestive, 
                   highlight = snps_to_highlight
                   )
  dev.off()
} 

ukb_qq_gwas_plots <- function(data, cohort, outcome, type, y_max = NULL) {
  y_limit <- max(ceiling(-log10(5e-8)), ceiling(-log10(min(data$P))))
  if (!is.null(y_max)) y_limit <- y_max
  png(paste0("UKBB_", cohort, "_", outcome, "_qqplot_", type, ".png"), 
      width = 1500, height = 1500, res = 120)
  qqman::qq(data$P, 
            col = "blue",
            ylim = c(0, y_limit)
            )
  dev.off()
}


# 5. lambda function
# -------------------
# Aim: compute the genomic inflation factor (lambda).     
# Input(s): gwas results
# Output: lambda.
ukb_lambda <- function(data) {
  chisq = qchisq(data$P, 1, lower.tail = FALSE)
  return(median(chisq) / qchisq(0.5, 1))
}


# 6. Pre-adjusting phenotype
# --------------------------
# Aim: pre-adjust the phenotype based on age, sex and genetic principal components
# Input(s): dataset, outcome, sd (for removing outliers)
# Output: pre-adjusted phenotype (residuals)

# Primary analysis
ukb_model_outcome <- function(data, outcome, sd_number = 5) {
  dataQ <- data %>%
    select(-glucose, -urate, -potassium, -sodium)
  dataQ$y <- data[[outcome]]
  model_res <- paste("C", 1:10, sep = "", collapse = " + ") %>%
    paste0("y ~ age +", .) %>%
    as.formula %>%
    lm(., data = dataQ) %>%
    residuals()
  max_res <- mean(model_res) + 5 * sd(model_res)
  min_res <- mean(model_res) - 5 * sd(model_res)
  
  dataQ <- left_join(dataQ[,c("f.eid", "sex")], data.frame(f.eid = dataQ$f.eid[complete.cases(dataQ)], residuals = model_res))
  # filter out values less than 5 times or more than 5 times the sd
  dataQ_sd <- dataQ %>%
    filter(residuals >= min_res & residuals <= max_res)
  # standardized to z scores with  mean 0 and variance 1 in each gender group
  dataQ_sd$std_res <- NA
  dataQ_sd$std_res[dataQ_sd$sex == "Male"] <- scale(dataQ_sd$residuals[dataQ_sd$sex == "Male"])
  dataQ_sd$std_res[dataQ_sd$sex == "Female"] <- scale(dataQ_sd$residuals[dataQ_sd$sex == "Female"])
  
  dataQ <- dataQ %>%
    select(f.eid) %>%
    left_join(., dataQ_sd)
  
  return(dataQ$std_res)
}

# Sebsitivity analysis
ukb_model_outcome2 <- function(data, outcome, sd_number = 5) {
  dataQ <- data %>%
    select(-glucose, -urate, -potassium, -sodium)
  dataQ$y <- data[[outcome]]
  model_res <- paste("C", 1:40, sep = "", collapse = " + ") %>%
    paste0("y ~ age + bmi + smoking + alcohol + pa + townsend + array + losartan + raas + bb + ccb + allo +", .) %>%
    as.formula %>%
    lm(., data = dataQ) %>%
    residuals()
  max_res <- mean(model_res) + 5 * sd(model_res)
  min_res <- mean(model_res) - 5 * sd(model_res)
  
  dataQ <- left_join(dataQ[,c("f.eid", "sex")], data.frame(f.eid = dataQ$f.eid[complete.cases(dataQ)], residuals = model_res))
  # filter out values less than 5 times or more than 5 times the sd
  dataQ_sd <- dataQ %>%
    filter(residuals >= min_res & residuals <= max_res)
  # standardized to z scores with  mean 0 and variance 1 in each gender group
  dataQ_sd$std_res <- NA
  dataQ_sd$std_res[dataQ_sd$sex == "Male"] <- scale(dataQ_sd$residuals[dataQ_sd$sex == "Male"])
  dataQ_sd$std_res[dataQ_sd$sex == "Female"] <- scale(dataQ_sd$residuals[dataQ_sd$sex == "Female"])
  
  dataQ <- dataQ %>%
    select(f.eid) %>%
    left_join(., dataQ_sd)
  
  return(dataQ$std_res)
}


# 7. Recode "Prefer not to answer" function
# ------------------------------------------
# Aim: change "Prefer not to answer" to NA.      
# Input(s): factor with or without the level "Prefer not to answer".
# Output: factor with "Prefer not to answer" entries changed to NA.
ukb_recode_factor_na <- function(x) {
  for (i in seq_along(x)) {
    x[[i]] <- x[[i]] %>% recode_factor("Prefer not to answer" = NA_character_)
  }
  return(x)
}


# 8. Image annotation
# --------------------
# Aim: annotate images, specifically manhattan plots.     
# Input(s): .png image.
# Output: annotated image.
ukb_image_annotate <- function(image, text, x = 0, y = 0,  size = 20, color = "black", 
                               boxcolor = "transparent", degrees = 0, font = 'Times',
                               style = "italic", bold = FALSE, decoration = NULL, kerning = 0)
{
  if (bold == FALSE) weight = 400 else weight = 700 
  image_annotate(image = image, text = text, size = size, color = color, boxcolor = boxcolor, degrees = degrees, 
                 location = paste0("+", x, "+", y), font = font, style = style, 
                 weight = weight, decoration = decoration, kerning = kerning)
}


# 9. Genes/functional significance function
# ------------------------------------------
# Aim: Obtain the genes associated with selected SNPs from the dbSNP database. 
#      (https://www.ncbi.nlm.nih.gov/snp/)
# Input(s): cohort, outcome and n_digits (assumes a .csv file containing these snps)
#           already exists.
# Output: data frame (.csv and .rds files) that contains associated genes as well 
#         as SNP locations (e.g. introns) 
ukb_dbsnp <- function(snps_list, n_digits = 3) {
  file <- tibble(SNP = rep(NA, length(snps_list)),
                 Gene = rep(NA, length(snps_list)),
                 `Functional consequence` = rep(NA, length(snps_list)),
                 )
  for (k in seq_along(snps_list)) {
    snp <- snps_list[[k]]
    file$SNP[k] <- snp
    if (str_detect(snp, ":") == FALSE) {
      SNP_info <- paste0("https://www.ncbi.nlm.nih.gov/snp/?term=", snp) %>%
        read_html() %>%
        html_nodes("dl.snpsum_dl_left_align") %>%
        html_text() %>%
        gsub("[\n() ]", "", .)
      file$Gene[k] <- SNP_info[1] %>% 
        str_extract("Gene:.*Varview") %>%
        gsub("Gene:|Varview", "", .) %>%
        gsub(",", ", ", .)
      file$`Functional consequence`[k] <- SNP_info[1] %>% 
        str_extract("FunctionalConsequence:.*Validated") %>%
        gsub("FunctionalConsequence:|Validated", "", .) %>%
        gsub(",", ", ", .) %>%
        gsub("_|:", " ", .) %>%
        gsub("non coding", "non-coding", .) %>%
        gsub("Clinicalsignificance", "; Clinical significance:", .) %>%
        gsub(" transcript| variant", "", .) %>%
        gsub("prime", "-prime", .)
    }
    message(paste(round(k / length(snps_list) * 100, 3), "% complete", sep = ""))
    if (k %% 100 == 0) {Sys.sleep(10)} # Pause for 10 seconds after every 100 requests
  }
  return(file)
}

# 10. Additive coding function
# ----------------------------
# Aim: Code genotypes using the additive mode of inheritance
# Input(s): vector of genotypes.
# Output: Coded genotypes, with wild-type homozygotes = 0, heterozygotes = 1, and mutant-type homozygotes = 2. 
additive_coding_fn <- function(x) {
  x <- gsub(" ", "", x) 
  x[x == "00"] <- NA
  y <- table(x)
  y <- names(sort(-y)) # Sort in descending order, assume mutant type is the least abundant
  if (length(y) == 3) x <- gsub(y[[3]], "2", x) 
  if (length(y) > 1) x <- gsub(y[[2]], "1", x) 
  x <- gsub(y[[1]], "0", x)
  return(as.numeric(x))
}

# 10. Minor allele frequency (MAF)
# -------------------------------
# Aim: Obtain MAFs
# Input(s): vector of genotypes.
# Output: MAF 
maf_fn <- function(x) {
  x <- table(x)
  if (length(x) == 1) maf <- 1 
  if (length(x) == 2) maf <- as.numeric(x["1"] / (2 * ( x["0"] +  x["1"])))
  if (length(x) == 3) maf <- as.numeric((x["1"] + 2 *x["2"])/ (2 * ( x["0"] +  x["1"] +  x["2"])))
  return(maf)
}