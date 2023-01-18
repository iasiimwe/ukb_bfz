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


# 2. Minimum sample size function
# -------------------------------
# Aim: Compute minimum sample-size (continuous outcomes), requires the 'pwr' package.      
# Input(s): required effect size (obtained from literature and/or expert opinion).
#         : known/extrapolated outcome standard deviation for the study population.
#         : statistical significance (p-value) threshold (default for GWAS is 5e-8).
#         : required power (default is 80%).
#         : assumed minor allele frequency (MAF, default is 20%).
#         : n_times, number of times the minimum sample size is assumed to decrease with MAF = 50%.
# Output: colour-formatted rmarkdown text.
ukb_min_sample <- function(effect, std_ev, p = 5e-8, power = 0.8, maf = 0.2, n_times = 2) {
  equal_sizes <- ceiling((pwr.t.test(d = effect / std_ev, 
                                     sig.level = p, 
                                     alternative = "two.sided", 
                                     power = power)$n)
                         * 2) # sample size assuming equally-sized groups (i.e. MAF = 50%)
  unequal_sizes <- equal_sizes * n_times # assuming having unequally-sized group will increase sample size by at most n_times
  ukb_power_diff <- function(sample_size) {
    current_power <- (pwr.t2n.test(n1 = maf * sample_size, 
                                   n2 = (1 - maf) * sample_size, 
                                   d = effect / std_ev, 
                                   sig.level = p, 
                                   alternative = "two.sided"))$power
    abs(current_power - power)
  }
  min_sample_size <- ceiling(optimize(ukb_power_diff, c(equal_sizes:unequal_sizes))$minimum)
  if (min_sample_size < floor(unequal_sizes)) return(min_sample_size) 
  else return("Assumed maximum sample size likely low, increase n_times by 1")
}


# 3. Relatedness function
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


# 4. Statistics summary function
# ------------------------------
# Aim: to obtain summary statistics.      
# Input(s): numeric vectors (for mean and median) or factors (for percentages).
#         : functions, namely "mean", "median" and "percentage".
#         : number of required decimal places (default is 1).
# Output: formatted means (with sds), medians (with IQRs and ranges), and percentages. If missing 
#         data exists, the percentage of missing data is also outputted. 
ukb_statistics <- function(dataQ, f, round_digits = 1) {
  if (f == "mean") {
    statistic <- paste(format(round(mean(dataQ, na.rm = TRUE), round_digits), nsmall = round_digits), 
                       " (", 
                       format(round(sd(dataQ, na.rm = TRUE), round_digits), nsmall = round_digits), 
                       ")", 
                       sep = ""
    )
    return(statistic)
  } else if (f == "median") {
    quantiles <- format(round(quantile(dataQ, na.rm = TRUE), digits = round_digits), 
                        nsmall = round_digits)
    statistic <- paste(quantiles[[3]],
                       " (",
                       quantiles[[2]],
                       "–",
                       quantiles[[4]],
                       "; ",
                       quantiles[[1]],
                       "–",
                       quantiles[[5]],
                       ")",
                       sep = ""
    )
    if (TRUE %in% is.na(dataQ)) {
      percent_na <- paste(sum(is.na(dataQ)),
                          " (",
                          format(round(sum(is.na(dataQ)) * 100 / length(dataQ), round_digits), 
                                 nsmall = round_digits),
                          "%)",
                          sep = ""
      )
      statistic <- c(statistic, percent_na)
    }
    return(statistic)
  } else if (f == "percentage") {
    out <- vector("double", length(levels(dataQ)))
    for (i in seq_along(levels(dataQ))) {
      out[[i]] <- paste(table(dataQ)[i],
                        " (",
                        format(round(table(dataQ)[i] * 100 / length(dataQ), 
                                     round_digits), 
                               nsmall = round_digits),
                        "%)",
                        sep = ""
      )
    }
    if (TRUE %in% is.na(dataQ)) {
      percent_na <- paste(sum(is.na(dataQ)), 
                          " (",
                          format(round(sum(is.na(dataQ)) * 100 / length(dataQ), 
                                       round_digits), 
                                 nsmall = round_digits),
                          "%)", 
                          sep = ""
      )
      out <- c(out, percent_na)
    }
    return(out)
  } else {
    return("function not specified")
  }
}


# 5. Descriptive summary function
# -------------------------------
# Aim: Produce summary statistics for a dataset.      
# Input(s): dataset with numerical vectors and/or factors.
#         : number of requires decimal places (default 1).
# Output: dataset with three columns (variable, statistic, and value).
ukb_summary <- function(data, round_digits = 1) {
  names_list <- colnames(data)
  output <- vector("list", length(names_list))
  for (i in seq_along(names_list)) {
    name_list <- names_list[[i]]
    dataQ <- data[[name_list]] # or: data %>% select(name_list)
    if (is.numeric(dataQ) == TRUE) {
      values1 <- ukb_statistics(dataQ, f = "mean", round_digits)
      values2 <- ukb_statistics(dataQ, f = "median", round_digits)
      name_levels <- rep(name_list, 2)
      statistic <- c("mean (sd)", "median (IQR; range)")
      if (TRUE %in% is.na(dataQ)) {
        name_levels <- c(name_levels, paste(name_list, "missing", sep = ", "))
        statistic <- c(statistic, "n (%)")
      }
      variable_summary <- tibble(variable = name_levels,
                                 statistic = statistic,
                                 value = c(values1, values2)
      )
    } else if (is.factor(dataQ) == TRUE) {
      name_levels <- levels(dataQ)
      if (TRUE %in% is.na(dataQ)) name_levels <- c(name_levels, "missing")
      values <- ukb_statistics(dataQ, f = "percentage", round_digits)
      variable_summary <- tibble(variable = paste(name_list, 
                                                  name_levels, sep = ", "),
                                 statistic = rep("n (%)", length(name_levels)),
                                 value = values
      )
    } else {
      variable_summary <- tibble(variable = name_list,
                                 statistic = "na (neither numeric nor factor)",
                                 value = "na (neither numeric nor factor)"
      )
    }
    output[[i]] <- variable_summary
  }
  
  for (i in seq_along(output)) {
    if (i == 1) {
      out <- output[[i]]
    } else {
      out <- rbind(out, output[[i]])
    }
  }
  return(out)
}


# 6. Label and title names functions
# ----------------------------------
# Aim: Provide labels and titles for plotting graphs, limited to the variables used during analysis, 
#      but can be accordingly expanded.
# Input(s): the variable being plotted.
# Output: label to use in a plot.
ukb_label_names <- function(covariate) {
  switch(covariate,
         age = "Age (years)",
         sex = "Sex",
         glucose = "mmol/L",
         urate = expression(paste(mu, "mol/L", sep = "")),
         potassium = "mmol/L",
         sodium = "mmol/L",
         thiazide = "Thiazide status",
         stop("Unknown covariate!"))
}

ukb_label_names_2 <- function(covariate) {
  switch(covariate,
         age = "Age (years)",
         sex = "Sex",
         glucose = "mmol/L",
         urate = "umol/L",
         potassium = "mmol/L",
         sodium = "mmol/L",
         thiazide = "Thiazide status",
         stop("Unknown covariate!"))
}

ukb_title_names <- function(covariate) {
  switch(covariate,
         age = "Age",
         sex = "Sex",
         glucose = "Blood glucose",
         urate = "Serum urate",
         potassium = "Urine potassium",
         sodium = "Urine sodium",
         thiazide = "Thiazide",
         cohort = "Cohort",
         stop("Unknown covariate!"))
}

ukb_reference_levels <- function(covariate) {
  switch(covariate,
         age = NULL,
         sex = NULL,
         glucose = "Normal fasting range: 3.9 to 5.6",
         urate = "Reference range: Males, 200 to 420;\nFemales, 140 to 360",
         potassium = NULL,
         sodium = NULL,
         thiazide = NULL,
         stop("Unknown covariate!"))
}


# 7. Plotting functions
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
         y = ukb_label_names(covariate),
         caption = ukb_reference_levels(covariate)
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

ukb_qq_plots <- function(data, cohort = NULL) {
  for (i in seq_along(colnames(data))) {
    data_covariate <- data[,i]
    data_covariate <- data_covariate[!is.na(data_covariate)]
    qqnorm(data_covariate, pch = 1, frame = FALSE,
           main = paste(cohort, ": ", tolower(ukb_title_names(colnames(data)[i])), sep = "")
    )
    qqline(data_covariate, col = "steelblue", lwd = 2)
  }
}

ukb_manhattan_plots <- function(data, cohort, outcome, y_max = NULL) {
  y_limit <- max(ceiling(-log10(5e-8)), ceiling(-log10(min(data$P))))
  if (!is.null(y_max)) y_limit <- y_max
  png(paste("UKBB_", cohort, "_", outcome, "_manhattan.png", sep = ""), 
      width = 1500, height = 800, res = 120)
  qqman::manhattan(data,
                   genomewideline = -log10(5e-8),
                   col = c("blue", "deepskyblue3"),
                   ylim = c(0, y_limit)
                   )
  dev.off()
} 

ukb_qq_gwas_plots <- function(data, cohort, outcome, y_max = NULL) {
  y_limit <- max(ceiling(-log10(5e-8)), ceiling(-log10(min(data$P))))
  if (!is.null(y_max)) y_limit <- y_max
  png(paste("UKBB_", cohort, "_", outcome, "_qqplot.png", sep = ""), 
      width = 1500, height = 1500, res = 120)
  qqman::qq(data$P, 
            col = "blue",
            ylim = c(0, y_limit)
            )
  dev.off()
}


# 8. Recode to SNPTest format
# ----------------------------
# Aim: code categorical variables to either "1" or "2"    
# Input(s): factors with two levels.
# Output: vector with 1s and 2s.
ukb_recode_snptest <- function(x) {
  for (y in names(x)) {
    x[y] <- ifelse (x[[y]] == 0, 1, 2)
  }
  return(x)
}


# 9. lambda function
# -------------------
# Aim: compute the genomic inflation factor (lambda).     
# Input(s): gwas results
# Output: lambda.
ukb_lambda <- function(data) {
  chisq = qchisq(data$P, 1, lower.tail = FALSE)
  return(median(chisq) / qchisq(0.5, 1))
}