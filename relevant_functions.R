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
      select(contains(c("f.eid", str_c("f.", fields[k], sep="")))) %>%
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

# 4. Recode "Prefer not to answer" function
# -----------------------------------------
# Aim: change "Prefer not to answer" to NA.      
# Input(s): factor with or without the level "Prefer not to answer".
# Output: factor with "Prefer not to answer" entries changed to NA.
ukb_recode_factor_na <- function(x) {
  for (i in seq_along(x)) {
    x[[i]] <- x[[i]] %>% recode_factor("Prefer not to answer" = NA_character_)
  }
  return(x)
}

# 5. Statistics summary function
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

# 6. Descriptive summary function
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

# 7. P-value function
# -------------------
# Aim: To derive p-values.    
# Input(s): an 'all_data' summary that combines outputs of the descriptive summary function (has 
#           four columns: variable, statistic, included [summary], excluded [summary]).
#         : 'included' dataset that contains data for participants included in a specific analysis.
#         : 'excluded' dataset that contains data for participants excluded from a specific analysis.
#         : number of requires decimal places (default is 4).
# Output: T-test for means; Mann Whitney for medians; and, Fisher's Exact Test (N < 1000) and 
#         Chi-square (N ≥ 1000) for categorical variables.
ukb_p_value <- function(all_data, included, excluded, round_digits = 4) {
  out <- vector("list", length(included))
  for (x in seq_along(included)) {
    variable <- colnames(included)[[x]]
    index <- sapply(str_split(all_data$variable, ", "), `[`, 1) == variable
    variable_data <- all_data[index, ]
    if (is.numeric(included[[variable]]) == TRUE) {
      t_test_mean <- round(t.test(included[[variable]],
                                  excluded[[variable]])$p.value, 
                           round_digits)
      mann_whitney_test_median <- round(wilcox.test(included[[variable]],
                                                    excluded[[variable]])$p.value, 
                                        round_digits)
      out[[x]] <- c(t_test_mean, mann_whitney_test_median)
    } else if (is.factor(included[[variable]]) == TRUE) {
      analysis_matrix <- matrix(c(as.numeric(
        sapply(str_split(variable_data$included, " "), `[`, 1)),
        as.numeric(sapply(str_split(variable_data$excluded, " "), `[`, 1))), 
        ncol = 2)
      if (sum(analysis_matrix) < 1000 ) {
        fisher_chi_test <- round(fisher.test(analysis_matrix)$p.value, 
                                 round_digits)
      } else {
        fisher_chi_test <- round(chisq.test(analysis_matrix)$p.value, 
                                 round_digits)
      }
      out[[x]] <- c(fisher_chi_test, rep(NA, times = nrow(variable_data) - 1))
      # Fisher's Exact Test (equivalent to Chi-square with large sample sizes, more accurate with 
      # small sample sizes, PMID: 28503482). http://www.biostathandbook.com/small.html recommends
      # always using an exact test if the total sample size is less than 1000.
    } else {
      out[[x]] <- "na (neither numeric nor factor)"
    }
  }
  return(unlist(out))
}

# 8. Label and title names functions
# ----------------------------------
# Aim: Provide labels and titles for plotting graphs, limited to the variables used during analysis, 
#      but can be accordingly expanded.
# Input(s): the variable being plotted.
# Output: label to use in a plot.
ukb_label_names <- function(covariate) {
  switch(covariate,
         age = "Age (years)",
         sex = "Sex",
         bmi = expression(Kg/m^2),
         smoking_status = "Smoking status",
         alcohol_status = "Alcohol status",
         townsend_index = "Townsend index score",
         pa = "Types of physical activity in last 4 weeks",
         glucose = "mmol/L",
         urate = expression(paste(mu, "mol/L", sep = "")),
         sbp = "mmHg",
         dbp = "mmHg",
         stop("Unknown covariate!"))
}

ukb_label_names_2 <- function(covariate) {
  switch(covariate,
         age = "Age (years)",
         sex = "Sex",
         bmi = expression(Kg/m^2),
         smoking_status = "Smoking status",
         alcohol_status = "Alcohol status",
         townsend_index = "Townsend index score",
         pa = "Types of physical activity in last 4 weeks",
         glucose = "mmol/L",
         urate = "umol/L",
         sbp = "mmHg",
         dbp = "mmHg",
         stop("Unknown covariate!"))
}

ukb_title_names <- function(covariate) {
  switch(covariate,
         age = "Age",
         sex = "Sex",
         bmi = "BMI",
         smoking_status = "Smoking status",
         alcohol_status = "Alcohol status",
         townsend_index = "Townsend index",
         pa = "Physical activity",
         glucose = "Blood glucose",
         urate = "Blood urate",
         sbp = "Systolic blood pressure",
         dbp = "Diastolic blood pressure",
         stop("Unknown covariate!"))
}

ukb_reference_levels <- function(covariate) {
  switch(covariate,
         age = NULL,
         sex = NULL,
         bmi = "Ideal range: 18.5 to 24.9",
         smoking_status = NULL,
         alcohol_status = NULL,
         townsend_index = NULL,
         pa = NULL,
         glucose = "Normal fasting range: 3.9 to 5.6",
         urate = "Reference range: Males, 200 to 420;\nFemales, 140 to 360",
         sbp = "Ideal range: 90 to 120", # from https://www.nhs.uk/common-health-questions/lifestyle/what-is-blood-pressure/
         dbp = "Ideal range: 60 to 80",
         stop("Unknown covariate!"))
}

# 9. Plotting functions
# ---------------------
# Aim: Functions to make project-specific plotting easier and consistent.      
# Input(s): a dataset containing individuals with the desired phenotype e.g. bd_BFZ.df contains 
#           bendroflumethiazide-treated participants.
#         : the variable to be plotted (continuous for box/qq plots; categorical for bar plots).
# Output: variable-specific box, bar, manhattan or quantile-quantile (qq) plots.
ukb_box_plots <- function(bd_BFZ, covariate) {
  ggplot(data = bd_BFZ, 
         mapping = aes(x = race, y = .data[[covariate]], fill = race)
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

ukb_bar_plots <- function(bd_BFZ, covariate) {
  ggplot(data = bd_BFZ, mapping = aes(x = race, 
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

ukb_qq_plots <- function(data, race = NULL) {
  for (i in seq_along(colnames(data))) {
    data_covariate <- data[,i]
    data_covariate <- data_covariate[!is.na(data_covariate)]
    qqnorm(data_covariate, pch = 1, frame = FALSE,
           main = paste(race, ": ", tolower(ukb_title_names(colnames(data)[i])), sep = "")
    )
    qqline(data_covariate, col = "steelblue", lwd = 2)
  }
}

ukb_manhattan_plots <- function(data, race, outcome, y_max = NULL) {
  y_limit <- max(ceiling(-log10(5e-8)), ceiling(-log10(min(data$P))))
  if (!is.null(y_max)) y_limit <- y_max
  png(paste("UKBB_", race, "_", outcome, "_manhattan.png", sep = ""), 
      width = 1500, height = 800, res = 120)
  qqman::manhattan(data,
                   genomewideline = -log10(5e-8),
                   col = c("blue4", "orange3"),
                   ylim = c(0, y_limit)
                   )
  dev.off()
} 

ukb_qq_gwas_plots <- function(data, race, outcome, y_max = NULL) {
  y_limit <- max(ceiling(-log10(5e-8)), ceiling(-log10(min(data$P))))
  if (!is.null(y_max)) y_limit <- y_max
  png(paste("UKBB_", race, "_", outcome, "_qqplot.png", sep = ""), 
      width = 1500, height = 1500, res = 120)
  qqman::qq(data$P, 
            col = "blue4",
            ylim = c(0, y_limit)
            )
  dev.off()
}

# 10. Recode to SNPTest format
# ----------------------------
# Aim: code categorical variables to either "1" or "2", a format that SNPTest requires.    
# Input(s): factors with two levels.
# Output: vector with 1s and 2s.
ukb_recode_snptest <- function(x) {
  for (y in names(x)) {
    x[y] <- ifelse (x[[y]] == 0, 1, 2)
  }
  return(x)
}

# 11. Comma function
# ------------------
# Aim: display large numbers separated with commas.
# Input(s): numeric vector.
# Output: character vector.
ukb_comma <- function(x) format(x, big.mark = ",")

# 12. Colour format function
# --------------------------
# Aim: change text colour in rmarkdown using raw HTML or LaTeX code.
#      This combines HTML and LaTeX syntax:
#        HTML: <span style="color: red;">text</span>.
#        PDF: \textcolor{}{}.      
# Input(s): rmarkdown text.
# Output: colour-formatted rmarkdown text.
ukb_colour_format <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color,
            x)
  } else x
}

# 13. Huxtable functions
# ----------------------
# Aim: to produce project-specific huxtables for inclusion in rmarkdown output documents.
# Input(s): data frame having summary statistics for the three main race categories (columns are 
#           "Independent variables", "statistic", "White", "Black" and "Asian").
#         : list of covariates e.g. "age" to include in the table (covariate_list).
#         : labels e.g. "Age (years)" that correspond to the above covariate list (covariate_list_2)
#         : the numbers of analyzed participants (analyzed_white, analyzed_black and analyzed_asians)
#           are assumed to have been computed earlier. If not, include them in the environment.
# Output: ukb_hux outputs a huxtable for participant characteristics stratified by the three races 
#         while ukb_hux_2 shows those included and excluded from a specific analysis (race-stratified).
ukb_set_contents <- function(ht, ht_row, ht_cols, value) {
  i <- 1
  while (i <= length(ht_cols)) {
    ht <- ht %>% set_contents(ht_row, ht_cols[[i]], value[[i]])
    i <- i + 1
  }
  return(ht)
}

ukb_set_merge_contents <- function(ht, ht_rows, ht_cols) {
  for (i in seq_along(ht_cols)) {
    ht_col <- ht_cols[[i]]
    nrows <- max(ht_rows) - min(ht_rows) + 1
    ht <- ht %>% 
      set_contents(ht_rows, ht_col,
                   gsub("NA", "", paste0(.[ht_rows, ht_col][[1]],
                                         c(rep("\n", times = (nrows - 1)),
                                           ""), collapse = "")
                        )
                   ) %>% 
      merge_cells(ht_rows, ht_col)
    if (i == 1) ht_all <- ht[, 1:min(ht_col)] else ht_all <- add_columns(ht_all, ht[, ht_col])
  }
  return(ht_all)
}

ukb_hux <- function(race_summary_all, covariate_list, covariate_list_2) {
  for (x in seq_along(covariate_list)) {
    variable <- covariate_list[[x]]
    index <- sapply(str_split(race_summary_all[[1]], ", "), `[`, 1) == variable
    variable_data <- race_summary_all[index, ]
    if (nrow(variable_data) > 0) {
      hux_values <- c("Variables", paste0(c("White", "Black", "Asian"),
                                          " (<i>N</i> = ",
                                          c(ukb_comma(analyzed_white),
                                            ukb_comma(analyzed_black),
                                            ukb_comma(analyzed_asian)
                                            ),
                                          ")"))
      hux_header <- as_hux(variable_data) %>%
        ukb_set_contents(1, c(1, 3:5), hux_values) %>%
        merge_cells(1, 1:2) %>%
        set_header_rows(1, TRUE) %>%
        style_headers(bold = TRUE, text_color = "black", 
                      background_color = "light blue") %>%
        set_all_borders(value = 0.8)
      if (nrow(variable_data) == 2) {
        hux_var_2 <- hux_header %>%
          merge_cells(2:3, 1) %>%
          merge_cells(2:3, 2) %>%
          set_contents(2, 1, covariate_list_2[[x]]) %>%
          ukb_set_merge_contents(2:3, 3:5)
        if (TRUE %in% str_detect(variable_data$statistic, "mean")) {
          hux_var <- hux_var_2 %>%
            set_contents(2, 2, "Mean (SD)\nMedian (IQR; range)")
          } else {
          factor_levels <- str_to_title(sapply(str_split(variable_data[[1]], ", "), `[`, 2))
          hux_var <- hux_var_2 %>%
            set_contents(2, 2, paste(factor_levels, collapse = "\n"))
          }
      } else{
        if (TRUE %in% str_detect(variable_data$statistic, "mean")) {
          hux_var <- hux_header %>%
            merge_cells(2:4, 1) %>%
            merge_cells(2:4, 2) %>%
            set_contents(2, 1, covariate_list_2[[x]]) %>%
            set_contents(2, 2, "Mean (SD)\nMedian (IQR; range)\nMissing, n (%)") %>%
            ukb_set_merge_contents(2:4, 3:5)
        } else {
          factor_levels <- str_to_title(sapply(str_split(variable_data[[1]], ", "), `[`, 2))
          hux_var <- hux_header %>%
            merge_cells(2:5, 1) %>%
            merge_cells(2:5, 2) %>%
            set_contents(2, 1, covariate_list_2[[x]]) %>%
            set_contents(2, 2, paste(factor_levels, collapse = "\n")) %>%
            ukb_set_merge_contents(2:5, 3:5)
        }
      }
      if (x == 1) hux_all <- hux_var else hux_all <- add_rows(hux_all, hux_var[-1,])
    }
    }
  return(hux_all)
}

ukb_hux_2 <- function(data, covariate_list, covariate_list_2, outcome) {
  data <- data[[outcome]]
  for (x in seq_along(covariate_list)) {
    variable <- covariate_list[[x]]
    index <- sapply(str_split(data[[1]], ", "), `[`, 1) == variable
    variable_data <- data[index, ]
    if (nrow(variable_data) > 0) {
      hux_values <- c("Variables", paste0(c("White", "Black", "Asian"),
                                          " (<i>N</i> = ",
                                          c(ukb_comma(analyzed_white),
                                            ukb_comma(analyzed_black),
                                            ukb_comma(analyzed_asian)),
                                          ")"))
      hux_values_2 <- c(paste0(c("Included (<i>N</i> = ",
                                 "Excluded (<i>N</i> = ",
                                 "<i>P</i> value<sup>a</sup>"),
                               c(ukb_comma(nrow_included_outcomes[[outcome]][1]),
                                ukb_comma(nrow_excluded_outcomes[[outcome]][1]), "",
                                ukb_comma(nrow_included_outcomes[[outcome]][2]), 
                                ukb_comma(nrow_excluded_outcomes[[outcome]][2]), "",
                                ukb_comma(nrow_included_outcomes[[outcome]][3]), 
                                ukb_comma(nrow_excluded_outcomes[[outcome]][3]), ""),
                               c(")", ")", "")))
      hux_header <- as_hux(variable_data) %>%
        insert_row(after = 0, colspan = length(.), fill = "") %>%
        ukb_set_contents(1, c(1, 3, 6, 9), hux_values) %>%
        merge_cells(1:2, 1:2) %>%
        ukb_set_contents(2, c(3:11), hux_values_2) %>%
        set_header_rows(1:2, TRUE) %>%
        style_headers(bold = TRUE, text_color = "black", 
                      background_color = "light blue") %>%
        set_all_borders(value = 0.8)
      if (nrow(variable_data) == 2) {
        hux_var_2 <- hux_header %>%
          merge_cells(3:4, 1) %>%
          merge_cells(3:4, 2) %>%
          set_contents(3, 1, covariate_list_2[[x]]) %>%
          ukb_set_merge_contents(3:4, 3:11)
        if (TRUE %in% str_detect(variable_data$statistic, "mean")) {
          hux_var <- hux_var_2 %>%
            set_contents(3, 2, "Mean (SD)\nMedian (IQR; range)")
        } else {
          factor_levels <- str_to_title(sapply(str_split(variable_data[[1]], ", "), `[`, 2))
          hux_var <- hux_var_2 %>%
            set_contents(3, 2, paste(factor_levels, collapse = "\n")) 
        }
      } else {
        factor_levels <- str_to_title(sapply(str_split(variable_data[[1]], ", "), `[`, 2))
        hux_var <- hux_header %>%
          merge_cells(3:5, 1) %>%
          merge_cells(3:5, 2) %>%
          set_contents(3, 1, covariate_list_2[[x]]) %>%
          set_contents(3, 2, paste(factor_levels, collapse = "\n")) %>%
          ukb_set_merge_contents(3:5, 3:11) 
      }
      if (x == 1) hux_all <- hux_var else hux_all <- add_rows(hux_all, hux_var[-1:-2,])
    }
  }
  hux_all <- hux_all %>%
    merge_cells(1, 3:5) %>%
    merge_cells(1, 6:8) %>%
    merge_cells(1, 9:11)
  return(hux_all)
}

# 14. Splitting text and tables
# -----------------------------
# Aim: to enable strings fit in huxtables, and tables on pages     
# Input(s): a string or long table.
# Output: string with "\n" inserted or split tables.
ukb_split_text <- function(text, n_digits = 5) {
  split_text <- strsplit(text, "")[[1]]
  parts <- ceiling(length(split_text)/n_digits)
  for (x in 1:parts){
    y <- n_digits
    start <- (y * (x - 1) + 1)
    end <- min((y * (x - 1) + y), length(split_text))
    if (x == parts) {
      comb_text <- paste0(split_text[start:end], collapse = "")
    } else {
      comb_text <- paste0(paste0(split_text[start:end], collapse = ""), "\n", collapse = "")
    }
    if (x == 1) all_text <- comb_text else all_text <- paste0(all_text, comb_text, collapse = "")
  }
  return(all_text)
}

ukb_split_across <- function(ht, n_rows = 50, print_hux = TRUE) {
  parts <- ceiling((nrow(ht) - 1) / n_rows)
  output <- vector("list", parts)
  for (x in 1:parts){
    y <- n_rows
    start <- (y * (x - 1) + 1) 
    end <- min((y * (x - 1) + y + 1), nrow(ht))
    if (x == 1) {
      ht_2 <- split_across(ht, after = end)[[1]]
    } else if (x == parts) {
      ht_2 <- split_across(ht, after = start)[[2]]
    } else {
      ht_2 <- split_across(ht, after = c(start, end))[[2]]
    }
    output[[x]] <- ht_2
    if (print_hux == TRUE) {
      quick_html(ht_2, file = paste0("hux_table_4", "_", x, "_", race, "_", outcome, ".html"), open = FALSE)
    }
  }
  return(output)
}

# 15. Subsetting function
# -----------------------
# Aim: to subset dataframe based on race or covariates and count number of participants.     
# Input(s): races ("White". "Black", "Asian") or covariates (e.g. "sbp", "glucose").
# Output: a count of analyzed participants.
ukb_subset_count <- function(bd_BFZ_df, specific_race, covariate = NULL) {
  if (is.null(covariate) == TRUE) return(bd_BFZ_df %>% filter(race == specific_race) %>% nrow()) 
  else return(bd_BFZ_df %>% filter(race == specific_race) %>% select(all_of(covariate)) %>% na.omit() %>% nrow())
}

# 16. lambda function
# -------------------
# Aim: compute the genomic inflation factor (lambda).     
# Input(s): gwas results
# Output: lambda.
ukb_lambda <- function(data) {
  chisq = qchisq(data$P, 1, lower.tail = FALSE)
  return(median(chisq) / qchisq(0.5, 1))
}

# 17. Genes/functional significance function
# ------------------------------------------
# Aim: Obtain the genes associated with selected SNPs from the dbSNP database. 
#      (https://www.ncbi.nlm.nih.gov/snp/)
# Input(s): race, outcome and n_digits (assumes a .csv file containing these snps)
#           already exists.
# Output: data frame (.csv and .rds files) that contains associated genes as well 
#         as SNP locations (e.g. introns) 
ukb_dbsnp <- function(race, outcome, n_digits = 3) {
  file <- fread(paste0("UKBB_", race, "_", outcome, ".csv")) 
  N <- unique(file$all_total)
  file <- file %>%
    mutate(Gene = NA,
           `Functional consequence` = NA,
           all_maf = str_pad(round(all_maf, n_digits), 
                             n_digits + 1 + str_count(all_maf, "[-.]"), 
                             "right", 
                             pad = "0"
                             ),
           `Beta (SE)` = paste0(str_pad(round(frequentist_add_beta_1, n_digits), 
                                        n_digits + 1 + str_count(frequentist_add_beta_1, "[-.]"), 
                                        "right", 
                                        pad = "0"),
                                " (",
                                str_pad(round(frequentist_add_se_1, n_digits), 
                                        n_digits + 1 + str_count(frequentist_add_se_1, "[-.]"), 
                                        "right", 
                                        pad = "0"),
                                ")"
                                ),
           info = round(info, n_digits),
           `#` = seq.int(nrow(.))
    ) %>%
    rename(Chr = CHR,
           Position = BP,
           `P value` = P,
           `Ref allele` = alleleA,
           `Alt allele` = alleleB,
           `Info score` = info
    ) %>%
    select(`#`, SNP:Position, Gene, `Functional consequence`, `Ref allele`:`Info score`,
           all_maf, `Beta (SE)`, `P value`)
  colnames(file)[which(colnames(file) == "all_maf")] <- paste0("MAF (N = ", N, ")")
  for (k in 1:nrow(file)) {
    snp <- file$SNP[k]
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
    message(paste(round(k / nrow(file) * 100, 3), "% complete", sep = ""))
    if (k %% 100 == 0) {Sys.sleep(10)} # Pause for 10 seconds after every 100 requests
  }
  write.csv(file,  file = paste0("UKBB_", race, "_", outcome, "_final.csv"), row.names = FALSE)
  write_rds(file,  file = paste0("UKBB_", race, "_", outcome, "_final.rds"))
}

# 18. Getting gene summary
# ------------------------
# Aim: to obtain gene summary from https://www.ncbi.nlm.nih.gov/gene.       
# Input(s): vectors of races and outcomes. 
#           P-value threshold (default 5e-8) for SNPs whose associated genes are 
#           required.
#           Assumes at least one race- and outcome-specific file (.csv output of 
#           ukb_dbsnp()) is available in the environment.
# Output: summary e.g. gene of the genes associated with the top SNPs.
ukb_genes_summary <- function(races, outcomes, p_value = 5e-8) {
  for(i in seq_along(races)) {
    race <- races[[i]]
    for(j in seq_along(outcomes)){
      outcome <- outcomes[[j]]
      genes_race_outcome <- as_tibble(fread(paste0("UKBB_", race, "_", outcome, "_final.csv"))) %>%
        filter (`P value` < p_value) %>%
        select(Gene) %>%
        na.omit() %>%
        unique() 
      genes_race_outcome <- unique(unlist(str_split(genes_race_outcome$Gene, ", ")))
      if(j == 1) genes_race <- genes_race_outcome 
      else genes_race <- unique(bind_rows(genes_race, genes_race_outcome))
    }
    if(i == 1) genes <- genes_race 
    else genes <- unique(bind_rows(genes, genes_race))
  }
  gene_summary <- tibble(Gene = rep("NA", length(genes)),
                         ID = rep("NA", length(genes)),
                         Symbol = rep("NA", length(genes)), 
                         `Gene name` = rep("NA", length(genes)),
                         `Also known as` = rep("NA", length(genes)),
                         Type = rep("NA", length(genes)),
                         Summary = rep("NA", length(genes)),
                         Expression = rep("NA", length(genes)),
                         `Annotation info` = rep("NA", length(genes))
  )
  for(k in seq_along(genes)) {
    gene <- genes[[k]]
    gene_summary$Gene[[k]] <- gene
    
    # obtain gene id
    gene_id <- paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", gene) %>%
      read_html() %>%
      html_nodes(xpath = '//*[@id="padded_content"]') %>%
      html_text() %>%
      gsub("[\n()]|  ", "", .) %>% 
      str_extract("Gene ID: .*PubMed") %>%
      str_extract("Gene ID: .*RefSeq transcripts") %>%
      gsub("Gene ID: |RefSeq transcripts", "", .)
    if (is.na(gene_id) == TRUE && str_detect(gene, "LOC") == TRUE) 
      gene_id <- gsub("LOC", "", gene)
    gene_summary$ID[[k]] <- gene_id
    
    # use gene id to obtain info (preferred to name as names could have changed)
    gene_info <- paste0("https://www.ncbi.nlm.nih.gov/gene/", gene_id) %>%
      read_html() 
    
    # obtain official symbol
    gene_summary$Symbol[[k]] <- gene_info %>%
      html_elements(xpath = '//*[@id="summaryDl"]/dd[1]/text()') %>%
      html_text() 
    
    # obtain gene_name
    gene_name <- gene_info %>%
      html_elements(xpath = '//*[@id="summaryDl"]/dd[2]/text()') %>%
      html_text() 
    gene_summary$`Gene name`[[k]] <- gene_name
    
    gene_info_2 <- gene_info %>%
      html_elements(xpath = '//*[@id="summaryDl"]') %>%
      html_text() %>%
      gsub("[\n]|  ", "", .)
    
    # obtain gene alias
    gene_summary$`Also known as`[[k]] <- gene_info_2 %>% 
      str_extract("Also known as.*Summary") %>%
      gsub("Also known as|Summary", "", .)
    
    # obtain gene type
    gene_summary$Type[[k]] <- gene_info_2 %>% 
      str_extract("Gene type.*RefSeq ") %>%
      gsub("Gene type|RefSeq ", "", .)
    
    # obtain gene summary
    gene_summary$Summary[[k]] <- gene_info_2 %>% 
      str_extract("Summary.*\\[provided by") %>%
      gsub("Summary|\\[provided by", "", .)
    
    # obtain gene expression
    gene_expression <- gene_info_2 %>% 
      str_extract("Expression.*Orthologs") %>%
      gsub("Expression|Orthologs", "", .)
    if (TRUE %in% str_detect(gene_expression, "See more")) {
      gene_expression <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=",
                                gene_id,
                                "&retmode=xml") %>%
        read_html() %>%
        html_nodes(xpath = '/html/body') %>%
        html_text() %>%
        gsub("[\n]|  ", "", .) %>%
        str_extract("Expression[0-9]*Text Summary.*[0-9]*Category") %>%
        gsub("Expression[0-9]*Text Summary|[0-9]*Category", "", .) %>% 
        gsub("254", ". ", .) %>%
        gsub("Tissue List", "Tissue List: ", .)
    }
    gene_summary$Expression[[k]] <- gene_expression
    
    # obtain gene annotation information
    gene_summary$`Annotation info`[[k]] <- gene_info_2 %>% 
      str_extract("Annotation information.*Expression") %>%
      gsub("Annotation information|Expression", "", .)
    
    message(paste(round(k / length(genes) * 100, 3), "% complete", sep = ""))
    if(k %% 100 == 0) {Sys.sleep(10)} # Pause for 10 seconds after every 100 requests
  }
  return(gene_summary)
}

# 19. Image annotation
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

# 20. Pathway analysis 
# --------------------
# Aim: facilitate pathway analysis using Reactome (version 82, https://reactome.org/).
#      Requires the "rbioapi" package.
# Input(s): list of genes associated with the top SNPs, number of required pathways
#           and whether to download a full pdf report or not.
# Output: gene summaries based on the Reactome database (ukb_gene_reactome()).
#         gene-related pathways (ukb_gene_pathways()).
ukb_gene_reactome <- function(genes) {
  gene_reactome <- tibble(query_name = rep("NA", length(genes)),
                          db_id = rep("NA", length(genes)),
                          gene_name = rep("NA", length(genes)),
                          gene_symbol = rep("NA", length(genes)), 
                          gene_description = rep("NA", length(genes)),
                          gene_summary = rep("NA", length(genes))
  )
  for(k in seq_along(genes)) {
    gene <- genes[[k]]
    gene_reactome$query_name[k] <- gene
    
    not_in_database <- tryCatch(
      xref_gene <- rba_reactome_xref(gene),
      error = function(e) e
    )
    
    if(inherits(not_in_database, "error")) next
    
    gene_reactome$db_id[k] <- xref_gene$dbId
    gene_reactome$gene_name[k] <- xref_gene$name[[1]]
    gene_reactome$gene_symbol[k] <- paste0(xref_gene$geneName[[1]], collapse = ", ")
    gene_reactome$gene_description[k] <- gsub("recommendedName: ", "", xref_gene$description[[1]])
    gene_reactome$gene_summary[k] <- gsub("FUNCTION ", "", xref_gene$comment[[1]])
    
    message(paste(round(k / length(genes) * 100, 3), "% complete", sep = ""))
    if(k %% 100 == 0) {Sys.sleep(10)} # Pause for 10 seconds after every 100 requests
  }
  return(gene_reactome)
}

ukb_gene_pathways <- function(genes, top_pathways = 25, download_report = TRUE){
  # Analyze all genes
  analyzed_all_genes <- rba_reactome_analysis(input = genes, species = "Homo sapiens")
  
  # Download a full pdf report
  if (download_report == TRUE) {
    file_name <- "rba_reactome"
    rba_reactome_analysis_pdf(token = analyzed_all_genes$summary$token,
                              species = 9606, 
                              number = top_pathways,
                              save_to = paste0(file_name, ".pdf")
    )
  }
  
  # Obtain a summary of the top pathways
  top_summary <- as_tibble(analyzed_all_genes$pathways[1:top_pathways, ]) %>%
    mutate(`Pathway name` = name,
           reactions.found = paste0(reactions.found, "/", reactions.total),
           entities.found = paste0(entities.found, "/", entities.total)
           ) %>%
    select(stId, `Pathway name`, entities.found, entities.ratio, entities.pValue, 
           entities.fdr, reactions.found, reactions.ratio)
  
  # Check if a gene is related to a pathway
  ukb_gene_paths <- function(genes) {
    gene_pathways <- tibble(query_name = rep("NA", length(genes)),
                            pathways = rep("NA", length(genes))
                            )
    for(k in seq_along(genes)) {
      gene <- genes[[k]]
      gene_pathways$query_name[k] <- gene
      gene_pathways$pathways[k] <- paste0(rba_reactome_analysis(input = gene, 
                                                                species = "Homo sapiens"
      )$pathways[[1]],
      collapse = ", "
      )
      message(paste(round(k / length(genes) * 100, 3), "% complete", sep = ""))
      if(k %% 100 == 0) {Sys.sleep(10)} # Pause for 10 seconds after every 100 requests
    }
    return(gene_pathways)
  }
  gene_pathways <- ukb_gene_paths(genes)
  
  # Check if a pathway is related to a gene
  ukb_check_pathway <- function(pathway, gene_pathways){
    return(gene_pathways$query_name[str_detect(gene_pathways$pathways, pathway)])
  }
  
  # Check genes related to included pathways
  included_genes <- vector("character", length(top_summary$stId))
  for(i in seq_along(top_summary$stId)){
    included_genes[i] <- paste0(ukb_check_pathway(top_summary$stId[[i]], 
                                                  gene_pathways), collapse = ", ")
  }
  
  # Add the related genes to the summary
  top_summary <- top_summary %>%
    mutate(included.genes = included_genes)
  
  return(top_summary)
}