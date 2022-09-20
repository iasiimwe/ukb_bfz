# Below are project-specific functions used during analysis.


# 1. Change to long format function
#----------------------------------
# Aim: uses the pivot_longer() function {tidyr} to "lengthen" UKBiobank datasets.
# Input(s): a UKBiobank dataset (bd_enc.csv) converted (using ukbconv) - one participant per row,
#           the same variable can have multiple columns.  
#         : a dictionary to identify the Field IDs (Data_Dictionary_Showcase.csv).
# Output: a UKBiobank dataset with more rows (one participant has multiple rows that correspond to 
#           the different instances/assessments) and fewer columns (each variable is present in its
#           own column).
ukb_reshape_long <- function(data, dict){
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
  for(k in seq_along(fields)){
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
    if(k == 1) 
      data_long <- dataQ 
    else 
      data_long <- left_join(data_long, dataQ)   #Join by "f.eid" and "I"
  }
  return(data_long)
}

# 2. Minimum sample size function
#--------------------------------
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
  ukb_power_diff <- function(sample_size){
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
#------------------------
# Aim: from highly related individuals (Kinship > 0.0884, greater than 3rd-degree relatedness), 
#      to exclude those with lower call rates/higher missingness levels.
# Input(s): relatedness dataset,containing the headers ID1 and ID2 (HetHet, IBS0, and Kinship removed
#           after the related individuals are identified).           
#         : missingness dataset.
#         : a dataset containing individuals with the desired phenotype e.g. bd_BFZ.df contains 
#           bendroflumethiazide-treated participants.
# Output: a dataset containing IDs to be included (of the related individuals, those with lower 
#         missingness levels).
ukb_related <- function(related.df, missingness.df, bd_BFZ.df){
  ukb_missingness <- function(id, missingness.df){
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
    for(i in seq_along(ids)) output[[i]] <- ukb_missingness(ids[[i]], missingness.df)
    return(as.character(ids[[which.min(output)]]))
  } else {
    return("not_in_bd_BFZ")
  }
}

# 4. Recode "Prefer not to answer" function
#------------------------------------------
# Aim: change "Prefer not to answer" to NA.      
# Input(s): factor with or without the level "Prefer not to answer".
# Output: factor with "Prefer not to answer" entries changed to NA.
ukb_recode_factor_na <- function(x){
  for (i in seq_along(x)) {
    x[[i]] <- x[[i]] %>% recode_factor("Prefer not to answer" = NA_character_)
  }
  return(x)
}

# 5. Statistics summary function
#-------------------------------
# Aim: to obtain summary statistics.      
# Input(s): numeric vectors (for mean and median) or factors (for percentages).
#         : functions, namely "mean", "median" and "percentage".
#         : number of required decimal places (default is 1).
# Output: formatted means (with sds), medians (with IQRs and ranges), and percentages. If missing 
#         data exists, the percentage of missing data is also outputted. 
ukb_statistics <- function(dataQ, f, round_digits = 1) {
  if(f == "mean") {
    statistic <- paste(format(round(mean(dataQ, na.rm = TRUE), round_digits), nsmall = round_digits), 
                       " (", 
                       format(round(sd(dataQ, na.rm = TRUE), round_digits), nsmall = round_digits), 
                       ")", 
                       sep = ""
    )
    return(statistic)
  } else if(f == "median") {
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
    if(TRUE %in% is.na(dataQ)) {
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
  } else if(f == "percentage") {
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
    if(TRUE %in% is.na(dataQ)) {
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
#--------------------------------
# Aim: Produce summary statistics for a dataset.      
# Input(s): dataset with numerical vectors and/or factors.
#         : number of requires decimal places (default 1).
# Output: dataset with three columns (variable, statistic, and value).
ukb_summary <- function(data, round_digits = 1){
  names_list <- colnames(data)
  output <- vector("list", length(names_list))
  for(i in seq_along(names_list)) {
    name_list <- names_list[[i]]
    dataQ <- data[[name_list]] # or: data %>% select(name_list)
    if(is.numeric(dataQ) == TRUE){
      values1 <- ukb_statistics(dataQ, f = "mean", round_digits)
      values2 <- ukb_statistics(dataQ, f = "median", round_digits)
      name_levels <- rep(name_list, 2)
      statistic <- c("mean (sd)", "median (IQR; range)")
      if(TRUE %in% is.na(dataQ)) {
        name_levels <- c(name_levels, paste(name_list, "missing", sep = ", "))
        statistic <- c(statistic, "n (%)")
      }
      variable_summary <- tibble(variable = name_levels,
                                 statistic = statistic,
                                 value = c(values1, values2)
      )
    } else if(is.factor(dataQ) == TRUE){
      name_levels <- levels(dataQ)
      if(TRUE %in% is.na(dataQ)) name_levels <- c(name_levels, "missing")
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
  
  for(i in seq_along(output)){
    if(i == 1) {
      out <- output[[i]]
    } else {
      out <- rbind(out, output[[i]])
    }
  }
  return(out)
}

# 7. P-value function
#--------------------
# Aim: To derive p-values.    
# Input(s): an 'all_data' summary that combines outputs of the descriptive summary function (has 
#           four columns: variable, statistic, included [summary], excluded [summary]).
#         : 'included' dataset that contains data for participants included in a specific analysis.
#         : 'excluded' dataset that contains data for participants excluded from a specific analysis.
#         : number of requires decimal places (default is 4).
# Output: T-test for means; Mann Whitney for medians; and, Fisher's Exact Test (N < 1000) and 
#         Chi-square (N ≥ 1000) for categorical variables.
ukb_p_value <- function(all_data, included, excluded, round_digits = 4){
  out <- vector("list", length(included))
  for (x in seq_along(included)){
    variable <- colnames(included)[[x]]
    index <- sapply(str_split(all_data$variable, ", "), `[`, 1) == variable
    variable_data <- all_data[index, ]
    if(is.numeric(included[[variable]]) == TRUE) {
      t_test_mean <- round(t.test(included[[variable]],
                                  excluded[[variable]])$p.value, 
                           round_digits)
      mann_whitney_test_median <- round(wilcox.test(included[[variable]],
                                                    excluded[[variable]])$p.value, 
                                        round_digits)
      out[[x]] <- c(t_test_mean, mann_whitney_test_median)
    } else if(is.factor(included[[variable]]) == TRUE) {
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
#-----------------------------------
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
#----------------------
# Aim: Functions to make project-specific plotting easier and consistent.      
# Input(s): a dataset containing individuals with the desired phenotype e.g. bd_BFZ.df contains 
#           bendroflumethiazide-treated participants.
#         : the variable to be plotted (continuous for box/qq plots; categorical for bar plots).
# Output: variable-specific box, bar or quantile-quantile (qq) plots.
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

ukb_qq_plots <- function(data, race = NULL){
  for (i in seq_along(colnames(data))){
    data_covariate <- data[,i]
    data_covariate <- data_covariate[!is.na(data_covariate)]
    qqnorm(data_covariate, pch = 1, frame = FALSE,
           main = paste(race, ": ", tolower(ukb_title_names(colnames(data)[i])), sep = "")
    )
    qqline(data_covariate, col = "steelblue", lwd = 2)
  }
}

# 10. Recode to SNPTest format
#-----------------------------
# Aim: code categorical variables to either "1" or "2", a format that SNPTest requires.    
# Input(s): factors with two levels.
# Output: vector with 1s and 2s.
ukb_recode_snptest <- function(x){
  for (y in names(x)) {
    x[y] <- ifelse(x[[y]] == 0, 1, 2)
  }
  return(x)
}

# 11. Comma function
#-------------------
# Aim: display large numbers separated with commas.
# Input(s): numeric vector.
# Output: character vector.
ukb_comma <- function(x) format(x, big.mark = ",")

# 12. Colour format function
#---------------------------
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
#-----------------------
# Aim: to produce project-specific huxtables for inclusion in rmarkdown output documents.
# Input(s): data frame having summary statistics for the three main race categories (columns are 
#           "Independent variables", "statistic", "White", "Black" and "Asian").
#         : list of covariates e.g. "age" to include in the table (covariate_list).
#         : labels e.g. "Age (years)" that correspond to the above covariate list (covariate_list_2)
#         : the numbers of analyzed participants (analyzed_white, analyzed_black and analyzed_asians)
#           are assumed to have been computed earlier. If not, include them in the environment.
# Output: ukb_hux outputs a huxtable for participant characteristics stratified by the three races 
#         while ukb_hux_2 shows those included and excluded from a specific analysis (race-stratified).
ukb_hux <- function(race_summary_all, covariate_list, covariate_list_2){
  for (x in seq_along(covariate_list)){
    variable <- covariate_list[[x]]
    index <- sapply(str_split(race_summary_all[[1]], ", "), `[`, 1) == variable
    variable_data <- race_summary_all[index, ]
    if(nrow(variable_data) > 0){
      if(TRUE %in% str_detect(variable_data$statistic, "mean")) {
        if (nrow(variable_data) == 2) {
          hux_var <- as_hux(variable_data) %>%
            set_contents(1, 1, "Variables") %>%
            set_contents(1, 3, paste0(.[1, 3], " (<i>N</i> = ", ukb_comma(analyzed_white), ")")) %>%
            set_contents(1, 4, paste0(.[1, 4], " (<i>N</i> = ", ukb_comma(analyzed_black), ")")) %>%
            set_contents(1, 5, paste0(.[1, 5], " (<i>N</i> = ", ukb_comma(analyzed_asian), ")")) %>%
            merge_cells(1, 1:2) %>%
            set_header_rows(1, TRUE) %>%
            style_headers(bold = TRUE, text_color = "black", 
                          background_color = "light blue") %>%
            set_all_borders(value = 0.8) %>%
            merge_cells(2:3, 1) %>%
            merge_cells(2:3, 2) %>%
            set_contents(2, 1, covariate_list_2[[x]]) %>%
            set_contents(2, 2, "Mean (SD)\nMedian (IQR; range)") %>%
            set_contents(2, 3, paste0(.[2, 3], "\n", .[3, 3])) %>%
            merge_cells(2:3, 3) %>%
            set_contents(2, 4, paste0(.[2, 4], "\n", .[3, 4])) %>%
            merge_cells(2:3, 4) %>%
            set_contents(2, 5, paste0(.[2, 5], "\n", .[3, 5])) %>%
            merge_cells(2:3, 5)
        } else {
          hux_var <- as_hux(variable_data) %>%
            set_contents(1, 1, "Variables") %>%
            set_contents(1, 3, paste0(.[1, 3], " (N = ", ukb_comma(analyzed_white), ")")) %>%
            set_contents(1, 4, paste0(.[1, 4], " (N = ", ukb_comma(analyzed_black), ")")) %>%
            set_contents(1, 5, paste0(.[1, 5], " (N = ", ukb_comma(analyzed_asian), ")")) %>%
            merge_cells(1, 1:2) %>%
            set_header_rows(1, TRUE) %>%
            style_headers(bold = TRUE, text_color = "black", 
                          background_color = "light blue") %>%
            set_all_borders(value = 0.8) %>%
            merge_cells(2:4, 1) %>%
            merge_cells(2:4, 2) %>%
            set_contents(2, 1, covariate_list_2[[x]]) %>%
            set_contents(2, 2, 
                         "Mean (SD)\nMedian (IQR; range)\nMissing, n (%)") %>%
            set_contents(2, 3, paste0(.[2, 3], "\n", .[3, 3], "\n", .[4, 3])) %>%
            merge_cells(2:4, 3) %>%
            set_contents(2, 4, paste0(.[2, 4], "\n", .[3, 4], "\n", .[4, 4])) %>%
            merge_cells(2:4, 4) %>%
            set_contents(2, 5, paste0(.[2, 5], "\n", .[3, 5], "\n", .[4, 5])) %>%
            merge_cells(2:4, 5)
        }
        if(x == 1) hux_all <- hux_var else hux_all <- add_rows(hux_all, hux_var[-1,])
      } else {
        factor_levels <- str_to_title(sapply(str_split(variable_data[[1]], ", "), `[`, 2))
        if (nrow(variable_data) == 2) {
          hux_var <- as_hux(variable_data) %>%
            set_contents(1, 1, "Variables") %>%
            set_contents(1, 3, paste0(.[1, 3], " (N = ", ukb_comma(analyzed_white), ")")) %>%
            set_contents(1, 4, paste0(.[1, 4], " (N = ", ukb_comma(analyzed_black), ")")) %>%
            set_contents(1, 5, paste0(.[1, 5], " (N = ", ukb_comma(analyzed_asian), ")")) %>%
            merge_cells(1, 1:2) %>%
            set_header_rows(1, TRUE) %>%
            style_headers(bold = TRUE, text_color = "black", 
                          background_color = "light blue") %>%
            set_all_borders(value = 0.8) %>%
            merge_cells(2:3, 1) %>%
            merge_cells(2:3, 2) %>%
            set_contents(2, 1, covariate_list_2[[x]]) %>%
            set_contents(2, 2, paste(factor_levels, collapse = "\n")) %>%
            set_contents(2, 3, paste0(.[2, 3], "\n", .[3, 3])) %>%
            merge_cells(2:3, 3) %>%
            set_contents(2, 4, paste0(.[2, 4], "\n", .[3, 4])) %>%
            merge_cells(2:3, 4) %>%
            set_contents(2, 5, paste0(.[2, 5], "\n", .[3, 5])) %>%
            merge_cells(2:3, 5)
        } else {
          hux_var <- as_hux(variable_data) %>%
            set_contents(1, 1, "Variables") %>%
            set_contents(1, 3, paste0(.[1, 3], " (N = ", ukb_comma(analyzed_white), ")")) %>%
            set_contents(1, 4, paste0(.[1, 4], " (N = ", ukb_comma(analyzed_black), ")")) %>%
            set_contents(1, 5, paste0(.[1, 5], " (N = ", ukb_comma(analyzed_asian), ")")) %>%
            merge_cells(1, 1:2) %>%
            set_header_rows(1, TRUE) %>%
            style_headers(bold = TRUE, text_color = "black", 
                          background_color = "light blue") %>%
            set_all_borders(value = 0.8) %>%
            merge_cells(2:5, 1) %>%
            merge_cells(2:5, 2) %>%
            set_contents(2, 1, covariate_list_2[[x]]) %>%
            set_contents(2, 2, paste(factor_levels, collapse = "\n")) %>%
            set_contents(2, 3, paste0(.[2, 3], "\n", .[3, 3], "\n", .[4, 3], "\n", .[5, 3])) %>%
            merge_cells(2:5, 3) %>%
            set_contents(2, 4, paste0(.[2, 4], "\n", .[3, 4], "\n", .[4, 4], "\n", .[5, 4])) %>%
            merge_cells(2:5, 4) %>%
            set_contents(2, 5, paste0(.[2, 5], "\n", .[3, 5], "\n", .[4, 5], "\n", .[5, 5])) %>%
            merge_cells(2:5, 5)
        }
        if(x == 1) hux_all <- hux_var else hux_all <- add_rows(hux_all, hux_var[-1,])
      }
    }
  }
  return(hux_all)
}

ukb_hux_2 <- function(data, covariate_list, covariate_list_2, outcome){
  data <- data[[outcome]]
  for (x in seq_along(covariate_list)){
    variable <- covariate_list[[x]]
    index <- sapply(str_split(data[[1]], ", "), `[`, 1) == variable
    variable_data <- data[index, ]
    if(nrow(variable_data) > 0){
      if(TRUE %in% str_detect(variable_data$statistic, "mean")) {
        hux_var <- as_hux(variable_data) %>%
          insert_row(after = 0, colspan = length(.), fill = "") %>%
          set_contents(1, 1, "Variables") %>%
          set_contents(1, 3, paste0("White (<i>N</i> = ", ukb_comma(analyzed_white), ")")) %>%
          set_contents(1, 6, paste0("Black (<i>N</i> = ", ukb_comma(analyzed_black), ")")) %>%
          set_contents(1, 9, paste0("Asian (<i>N</i> = ", ukb_comma(analyzed_asian), ")")) %>%
          merge_cells(1:2, 1:2) %>%
          merge_cells(1, 3:5) %>%
          merge_cells(1, 6:8) %>%
          merge_cells(1, 9:11) %>%
          set_contents(2, 3, paste0("Included (<i>N</i> = ", 
                                    ukb_comma(nrow_included_outcomes[[outcome]][1]), ")")) %>%
          set_contents(2, 4, paste0("Excluded (<i>N</i> = ", 
                                    ukb_comma(nrow_excluded_outcomes[[outcome]][1]), ")")) %>%
          set_contents(2, 5, paste0("<i>P</i> value<sup>a</sup>")) %>%
          set_contents(2, 6, paste0("Included (<i>N</i> = ", 
                                    ukb_comma(nrow_included_outcomes[[outcome]][2]), ")")) %>%
          set_contents(2, 7, paste0("Excluded (<i>N</i> = ", 
                                    ukb_comma(nrow_excluded_outcomes[[outcome]][2]), ")")) %>%
          set_contents(2, 8, paste0("<i>P</i> value<sup>a</sup>")) %>%
          set_contents(2, 9, paste0("Included (<i>N</i> = ", 
                                    ukb_comma(nrow_included_outcomes[[outcome]][3]), ")")) %>%
          set_contents(2, 10, paste0("Excluded (<i>N</i> = ", 
                                     ukb_comma(nrow_excluded_outcomes[[outcome]][3]), ")")) %>%
          set_contents(2, 11, paste0("<i>P</i> value<sup>a</sup>")) %>%
          set_header_rows(1:2, TRUE) %>%
          style_headers(bold = TRUE, text_color = "black", 
                        background_color = "light blue") %>%
          set_all_borders(value = 0.8) %>%
          merge_cells(3:4, 1) %>%
          merge_cells(3:4, 2) %>%
          set_contents(3, 1, covariate_list_2[[x]]) %>%
          set_contents(3, 2, "Mean (SD)\nMedian (IQR; range)") %>%
          set_contents(3, 3, paste0(.[3, 3], "\n", .[4, 3])) %>%
          merge_cells(3:4, 3) %>%
          set_contents(3, 4, paste0(.[3, 4], "\n", .[4, 4])) %>%
          merge_cells(3:4, 4) %>%
          set_contents(3, 5, paste0(.[3, 5], "\n", .[4, 5])) %>%
          merge_cells(3:4, 5)  %>%
          set_contents(3, 6, paste0(.[3, 6], "\n", .[4, 6])) %>%
          merge_cells(3:4, 6) %>%
          set_contents(3, 7, paste0(.[3, 7], "\n", .[4, 7])) %>%
          merge_cells(3:4, 7) %>%
          set_contents(3, 8, paste0(.[3, 8], "\n", .[4, 8])) %>%
          merge_cells(3:4, 8)  %>%
          set_contents(3, 9, paste0(.[3, 9], "\n", .[4, 9])) %>%
          merge_cells(3:4, 9) %>%
          set_contents(3, 10, paste0(.[3, 10], "\n", .[4, 10])) %>%
          merge_cells(3:4, 10) %>%
          set_contents(3, 11, paste0(.[3, 11], "\n", .[4, 11])) %>%
          merge_cells(3:4, 11)
        if(x == 1) hux_all <- hux_var else hux_all <- add_rows(hux_all, hux_var[-1:-2,])
      } else {
        factor_levels <- str_to_title(sapply(str_split(variable_data[[1]], ", "), `[`, 2))
        if (nrow(variable_data) == 2) {
          hux_var <- as_hux(variable_data) %>%
            insert_row(after = 0, colspan = length(.), fill = "") %>%
            set_contents(1, 1, "Variables") %>%
            set_contents(1, 3, paste0("White (<i>N</i> = ", ukb_comma(analyzed_white), ")")) %>%
            set_contents(1, 6, paste0("Black (<i>N</i> = ", ukb_comma(analyzed_black), ")")) %>%
            set_contents(1, 9, paste0("Asian (<i>N</i> = ", ukb_comma(analyzed_asian), ")")) %>%
            merge_cells(1:2, 1:2) %>%
            merge_cells(1, 3:5) %>%
            merge_cells(1, 6:8) %>%
            merge_cells(1, 9:11) %>%
            set_contents(2, 3, paste0("Included (<i>N</i> = ", 
                                      ukb_comma(nrow_included_outcomes[[outcome]][1]), ")")) %>%
            set_contents(2, 4, paste0("Excluded (<i>N</i> = ", 
                                      ukb_comma(nrow_excluded_outcomes[[outcome]][1]), ")")) %>%
            set_contents(2, 5, paste0("<i>P</i> value<sup>a</sup>")) %>%
            set_contents(2, 6, paste0("Included (<i>N</i> = ", 
                                      ukb_comma(nrow_included_outcomes[[outcome]][2]), ")")) %>%
            set_contents(2, 7, paste0("Excluded (<i>N</i> = ", 
                                      ukb_comma(nrow_excluded_outcomes[[outcome]][2]), ")")) %>%
            set_contents(2, 8, paste0("<i>P</i> value<sup>a</sup>")) %>%
            set_contents(2, 9, paste0("Included (<i>N</i> = ", 
                                      ukb_comma(nrow_included_outcomes[[outcome]][3]), ")")) %>%
            set_contents(2, 10, paste0("Excluded (<i>N</i> = ", 
                                       ukb_comma(nrow_excluded_outcomes[[outcome]][3]), ")")) %>%
            set_contents(2, 11, paste0("<i>P</i> value<sup>a</sup>")) %>%
            set_header_rows(1:2, TRUE) %>%
            style_headers(bold = TRUE, text_color = "black", 
                          background_color = "light blue") %>%
            set_all_borders(value = 0.8) %>%
            merge_cells(3:4, 1) %>%
            merge_cells(3:4, 2) %>%
            set_contents(3, 1, covariate_list_2[[x]]) %>%
            set_contents(3, 2, paste(factor_levels, collapse = "\n")) %>%
            set_contents(3, 3, paste0(.[3, 3], "\n", .[4, 3])) %>%
            merge_cells(3:4, 3) %>%
            set_contents(3, 4, paste0(.[3, 4], "\n", .[4, 4])) %>%
            merge_cells(3:4, 4) %>%
            set_contents(3, 5, paste0(.[3, 5])) %>%
            merge_cells(3:4, 5)  %>%
            set_contents(3, 6, paste0(.[3, 6], "\n", .[4, 6])) %>%
            merge_cells(3:4, 6) %>%
            set_contents(3, 7, paste0(.[3, 7], "\n", .[4, 7])) %>%
            merge_cells(3:4, 7) %>%
            set_contents(3, 8, paste0(.[3, 8])) %>%
            merge_cells(3:4, 8)  %>%
            set_contents(3, 9, paste0(.[3, 9], "\n", .[4, 9])) %>%
            merge_cells(3:4, 9) %>%
            set_contents(3, 10, paste0(.[3, 10], "\n", .[4, 10])) %>%
            merge_cells(3:4, 10) %>%
            set_contents(3, 11, paste0(.[3, 11])) %>%
            merge_cells(3:4, 11)
        } else {
          hux_var <- as_hux(variable_data) %>%
            insert_row(after = 0, colspan = length(.), fill = "") %>%
            set_contents(1, 1, "Variables") %>%
            set_contents(1, 3, paste0("White (<i>N</i> = ", ukb_comma(analyzed_white), ")")) %>%
            set_contents(1, 6, paste0("Black (<i>N</i> = ", ukb_comma(analyzed_black), ")")) %>%
            set_contents(1, 9, paste0("Asian (<i>N</i> = ", ukb_comma(analyzed_asian), ")")) %>%
            merge_cells(1:2, 1:2) %>%
            merge_cells(1, 3:5) %>%
            merge_cells(1, 6:8) %>%
            merge_cells(1, 9:11) %>%
            set_contents(2, 3, paste0("Included (<i>N</i> = ", 
                                      ukb_comma(nrow_included_outcomes[[outcome]][1]), ")")) %>%
            set_contents(2, 4, paste0("Excluded (<i>N</i> = ", 
                                      ukb_comma(nrow_excluded_outcomes[[outcome]][1]), ")")) %>%
            set_contents(2, 5, paste0("<i>P</i> value<sup>a</sup>")) %>%
            set_contents(2, 6, paste0("Included (<i>N</i> = ", 
                                      ukb_comma(nrow_included_outcomes[[outcome]][2]), ")")) %>%
            set_contents(2, 7, paste0("Excluded (<i>N</i> = ", 
                                      ukb_comma(nrow_excluded_outcomes[[outcome]][2]), ")")) %>%
            set_contents(2, 8, paste0("<i>P</i> value<sup>a</sup>")) %>%
            set_contents(2, 9, paste0("Included (<i>N</i> = ", 
                                      ukb_comma(nrow_included_outcomes[[outcome]][3]), ")")) %>%
            set_contents(2, 10, paste0("Excluded (<i>N</i> = ", 
                                       ukb_comma(nrow_excluded_outcomes[[outcome]][3]), ")")) %>%
            set_contents(2, 11, paste0("<i>P</i> value<sup>a</sup>")) %>%
            set_header_rows(1:2, TRUE) %>%
            style_headers(bold = TRUE, text_color = "black", 
                          background_color = "light blue") %>%
            set_all_borders(value = 0.8) %>%
            merge_cells(3:5, 1) %>%
            merge_cells(3:5, 2) %>%
            set_contents(3, 1, covariate_list_2[[x]]) %>%
            set_contents(3, 2, paste(factor_levels, collapse = "\n")) %>%
            set_contents(3, 3, paste0(.[3, 3], "\n", .[4, 3], "\n", .[5, 3])) %>%
            merge_cells(3:5, 3) %>%
            set_contents(3, 4, paste0(.[3, 4], "\n", .[4, 4], "\n", .[5, 4])) %>%
            merge_cells(3:5, 4) %>%
            set_contents(3, 5, paste0(.[3, 5])) %>%
            merge_cells(3:5, 5)  %>%
            set_contents(3, 6, paste0(.[3, 6], "\n", .[4, 6], "\n", .[5, 6])) %>%
            merge_cells(3:5, 6) %>%
            set_contents(3, 7, paste0(.[3, 7], "\n", .[4, 7], "\n", .[5, 7])) %>%
            merge_cells(3:5, 7) %>%
            set_contents(3, 8, paste0(.[3, 8])) %>%
            merge_cells(3:5, 8)  %>%
            set_contents(3, 9, paste0(.[3, 9], "\n", .[4, 9], "\n", .[5, 9])) %>%
            merge_cells(3:5, 9) %>%
            set_contents(3, 10, paste0(.[3, 10], "\n", .[4, 10], "\n", .[5, 10])) %>%
            merge_cells(3:5, 10) %>%
            set_contents(3, 11, paste0(.[3, 11])) %>%
            merge_cells(3:5, 11)
        }
        if(x == 1) hux_all <- hux_var else hux_all <- add_rows(hux_all, hux_var[-1:-2,])
      }
    }
  }
  return(hux_all)
}

# 14. Subsetting function
#------------------------
# Aim: to subset dataframe based on race or covariates.     
# Input(s): races ("White". "Black", "Asian") or covariates (e.g. "sbp", "glucose").
# Output: Subsetted dataframe.
ukb_subset <- function(bd_BFZ_df, specific_race, covariate = NULL) {
  if (is.null(covariate) == TRUE) return(bd_BFZ_df %>% filter(race == specific_race) %>% nrow()) 
  else return(bd_BFZ_df %>% filter(race == specific_race) %>% select(all_of(covariate)) %>% na.omit() %>% nrow())
}