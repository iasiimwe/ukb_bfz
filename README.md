# ukb_bfz
A genome-wide association study of Bendroflumethiazide-related outcomes in the UK Biobank. 

# UK Biobank data
Details on accessing UK Biobank data are available at https://biobank.ctsu.ox.ac.uk/~bbdatan/Accessing_UKB_data_v2.3.pdf. Validating, decrypting and converting the downloaded UK Biobank file needs to be done using command line instructions in a command prompt in Windows or a terminal window in Linux. The main steps include:
1. Download the helper programs/file handlers (e.g. ukbmd5, ukbconv, ukbunpack) and the miscellaneous utility 'encoding.ukb' from https://biobank.ndph.ox.ac.uk/ukb/download.cgi and change the file handler permissions (e.g. using **chmod 755**).
2. Download the encrypted dataset (e.g. ukb12345.enc) through the UK Biobank Access Management System.
3. Verify the integrity of the downloaded file using a command similar to: **./ukbmd5 ukb12345.enc**
4. Decrypt the dataset using ukbunpack: **./ukbunpack ukb12345.enc k67890r12345.key** (this step requires a keyfile that is provided to approved users by the UK Biobank in a notification email).
5. Convert the dataset to a desired output (e.g. docs to create a data dictionary for the dataset or csv/txt/r/sas/stata to create files that can be used in the respective programs). 
    - To generate an R file for one field, run a command similar to: **./ukbconv ukb12345.enc_ukb r -s22009** (Field ID 22009 corresponds to principal components of genetic ancestry). 
    - For multiple fields, use a command like: **./ukbconv ukb12345.enc_ukb r -ifields.txt**, where fields.txt contains the required fields.
    - Two files are outputted and one (ukb12345.tab) contains the data in a tab-separated format. The second (ukb23456.R) can be opened and executed in an R environment to decode all categorical variables in ukb12345.tab. It can be edited to output a comma separated variable (bd_enc.csv) version of the data that is then used for downstream analysis.
    - Note that a .csv file could have been produced using the command: **./ukbconv ukb12345.enc_ukb csv -ifields.txt**, however, the data codings are retained, instead of being replaced by their meanings. 
    - For current analysis, the following fields were used:
       * Co-primary/efficacy outcomes: Systolic BP (automated, **Data-Field 4080**; manual, **Data-Field 93**), Diastolic BP (automated, **Data-Field 4079**; manual, **Data-Field 94**).
       * Secondary safety/electrolyte disorder outcomes: Blood glucose (<b>Data-Field 30740</b>), Serum urate (<b>Data-Field 30880</b>).
       * Genetic covariates: the first ten principal components of genetic ancestry (**Data-Field 22009**).
       * Non-genetic covariates: Age at baseline (years, <b>Data-Field 21003</b>), Sex (male vs female, **Data-Field 31**), Body mass index (Kg/m<sup>2</sup>, **Data-Field 21001**), Smoking status (current vs former vs never, **Data-Field 20116**), Alcohol drinker status (current vs former vs never, **Data-Field 20117**), Degree of physical activity (Types of physical activity in last 4 weeks, **Data-Field 6164**), and Townsend index reflecting socioeconomic status (a continuous score, **Data-Field 189**).
6. Details of downloading genetics data are included in section 4 of https://biobank.ctsu.ox.ac.uk/~bbdatan/Accessing_UKB_data_v2.3.pdf. In addition to the bulk genotyped and/or imputed genetic data, the following should be downloaded.
    - .fam files, obtained using **./gfetch 22418 -c1 -m -ak67890r12345.key** for genotype data (**Data-Field 22418**) or **./gfetch 22828 -c1 -m -k67890r12345.key** for imputed data (**Data-Field 22828**). 
    - Relatedness dataset, obtained using **./gfetch rel -ak67890r12345.key**.
    - Missingness dataset (**Data-Field 22005**).
    - Outliers for heterozygosity or missing rate (**Data-Field 22027**).
    - A dictionary (Data_Dictionary_Showcase.csv) to identify the Field IDs (downloaded from https://biobank.ctsu.ox.ac.uk/~bbdatan/Data_Dictionary_Showcase.csv).
    - A data codings file (available at: https://biobank.ctsu.ox.ac.uk/~bbdatan/Codings.csv).

# Downstream analysis
Before analysis, ensure to install the following R files using commands similar to:
   ```{r install required packages, include = FALSE}
   required <- c("data.table", "tidyverse", "lubridate", "mice", "qqman", "gridExtra", "grid", "psych", "rms", "pwr", "huxtable")
   install.packages(required)
   ```

For downstream analysis, two files are provided in this repository:
* A relevant_functions.R file that contains project-specific functions used during analysis.
* An Rmarkdown file (bfz_project_report.Rmd) that includes the project report and R code used during analysis. Before this work is published, the draft report will include only the code used during analysis, without the corresponding results.
