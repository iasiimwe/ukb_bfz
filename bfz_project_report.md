---
title: "ukb_bfz"
author: "Innocent G Asiimwe"
date: "20 Sep 2022"
output:
  html_document:
    keep_md: true
bibliography: references.bib
csl: frontiers-in-pharmacology.csl
---

<style type="text/css">
  body{
  font-size: 12pt;
}
</style>



### <ins>**Genome-wide association studies of Bendroflumethiazide-related outcomes in the UK Biobank**</ins>

<br>

#### **Introduction**
* Thiazide diuretics, including Bendroflumethiazide (BFZ), are among the most widely-used, most effective and least costly therapies for hypertension. 
* Many thiazide-treated patients do not attain adequate blood pressure (BP) response while many experience electrolyte disorders that can result in adverse outcomes, including mortality. 
* <span style='color: blue;'>*We will therefore conduct genome-wide association studies (GWASs) using BFZ-treated patients in the UK Biobank to identify the genomic factors that influence BP response and the risk of two electrolyte disorders (hyperglycaemia and hyperuricaemia).*</span>

<br>

\newpage
#### **Methods**
$\underline{Study\space design\space and\space setting}$

* UK Biobank, a population-based prospective cohort study [@RN1079].
  + Approximately half a million UK individuals, aged between 40–69 at recruitment that occurred between 2006–2010. 
  + Comprehensive genetic (>93 million markers) and phenotypic (including linked primary health care records) data. 

<br>

$\underline{Study\space population}$

* Prevalent and newly-initiated BFZ users who will be identified using data fields and coding algorithms provided by the UK Biobank.
  + Prevalent users (current analysis) identified using the UK Biobank **Data-Field 20003** (“Treatment/medication code” containing self-reported medications).
  + BFZ given singly or as part of a combination.

<br>

$\underline{Outcomes}$

* Co-primary/efficacy outcomes: 
  + Systolic BP (automated, **Data-Field 4080**; manual, **Data-Field 93**)
  + Diastolic BP (automated, **Data-Field 4079**; manual, **Data-Field 94**) 
* Secondary safety/electrolyte disorder outcomes: 
  + Blood glucose (<b>Data-Field 30740</b>)
  + Serum urate (<b>Data-Field 30880</b>)
  + Hyponatraemia (serum sodium) and hypokalaemia (serum potassium) not measured at baseline. We may later be able to include them using primary care records.
* Current analysis is continuous outcomes, to later consider dichotomous outcomes.

<br>

$\underline{GWAS\space covariates}$

* First ten principal components of genetic ancestry (**Data-Field 22009**)
* Non-genetic covariates:
  + Age at baseline (years, <b>Data-Field 21003</b>)
  + Sex (male vs female, **Data-Field 31**)
  + Body mass index (Kg/m^2^, **Data-Field 21001**)
  + Smoking status (current vs former vs never, **Data-Field 20116**)
  + Alcohol drinker status (current vs former vs never, **Data-Field 20117**)
  + Degree of physical activity (Types of physical activity in last 4 weeks, coded to heavy [a) strenuous sports and b) heavy DIY like weeding, lawn mowing, carpentry, digging], light [a) other exercises like swimming, cycling, keep fit, bowling, b) light DIY like pruning, watering the lawn, and c) walking for pleasure (not as a means of transport) and none, **Data-Field 6164**)
  + Townsend index reflecting socioeconomic status (a continuous score, **Data-Field 189**)

<br>

$\underline{Genotyping\space and\space QC}$

* Two similar (about 95% shared marker content) arrays: 
  + UK BiLEVE Axiom Array by Affymetrix (49 950 participants, 807 411 markers)
  + Applied Biosystems UK Biobank Axiom Array (438 427 participants, 825 927 markers)
  + Phasing (SHAPEIT3) and imputation (IMPUTE4) using UK10K, 1000 Genomes phase 3, and Haplotype Reference Consortium panels increased the number of markers to > 93 million.
* Per-SNP QC: SNPs included if MAF ≥ 1%, missingness ≤ 5%, HWE > 1E-06, and info score > 0.4.
* Per-sample QC: removed participants who:
  + Were outliers based on missing rate (>5%) and heterozygosity adjusted for population structure,
  + Had mismatches between self-reported and marker-inferred sex, and/or,
  + Had greater than 3rd-degree relatedness with another included participant.

<br>

$\underline{Minimum\space sample\space size}$

* Based on the primary outcomes of systolic and diastolic BPs.
* Assuming a MAF of 20%, 80% power and a significance threshold of 5 × 10-8:
  + Systolic BP: <span style='color: blue;'>909</span> participants, standard deviation (SD) of 19 mmHg for 475 540 UK Biobank participants and effect size 10 mmHg (obtained from previous studies)
  + Diastolic BP: <span style='color: blue;'>1213</span> participants, SD of 11 mmHg (UK Biobank) and effect size 5 mmHg (previous studies)
* These sample sizes achievable for Whites. Due to smaller sample sizes, the proposal is to use other races (Blacks and Asians) for replication.

<br>

$\underline{Statistical\space analysis}$

* Outcome transformation: log transformation for continuous outcomes.
* Predictor handling
  + Quantitative predictor variables neither transformed nor categorized. 
  + Nonlinearity between continuous predictors and continuous outcomes assessed using restricted cubic splines in the rms R package.
  + Categorical variables dummy coded in R.
* Missing data
  + Participants with missing outcome data excluded from the corresponding analyses.
  + No participant excluded on the basis of missing covariate data. Missing genotype data imputed using IMPUTE4, missing non-genetic covariates imputed using single imputation (MICE R package).
* Multivariable linear regression
  + Additive mode of inheritance assumption
  + Adjustment for the seven non-genetic covariates and the first 10 principal components of genetic ancestry.
  + Two significance thresholds: a genome-wide statistical significance threshold of p < 5 × 10-8, as well as a nominal significance threshold of p < 1 × 10-5. 
  + No further adjustment for multiplicity testing as this analysis is exploratory. 
  + All analyses will be undertaken using SNPTEST Version 2. 
* Further analysis using the Single Nucleotide Polymorphism (dbSNP) database, the Genotype-Tissue-Expression (GTEx) portal, LocusZoom and Haploview.

<br>

\newpage
#### **Results**













<br>


**Figure 1. Participant inclusion flow chart.** to be inserted here.



<br>



**Figure 2. Participant characteristics (continuous independent covariates plus outcomes).** to be inserted here.

<br>


**Figure 3. Participant characteristics (categorical independent covariates).** to be inserted here.

<br>




**Table 1** to be inserted here.

<br>




**Table 2** to be inserted here.

<br>




**Table 3** to be inserted here.

<br>


**Figure 4. QQ plots for non-transformed outcomes.** to be inserted here.

<br>


**Figure 5. QQ plots for log-transformed outcomes.** to be inserted here.

<br>


**Figure 6. An enhanced scatterplot matrix showing pairwise correlations between the variables.** to be inserted here.

<br>


**Figure 7. Strip plot showing imputed continuous covariates.** to be inserted here.

<br>


**Figure 8. Non-linearity assessmnet using restricted cubic splines.** to be inserted here.

<br>



<br>

\newpage
#### **References**

