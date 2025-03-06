#  Assessing Spatial Disparities: A Bayesian Linear Regression Approach

This repository contains code implementation for detecting spatial disparities in a BYM2 model using the methodology discussed in the reference below.  

Reference:
> Wu, K. L., and Banerjee, S., "Assessing Spatial Disparities: A Bayesian Linear Regression Approach." Manuscript in preparation.

A real data analysis example is included to showcase detecting disparities in estimates of lung cancer mortality rates in 2014 between neighboring US counties after accounting for estimates of total (non-daily and daily) smoking prevalence in 2012.

Data References:
> A. H. Mokdad, L. Dwyer-Lindgren, C. Fitzmaurice, R. W. Stubbs, A. Bertozzi-Villa, C. Morozoff,
R. Charara, C. Allen, M. Naghavi, and C. J. L. Murray. Trends and Patterns of Disparities in
Cancer Mortality Among US Counties, 1980-2014. JAMA, 317(4):388–406, Jan. 2017. ISSN
0098-7484. doi: 10.1001/jama.2016.20324.

> L. Dwyer-Lindgren, A. H. Mokdad, T. Srebotnjak, A. D. Flaxman, G. M. Hansen, and C. J. Murray.
Cigarette smoking prevalence in US counties: 1996-2012. Population Health Metrics, 12(1):5,
Mar. 2014. ISSN 1478-7954. doi: 10.1186/1478-7954-12-5.

US county-level mortality rate estimates are available from the Institute of Health Metrics and Evaluation (IHME) at:
https://ghdx.healthdata.org/record/ihme-data/united-states-cancer-mortality-rates-county-1980-2014

US county-level smoking prevalence estimates are available from the IHME at:
https://ghdx.healthdata.org/record/ihme-data/united-states-smoking-prevalence-county-1996-2012 

### Roadmap of the Repository

| Directory | Contents/Description |
| --- | --- |
| *src* | contains source code to setup simulation and data analysis, implement sampling algorithms, and analyze difference boundaries |
| *data* | contains data used in Section 5 of the manuscript |
| *test* | contains simulation and data analysis `R` scripts for all results and figures in the manuscript |
| *output* | contains all results and figures from simulation and data analysis scripts output from scripts in *test*|

# Reproducing Figures

Due to the IHME's data use agreement, we cannot share the data as part of the repository. To reproduce the figures, please download the datasets using the URLs listed above and save the main data files (containing the relevant estimates) at the following locations:

- Cancer mortality estimate data: `data/US_data/IHME_county_cancer_mortality.xlsx`
- Smoking prevalence data: `data/US_data/IHME_data/IHME_US_COUNTY_TOTAL_AND_DAILY_SMOKING_PREVALENCE_1996_2012.csv`

We provide several R files that reproduce the results and figures in the simulation examples and real data analysis example in the manuscript. For the simulation and real data analysis files, we describe the data used, the model utilized, and the corresponding section in the manuscript where the figures appear. For a majority of the files, we list R Markdown files that may be more accessible but produce equivalent results to the R scripts. 

R Markdown files:

1. `test/Rmd/pc_prior_selection.Rmd`
   - Selection of hyperparameter $\lambda_{\rho}$ for penalized complexity prior on $\rho$ for simulation example and real data analysis
   - R script: `test/pc_prior_selection.R`
2. `test/Rmd/BYM2_exact_simulation.Rmd`
   - Simulation over map of US counties
   - Exact conjugate BYM2 model
   - Section 4.1, Section S3
   - R script: `test/BYM2_exact_simulation.R`, outputs to `output/US_exact_sample_sim`
3. `test/Rmd/BYM2_MCMCsampler_test.Rmd`
   - Simulation over map of US counties
   - BYM2 model with penalized complexity prior on $\rho$ 
   - Section 4.2
   - R script: `test/BYM2_MCMCsampler_test.R`, outputs to `output/US_gibbs_sample_sim`
4. `test/Rmd/BYM2_poisson.Rmd`
   - Simulation over map of CA counties
   - BYM2 Poisson-GLM model
   - Section 4.3
   - R script: `test/bym2_poisson.R` outputs to `output/poisson_sim`
5. `test/Rmd/US_lung_RDA.Rmd`
   - Lung cancer mortality rate estimates on US county level in 2014, data setup in `src/R/RDA/US_data_setup.R`
   - BYM2 model with penalized complexity prior on $\rho$ 
   - Section 5, Section S7
   - R script: `test/US_lung_RDA.R` outputs to  `output/RDA/US_data`
6. `test/CA_sharp_boundary_test.R`
   - Simulations over map of CA counties with comparison to boundary detection using an areally-reference Dirichlet process
   - Section S2
   - Results saved to `output/CA_sharp_boundary_sim`, results assessed in `test/CA_sharp_boundary_assessment.R`

