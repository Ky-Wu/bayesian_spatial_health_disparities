# Bayesian Assessment of Spatial Health Disparities Using Disease Maps

This repository contains code implementation for detecting spatial disparities in a BYM2 model using the methodology discussed in the reference below.  

Reference:
> Wu, K. L., and Banerjee, S., "Bayesian Asessment of Spatial Health Disparities Using Disease Maps." Manuscript in preparation.

A real data analysis example is included to showcase detecting disparities in estimates of lung cancer mortality rates in 2014 between neighboring US counties after accounting for smoking prevalence.

References:
> A. H. Mokdad, L. Dwyer-Lindgren, C. Fitzmaurice, R. W. Stubbs, A. Bertozzi-Villa, C. Morozoff,
R. Charara, C. Allen, M. Naghavi, and C. J. L. Murray. Trends and Patterns of Disparities in
Cancer Mortality Among US Counties, 1980-2014. JAMA, 317(4):388–406, Jan. 2017. ISSN
0098-7484. doi: 10.1001/jama.2016.20324.

> L. Dwyer-Lindgren, A. H. Mokdad, T. Srebotnjak, A. D. Flaxman, G. M. Hansen, and C. J. Murray.
Cigarette smoking prevalence in US counties: 1996-2012. Population Health Metrics, 12(1):5,
Mar. 2014. ISSN 1478-7954. doi: 10.1186/1478-7954-12-5.

# Reproducing Figures

We provide several R Markdown files that reproduce the results and figures in the simulation examples and real data analysis example in the manuscript. For the simulation and real data analysis files, we describe the data used, the type of model utilized, and the corresponding section in the manuscript where the figures appear.

1. `test/Rmd/pc_prior_selection.Rmd`
   -Selection of hyperparameter $\lamba_\rho$ for penalized complexity prior on $\rho$ for simulation example and real data analysis
2. `test/Rmd/BYM2_exact_simulation.Rmd`
   -Simulation over map of US counties
   -Detecting spatial disparities under an exact conjugate model
   -Section 4.1
3. `test/Rmd/BYM2_MCMCsampler_test.Rmd`
   -Simulation over map of US counties
   -Detecting spatial disparities under a BYM2 model with penalized complexity prior on $\rho$ 
   -Section 4.2
4. `test/Rmd/US_lung_RDA.Rmd`
   -Lung cancer mortality rate estimates on US county level in 2014
   -Detecting spatial disparities under a BYM2 model with penalized complexity prior on $\rho$ 
   -Section 5

