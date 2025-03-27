---
title: "Stock & Miller (2021)."
date: "Mar 27, 2025"
output:
  ioslides_presentation:
    smaller: no
    widescreen: true
    transition: "faster" 
---

## Impact

* Cited 65 times. 

* Many citing papers have similar methodological aims. Some citing papers are data analysis.

* At the NOAA research center in Woods Hole, MA.

## TMB vs particle filtering

* "we note that the development of TMB has been a critical advance for fisheries assessment modeling frameworks such as WHAM, allowing us to rapidly fit models that treat population and environmental processes as time-varying random effects in a state-space framework."

* Is the magic of TMB due to the strength of the Laplace approximation, the use of autodiff, other software quality issues, or something else?


## Mohn's rho

[Mohn, R. (1999)](https://doi.org/10.1006/jmsc.1999.0481) The retrospective problem in sequential population analysis: An investigation using cod fishery and simulated data. ICES Journal of Marine Science, 56, 473–488.

* Apparently, model misspecification can lead to widespread incongruous results.

* E.g., failure to describe increases in skill at catching.

## Model comparison

* Appendix B deals with model specification. It appears to have overdispersion only in the measurement model.

* 2.1.2.2. Catch and index age composition. This explains why the log-normal is preferred for the measurment model, in order to be "self-weighting" and allow for correlations.

* Continuing work on evaluation and comparison of models: https://doi.org/10.1016/j.fishres.2024.106968

## Reproducibility

* "Documentation and tutorials for how to specify additional random effect structures in WHAM are available at https://timjmiller.github.io/wham/."

* "Code and data files to run the analysis presented here are available at https://github.com/brianstock-NOAA/wham-sim."

* The models and data are well specified, and within the reach of an (ambitious) 531/631 final project








