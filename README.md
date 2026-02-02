# [README.md](https://github.com/user-attachments/files/24998499/README.md)

Modelling the Impact of Wolbachia Replacement Strategies on Dengue Risk

**Overview**

This repository contains all scripts used to generate the results presented in the modelling study evaluating the impact of Wolbachia replacement strategies on dengue risk under varying baseline mosquito population sizes and Wolbachia blocking effects.

Baseline population size is defined as the initial wild-type mosquito population relative to the equilibrium carrying capacity. The study compares continuous release and threshold-based (80% goal) release strategies.

The modelling framework is a discrete-time, discrete-space population model, adapted from:
Turelli, M. (2020). Cytoplasmic incompatibility in populations with overlapping generations.https://doi.org/10.1111/j.1558-5646.2009.00822.x

All models and analyses are implemented in R.



**Model Description**

Model type: Discrete-time, discrete-space population model

Framework: Adapted from Turelli (2020)
Key components include wild-type and Wolbachia-infected mosquito populations, density-dependent population dynamics, cytoplasmic incompatibility, Wolbachia-induced blocking effects on dengue transmission


Strategies evaluated includes continuous release strategy, and threshold (80% goal) release strategy



**Core models**

turelli2020\_vc\_func\_20251103.R: Defines the core model for the continuous release strategy.
turelli2020\_vc\_80percentgoal\_func\_20251103.R: Defines the core model for the threshold (80% goal) release strategy

**Temporal**

temporal\_simulation\_20251128\_clean.R: Generates temporal simulation results for 10%, 20%, and 90% baseline population ratios under the continuous release strategy.
temporal\_simulation\_80goal\_20251128\_clean.R: Generates corresponding temporal simulation results under the threshold release strategy.

**Parameter sensitivity** 

varying\_para\_20251204\_clean.R: Explores simulation outcomes across varying baseline population ratios and Wolbachia blocking effects for the continuous release strategy.
varying\_para\_80goal\_20251204\_clean.R: Performs the same analyses for the threshold release strategy.

**Spatial**

turelli2020\_space\_runc\_20251208.R:Defines the spatial model framework.
spatial\_simulation\_20260130.R: Runs spatial simulations based on the defined framework.

**Figures**

results\_plot\_20260130\_clean.R:Compiles simulation outputs and generates all summary tables and figures presented in the study.





