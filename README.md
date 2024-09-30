# MCMCglmm
## Bayesian phylogenetic mixed models 

This repository describes the [R code](SppMCMCglmm_example_PG.R) that accompanies the paper titled [Host phylogeny matters: Examining sources of variation in infection risk by blood parasites across a tropical montane bird community in India](https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-020-04404-8) published in *Parasites and Vectors*. In this paper, I used a bayesian approach implemented in the R-package [MCMCglmm](https://rdrr.io/github/jarrodhadfield/MCMCglmm/) for fitting a generalized linear mixed model to our dataset by accounting for phylogenetic relationships among host species. The goal was to assess the association between parasite infection risk and host ecological and morphometric traits. Models were run separately for *Plasmodium* and *Haemoproteus* parasites. Host infection status was modelled as a binary response variable (0 for uninfected, 1 for infected) with a logit link, for *Plasmodium* and *Haemoproteus*, and different species- and individual-level ecological traits were used as predictor variables. To account for shared ancestry between host species, a variance-covariance matrix of phylogenetic distances between host species generated from the host phylogeny was fitted as a random effect.

For both *Plasmodium* and *Haemoproteus*, I started with a fully parameterized model including all predictors and then ran subsequent reduced models by excluding non-significant predictors, one at a time based on *P*-values. For more details on the dataset and MCMCglmm approach, please read the methods section of the published paper.

An example of the dataset (subset of the full dataset) used to run the code is given below
1. [Prevalence data](Prevalence_data_example.csv) - A dataframe containing parasite prevalence data and various ecological and morphometric traits tested.
2. [Phylogeny](species_tree.tre) - Host phylogeny file

## Fluctuating Asymmetry (FA) Calculation for Tarsus Measurements

The R script calculates fluctuating asymmetry (FA) for tarsus measurements in birds based on both left and right side measurements. The analysis is conducted using a linear mixed-effects model, controlling for individual variation, species, and trapping site differences.

### Methods:
Briefly, we first generate a variable `Side`, where values are set to `-1` for right-side tarsus measurements and `1` for left-side tarsus measurements. The fixed effects in the model include **Species**, **Trapping Site (Major.Loc)**, and **Side**. These fixed effects account for the average (intercept) values of the measurements for each species and site, while the **Side** variable estimates the level of **directional asymmetry** (DA) through the fixed intercept.

The random effects structure includes a random intercept for each individual, capturing deviations from the fixed intercept (i.e., individual-specific size differences). Additionally, a random slope for the **Side** variable within each individual is included. This random slope estimates the magnitude of **fluctuating asymmetry (FA)** for tarsus measurements within each individual, after controlling for other factors like individual-, species-, and site-specific differences in size and magnitude of directional asymmetry.

### Steps:
1. **Data Preparation:**
   - The dataset includes left (`Left.T1`, `Left.T2`, `Left.T3`) and right (`Right.T1`, `Right.T2`, `Right.T3`) tarsus measurements, along with species and trapping site data.
   - The data is reshaped from wide to long format, separating left and right tarsus measurements into separate rows with a `Side` variable indicating the measurement's side (left or right).

2. **Linear Mixed-Effects Model:**
   - A linear mixed-effects model is fitted using the `lme4` package.
   - **Fixed effects**: Species, trapping site (Major.Loc), and Side (to estimate directional asymmetry).
   - **Random effects**: A random intercept for each individual and a random slope for Side within each individual. The random slope of Side estimates fluctuating asymmetry (FA) for each individual.

3. **FA Calculation:**
   - The absolute value of the random slope for Side is extracted as the FA value for each individual, measuring the magnitude of asymmetry for tarsus measurements after controlling for species and site effects.
   - These FA values are merged back into the original dataset, adding an `FA` column to the original data.

4. **Output:**
   - The final dataset includes the original columns and a new `FA` column containing each individual's fluctuating asymmetry value for tarsus measurements.

### Requirements:
- **R Packages:** 
  - `lme4` for linear mixed-effects modeling
  - `dplyr` and `tidyr` for data manipulation

