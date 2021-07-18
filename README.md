# MCMCglmm
## Bayesian phylogenetic mixed models 

This repository describes the [R code](SppMCMCglmm_example_PG.R) that accompanies the paper titled [Host phylogeny matters: Examining sources of variation in infection risk by blood parasites across a tropical montane bird community in India](https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-020-04404-8) published in *Parasites and Vectors*. In this paper, I used a bayesian approach implemented in the R-package [MCMCglmm](https://rdrr.io/github/jarrodhadfield/MCMCglmm/) for fitting a generalized linear mixed model to our dataset by accounting for phylogenetic relationships among host species. The goal was to assess the association between parasite infection risk and host ecological and morphometric traits. Models were run separately for *Plasmodium* and *Haemoproteus* parasites. Host infection status was modelled as a binary response variable (0 for uninfected, 1 for infected) with a logit link, for *Plasmodium* and *Haemoproteus*, and different species- and individual-level ecological traits were used as predictor variables. To account for shared ancestry between host species, a variance-covariance matrix of phylogenetic distances between host species generated from the host phylogeny was fitted as a random effect.

For both *Plasmodium* and *Haemoproteus*, I started with a fully parameterized model including all predictors and then ran subsequent reduced models by excluding non-significant predictors, one at a time based on *P*-values. For more details on the MCMCglmm approach, please read the methods section of the published paper.

An example of the dataset (subset of the full dataset) used to run the code is given below
1. [Prevalence data](IndSpGeoData_example.csv) - A dataframe containing parasite prevalence data and various ecological and morphometric traits tested.
2. [Phylogeny](sp_ultraTree.tre) - Host phylogeny file
