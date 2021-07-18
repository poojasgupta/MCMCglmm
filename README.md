# mcmcglmm

The code described here accompanies the paper titled ["Host phylogeny matters: Examining sources of variation in infection risk by blood parasites across a tropical montane bird community in India" published in Parasites and Vectors](https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-020-04404-8). In this paper, I used a bayesian approach implemented in the R-package ["MCMCglmm"](https://rdrr.io/github/jarrodhadfield/MCMCglmm/) for fitting a generalized linear mixed model to our dataset by accounting for phylogenetic relationships among host species. The goal was to assess the association between parasite infection risk and host ecological and morphometric traits. Models were run separately for *Plasmodium* and *Haemoproteus* parasites.

An example of the expected format for the code to work is given [here] (../master/IndSpGeoData_example.csv)
