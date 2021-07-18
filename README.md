# mcmcglmm

The code described here accompanies the paper titled ["Host phylogeny matters: Examining sources of variation in infection risk by blood parasites across a tropical montane bird community in India" published in Parasites and Vectors](https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-020-04404-8). In this paper, I used a bayesian approach implemented in the R-package ["MCMCglmm"](https://rdrr.io/github/jarrodhadfield/MCMCglmm/) for fitting a generalized linear mixed model to our dataset by accounting for phylogenetic relationships among host species. The goal was to assess the association between parasite infection risk and host ecological and morphometric traits. Models were run separately for *Plasmodium* and *Haemoproteus* parasites.

An example of the dataset used to run the code is given 
* [Prevalence data](IndSpGeoData_example.csv) - A dataframe containing parasite prevalence data and various ecological and morphometric traits tested.
* [Phylogeny](sh-files/run_bwa_mem.sh) - Host phylogeny file 
