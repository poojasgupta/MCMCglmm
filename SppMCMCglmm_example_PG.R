library(MCMCglmm); library(ape)
library(lattice); library(plotMCMC); library(scape)
library(coda); library(ggplot2); library(dplyr)
library(MuMIn); library(gdata)
library(RColorBrewer);library(png);library(TeachingDemos)
library(xlsx); library(broom.mixed)
#Run MCMCglmm_functions.R before running this code.

resMCMC.List=list(P=list(),H=list())

##################################################################################
##Get input data
##################################################################################
mData <- read.table("IndSpGeoData_example.csv", header=TRUE, sep=",")

mData$Inf_Haem <- as.factor(mData$Inf_Haem)
mData$Inf_Plas <- as.factor(mData$Inf_Plas)

GEN.COLS = c("PhyloCode","Species","Lab.ID","Major.Loc","Minor.Loc",
             "Extracted", "Infected.N","G_Lineage.name","Inf_Plas","Inf_Haem","Phylo.Lineage.Name")

ECO.COLS=c("ECO_feedStrata","ECO_socialRoost","ECO_sdm",
           "ECO_habitat","ECO_connect", "ECO_SppWt.div10","ECO_minElev.div100")

IND.COLS=c("IND_Cen.BodyCI2.div10") 

##---------------------------------------------------------------------------------
#Getting host species tree
spTree <- read.tree("sp_ultraTree.tre")
# Transforming species phylo tree into a usable ginverse matrix
Ainv.m1 <-inverseA(spTree, nodes="ALL")$Ainv

#parasite infection data matrix
pData <- mData[,c(GEN.COLS,ECO.COLS,IND.COLS)]

####################################################################
#MCMCglmm priors and run settings
#------------------------------------------------------------------
NITT=2e+06;BITT=100000;THIN=1000

prior1 = list(R = list(V = 10,fix=1),
              G=list(G1=list(V=1, nu=1000, alpha.mu=0, alpha.V=1),
                     G2=list(V=1, nu=1000, alpha.mu=0, alpha.V=1)))

####################################################################
####################################################################
##----------------------------------------------------------------
## Plasmodium
##----------------------------------------------------------------
####################################################################
Pmod.List=list()
##---------------------------------------
## FULL MODEL
##---------------------------------------
Analysis="P.fullspp"
Pmod1.full <- P.mcmcrun(pData, ECO.COLS)
Pmod.List[[Analysis]]=Pmod1.full
Pmod1.full.res <- summary(Pmod1.full)
Pmod1.full.res.df <- tidy(Pmod1.full, effects = c("fixed", "ran_pars"), ess=T, conf.int = T,conf.level = 0.95,
                          conf.method = "HPDinterval")

##--------------------------------------------------------------------------
##Reduced models by removing non-sginificant variables one at a time, p>0.1
##-------------------------------------------------------------------------- 

Analysis="P.redspp1"  #Dropping ECO_minElev.div100"
newVars1=c("ECO_feedStrata","ECO_socialRoost","ECO_sdm","ECO_habitat","ECO_connect","ECO_SppWt.div10") 
Pmod1.red1 <- P.mcmcrun(pData,newVars1)
  
Pmod.List[[Analysis]]=Pmod1.red1
Pmod1.red1.res <- summary(Pmod1.red1)
Pmod1.red1.res.df <- tidy(Pmod1.red1, effects = c("fixed", "ran_pars"), ess=T, conf.int = T,conf.level = 0.95,
                          conf.method = "HPDinterval")

#-------------------------------------
Analysis="P.redspp2"  #Dropping Spp.Wt
newVars2=c("ECO_feedStrata","ECO_socialRoost","ECO_sdm","ECO_habitat","ECO_connect")
Pmod1.red2 <- P.mcmcrun(pData,newVars2)
Pmod.List[[Analysis]]=Pmod1.red2
Pmod1.red2.res <- summary(Pmod1.red2) ####
Pmod1.red2.res.df <- tidy(Pmod1.red2, effects = c("fixed", "ran_pars"), ess=T, conf.int = T,conf.level = 0.95,
                          conf.method = "HPDinterval")

#--------------------------------------
Analysis="P.redspp3"  #Drop connectivity
newVars3=c("ECO_feedStrata","ECO_socialRoost","ECO_sdm","ECO_habitat")
Pmod1.red3 <- P.mcmcrun(pData,newVars3)
Pmod.List[[Analysis]]=Pmod1.red3
Pmod1.red3.res <- summary(Pmod1.red3)
Pmod1.red3.res.df <- tidy(Pmod1.red3, effects = c("fixed", "ran_pars"), ess=T, conf.int = T,conf.level = 0.95,
                          conf.method = "HPDinterval")

##-------------------------------------
Analysis="P.redspp4"  #Dropping habitat
newVars4=c("ECO_feedStrata","ECO_socialRoost","ECO_sdm")
Pmod1.red4 <- P.mcmcrun(pData,newVars4)
Pmod.List[[Analysis]]=Pmod1.red4
Pmod1.red4.res <- summary(Pmod1.red4)
Pmod1.red4.res.df <- tidy(Pmod1.red4, effects = c("fixed", "ran_pars"), ess=T, conf.int = T,conf.level = 0.95,
                          conf.method = "HPDinterval")


###-----------------------------------------------------------------------
##ECO_feedStrata,ECO_socialRoost,ECO_sdm are the significant variables and will be used for multiple runs
##------------------------------------------------------------------------
#MCMC RUN2 ##Run second chain to check for convergence and increase confidence in estimates
##------------------------------------------------------------------------
set.seed(3212)
Analysis="P.redspp4.run2"
Pmod2.red4 <- P.mcmcrun(pData,newVars4)
Pmod.List[[Analysis]]=Pmod2.red4
Pmod2.red4.res <- summary(Pmod2.red4) #MCMC RUN2
Pmod2.red4.res.df <- tidy(Pmod2.red4, effects = c("fixed", "ran_pars"), ess=T, conf.int = T,conf.level = 0.95,
                          conf.method = "HPDinterval")

##--------------------------------------------------------------------------
#MCMC RUN3 ##Run third chain to check for convergence and increase confidence in estimates
##--------------------------------------------------------------------------
set.seed(5839)
Analysis="P.redspp4.run3"
Pmod3.red4 <- P.mcmcrun(pData,newVars4)
Pmod.List[[Analysis]]=Pmod3.red4
Pmod3.red4.res <- summary(Pmod3.red4) #MCMC RUN3
Pmod3.red4.res.df <- tidy(Pmod3.red4, effects = c("fixed", "ran_pars"), ess=T, conf.int = T,conf.level = 0.95,
                          conf.method = "HPDinterval")

##-------------------------------------------------------------
##Summarizing results from all models in a list and a dataframe
##-------------------------------------------------------------
resMCMC.List[["P"]]=Pmod.List 

Pmod.all.res.df <- bind_rows(Pmod1.full.res.df, Pmod1.red1.res.df, Pmod1.red2.res.df, Pmod1.red3.res.df, 
                             Pmod1.red4.res.df, Pmod2.red4.res.df, Pmod3.red4.res.df, .id="Model.type")
write.csv(Pmod.all.res.df, file ="Pmod_allRes_df.csv")

#For running clean.mcmc and getting pvalues
Pmod.all.clean <- list(Pmod1.full, Pmod1.red1, Pmod1.red2, Pmod1.red3, 
                             Pmod1.red4, Pmod2.red4, Pmod3.red4)

###############################################################
##Post fit Analysis
###############################################################
#Merging chains of both MCMC runs
P.mcmcl.Sol <- mcmc.list(Pmod1.red4$Sol, Pmod2.red4$Sol, Pmod3.red4$Sol)
P.mcmcl.vcv <- mcmc.list(log10(Pmod1.red4$VCV), log10(Pmod2.red4$VCV), log10(Pmod3.red4$VCV))

#Plotting merged trace plots to visually check for stationarity and convergence
par(mfrow=c(8,2), mar=c(2, 2, 2, 1)) #size 643png
xyplot(P.mcmcl.Sol,  main="Plasmodium Trace Plot for Fixed Variables", col=c("red", "green", "blue"))
xyplot(P.mcmcl.vcv, main="Plasmodium Trace Plot for Random Variables", col=c("red", "green", "blue"))
dev.off()

##------------------------------------------------------
##COnvergence diagnostics
##------------------------------------------------------
Pdiag.List=list()

##1)Get effective sample sizes for variables. EFS >1,000 are good 
Analysis="ESS"
ESS.Sol <- as.data.frame(effectiveSize(P.mcmcl.Sol))
ESS.vcv <- as.data.frame(effectiveSize(P.mcmcl.vcv))
ESS <- list(ESS.Sol, ESS.vcv)
Pdiag.List[[Analysis]]=ESS

##2) Autocorrelation - correlation between samples at different lags. 
# all correlations should be less (in magnitude) than 0.1
Analysis="Autocorr"
Autocorr <- as.data.frame(list(autocorr.diag(P.mcmcl.Sol), autocorr.diag(P.mcmcl.vcv)))
Pdiag.List[[Analysis]]=Autocorr

autocorr.plot(P.mcmcl.Sol)
autocorr.plot(P.mcmcl.vcv)

#3)Potential scale reduction factors and a Gelman and Rubin diagnostic for
# comparing multiple chains. Values near 1 are good
Analysis="Gelman"
gelman <- gelman.diag(P.mcmcl.Sol, confidence = 0.95,transform = T)
Pdiag.List[[Analysis]]=gelman

par(mfrow=c(2,2), mar=c(2,2,1,2))
gelman.plot(P.mcmcl.Sol, auto.layout = F, autoburnin = T)
            
#4) Heidelberg and Welch diagnostic tests a null hypoTion
Analysis="Heidel"
heidel <- heidel.diag(P.mcmcl.Sol) #Passed
Pdiag.List[[Analysis]]=heidel

heidel.diag(P.mcmcl.vcv) #Passed

## From the individual level analysis (code not included here), there were no significant variables for P
##################################################################
##----------------------------------------------------------------
## Haemoprotues
##----------------------------------------------------------------
##################################################################
Hmod.List=list()
##---------------------------------------
## FULL MODEL
##---------------------------------------
Analysis="H.fullspp"
Hmod1.full <- H.mcmcrun(pData, ECO.COLS)
Hmod.List[[Analysis]]=Hmod1.full
Hmod1.full.res <- summary(Hmod1.full)
Hmod1.full.res.df <- tidy(Hmod1.full, effects = c("fixed", "ran_pars"), ess=T, conf.int = T, conf.level = 0.95,
                          conf.method = "HPDinterval")

##-----------------------------------------------------------------------------------
## REDUCED MODELS- Drop all non-significant variables, p>0.01 
##------------------------------------------------------------------------------------
Analysis="H.redspp1" #dropped SppWt
newVars1=c("ECO_feedStrata","ECO_socialRoost","ECO_sdm","ECO_habitat","ECO_connect","ECO_minElev.div100") #for Hmod1.red1
Hmod1.red1 <- H.mcmcrun(pData, newVars1)
Hmod.List[[Analysis]]=Hmod1.red1
Hmod1.red1.res <- summary(Hmod1.red1) #MCMC RUN1
Hmod1.red1.res.df <- tidy(Hmod1.red1, effects = c("fixed", "ran_pars"), ess=T, conf.int = T,conf.level = 0.95,
                          conf.method = "HPDinterval")

#--------------------------------------
Analysis="H.redspp2"  #Dropped sdm
newVars2=c("ECO_feedStrata","ECO_socialRoost","ECO_SppWt.div10","ECO_habitat","ECO_minElev.div100") #for Hmod1.red2
Hmod1.red2 <- H.mcmcrun(pData, newVars2)
Hmod.List[[Analysis]]=Hmod1.red2
Hmod1.red2.res <- summary(Hmod1.red2) #MCMC RUN1
Hmod1.red2.res.df <- tidy(Hmod1.red2, effects = c("fixed", "ran_pars"), ess=T, conf.int = T,conf.level = 0.95,
                          conf.method = "HPDinterval")

#--------------------------------------
Analysis="H.redspp3" #dropped habitat
newVars3=c("ECO_feedStrata","ECO_socialRoost","ECO_SppWt.div10","ECO_minElev.div100") #for Hmod1.red3
Hmod1.red3 <- H.mcmcrun(pData, newVars3)
Hmod.List[[Analysis]]=Hmod1.red3
Hmod1.red3.res <- summary(Hmod1.red3)
Hmod1.red3.res.df <- tidy(Hmod1.red3, effects = c("fixed", "ran_pars"), ess=T, conf.int = T,conf.level = 0.95,
                          conf.method = "HPDinterval")

#----------------------------------------
Analysis="H.redspp3b" ##dropped feeding strata
newVars3b=c("ECO_socialRoost","ECO_SppWt.div10","ECO_minElev.div100") #for Hmod1.red2 #for Hmod1.red3
Hmod1.red3b <- H.mcmcrun(pData, newVars3b)
Hmod.List[[Analysis]]=Hmod1.red3b
Hmod1.red3b.res <- summary(Hmod1.red3b)
Hmod1.red3b.res.df <- tidy(Hmod1.red3b, effects = c("fixed", "ran_pars"), ess=T, conf.int = T,conf.level = 0.95,
                           conf.method = "HPDinterval")

#----------------------------------------
Analysis="H.redspp4" #dropped SppWt
newVars4=c("ECO_socialRoost","ECO_minElev.div100") #for Hmod1.red4
Hmod1.red4 <- H.mcmcrun(pData, newVars4)
Hmod.List[[Analysis]]=Hmod1.red4
Hmod1.red4.res <- summary(Hmod1.red4)
Hmod1.red4.res.df <- tidy(Hmod1.red4, effects = c("fixed", "ran_pars"), ess=T, conf.int = T,conf.level = 0.95,
                          conf.method = "HPDinterval")

##------------------------------------------------------------------------------
## ECO_socialRoost, ECO_minElev.Num, are significant in the species level model
## From the individual level mcmcglmm analysis (code not shown here), IND_Cen.BodyCI2.div10 was significant.
## For the final Haemoproteus runs, merging species and Individual models and running
## 3 runs for model containing ECO_socialRoost, ECO_minElev.Num and IND_Cen.BodyCI2.div10
##------------------------------------------------------------------------------
						 
 
 H.combVars2=c("ECO_socialRoost","ECO_minElev.div100","IND_Cen.BodyCI2.div10")
 SppIndList.Res=list()
 
##------------------------------------------------------------------------------
#MCMC RUN1 ##Run first chain to check for convergence and increase confidence in estimates
##------------------------------------------------------------------------------
  Analysis="H.SppInd.run1"
  set.seed(9283)
  H.SppInd1 <- H.mcmcrun(pData, H.combVars2)
  H.SppInd1.res <- summary(H.SppInd1)
  SppIndList.Res[[Analysis]] <- H.SppInd1
  H.SppInd1.df <- tidy(H.SppInd1, effects = c("fixed", "ran_pars"), ess=T, conf.int = T,conf.level = 0.95,
                          conf.method = "HPDinterval")
##------------------------------------------------------------------------------
#MCMC RUN2 ##Run second chain to check for convergence and increase confidence in estimates
##------------------------------------------------------------------------------ 
  Analysis="H.SppInd.run2"
  set.seed(8372)
  H.SppInd2 <- H.mcmcrun(pData, H.combVars2)
  H.SppInd2.res <- summary(H.SppInd2)
  SppIndList.Res[[Analysis]] <- H.SppInd2
  H.SppInd2.df <- tidy(H.SppInd2, effects = c("fixed", "ran_pars"), ess=T, conf.int = T,conf.level = 0.95,
                          conf.method = "HPDinterval")
						  
##------------------------------------------------------------------------------
#MCMC RUN3 ##Run third chain to check for convergence and increase confidence in estimates
##------------------------------------------------------------------------------
  Analysis="H.SppInd.run3 "
  set.seed(37289)
  H.SppInd3 <- H.mcmcrun(pData, H.combVars2)
  H.SppInd3.res <- summary(H.SppInd3)
  SppIndList.Res[[Analysis]] <- H.SppInd3
  H.SppInd3.df <- tidy(H.SppInd3, effects = c("fixed", "ran_pars"), ess=T, conf.int = T,conf.level = 0.95,
                          conf.method = "HPDinterval")

##########################################################
##------------------------------------------------------
##Combine all model results in a list and a dataframe
##------------------------------------------------------
resMCMC.List[["H"]]=Hmod.List 

Hmod.all.res.df <- bind_rows(Hmod1.full.res.df, Hmod1.red1.res.df, Hmod1.red2.res.df, Hmod1.red3.res.df, Hmod1.red3b.res.df, 
                             Hmod1.red4.res.df,H.SppInd1.df ,H.SppInd2.df, H.SppInd3.df ,.id="Model.type")
write.csv(Hmod.all.res.df, file ="Hmod_allRes_df.csv")

#For running clean.mcmc and getting pvalues
Hmod.all.clean <- list(Hmod1.full,Hmod1.red1,Hmod1.red2,Hmod1.red3,Hmod1.red3b,
                       Hmod1.red4,H.SppInd1,H.SppInd2,H.SppInd3)

###############################################################
  ##Post fit Analysis
  ###############################################################
  #Merging chains of both MCMC runs
  H.mcmcl.Sol <- mcmc.list(H.SppInd1$Sol, H.SppInd2$Sol, H.SppInd3$Sol)
  H.mcmcl.vcv <- mcmc.list(log10(H.SppInd1$VCV), log10(H.SppInd2$VCV), log10(H.SppInd3$VCV))
  
  
  #Plotting merged trace plots to visually check for stationarity and convergence
  par(mfrow=c(4,2), mar=c(2, 2, 2, 1)) #size 643png
  xyplot(H.mcmcl.Sol, main="Haemoproteus Trace Plots for Fixed Variables", col=c("red", "green", "blue")); 
  xyplot(H.mcmcl.vcv, main="Haemoproteus Trace Plots for Random Variables", col=c("red", "green", "blue"))
  
  ##########################
  #COnvergence diagnostics
  ##########################
  H_SppInd_diag.List=list()
  
  ##1)Get effective sample sizes for variables. EFS >1,000 are good
  Analysis="ESS"
  ESS.Sol <- effectiveSize(H.mcmcl.Sol)
  ESS.vcv <- effectiveSize(H.mcmcl.vcv)
  
  ESS <- list(ESS.Sol, ESS.vcv)
  H_SppInd_diag.List[[Analysis]]=ESS
  
  ##2) Autocorrelation - correlation between samples at different lags. 
  # all correlations should be less (in magnitude) than 0.1
  Analysis="Autocorr"
  Autocorr <- as.data.frame(list(autocorr.diag(H.mcmcl.Sol), autocorr.diag(H.mcmcl.vcv)))
  H_SppInd_diag.List[[Analysis]]=Autocorr
  
  autocorr.plot(H.mcmcl.Sol)
  autocorr.plot(H.mcmcl.vcv)
  
  #3)Potential scale reduction factors and a Gelman and Rubin diagnostic for
  # comparing multiple chains. Values near 1 are good
  Analysis="Gelman"
  gelman <- gelman.diag(H.mcmcl.Sol, confidence = 0.95, transform = T)
  H_SppInd_diag.List[[Analysis]]=gelman
  
  par(mfrow=c(2,2), mar=c(2,2,1,2))
  gelman.plot(H.mcmcl.Sol, auto.layout = F, autoburnin = T, transform = T)
  
  #4) Heidelberg and Welch diagnostic tests a null hypothesis that 
  # the Markov chain is from a stationary distribution
  Analysis="Heidel"
  heidel <- heidel.diag(H.mcmcl.Sol) #Passed
  H_SppInd_diag.List[[Analysis]]=heidel

##################################################################
##################################################################
# Gather regression coefficients, HPD intervals and plot
##################################################################
##################################################################
#----------------------------------------------------------------
#Plasmodium
#----------------------------------------------------------------
vn <- varnames(P.mcmcl.Sol)  
P.mode <- posterior.mode(as.mcmc(rbind(Pmod1.red4$Sol, Pmod2.red4$Sol, Pmod3.red4$Sol)))
P.intervals <-HPDinterval(as.mcmc(rbind(Pmod1.red4$Sol, Pmod2.red4$Sol, Pmod3.red4$Sol)), prob=0.95)
P.estimates <-cbind(data.frame(Predictors= vn, Posterior.Mode= P.mode, P.intervals))
##----------------------
#Calculate effect sizes- exponentiate coefficients to back-transform from log link
##----------------------
P.exp.estimates= P.estimates[,2:4] <- exp(P.estimates[,2:4])
write.csv(P.exp.estimates, file="Plasm_expEstimates.csv")

#-----------------------------------------------------------------
#Haemoproteus
#-----------------------------------------------------------------
vn <- varnames(H.SppInd1$Sol)
  H.mode<- posterior.mode(as.mcmc(rbind(H.SppInd1$Sol, H.SppInd2$Sol,H.SppInd3$Sol)))
  H.intervals<-HPDinterval(as.mcmc(rbind(H.SppInd1$Sol, H.SppInd2$Sol,H.SppInd3$Sol)), prob=0.95)
  H.SppInd.estimates<-cbind(data.frame("Predictors"=vn, "Posterior.Mode"=H.mode, H.intervals))
  
 ##----------------------
 #Calculate effect sizes- exponentiate coefficients to back-transform from log link
 #-----------------------
  H.SppInd.exp.estimates= H.SppInd.estimates[,2:4] <- exp(H.SppInd.estimates[,2:4])
  write.csv(H.SppInd.exp.estimates, file="SppInd_expEstimates.csv")


##################################################################
###################################################################
#Getting MCMC model all outputs for P and H together
##################################################################
##################################################################
##Plasmodium 
##-------------------------
ModListNames <- list("Pmod1.full","Pmod1.red1","Pmod1.red2","Pmod1.red3","Pmod1.red4",
                     "Pmod2.red4","Pmod3.red4")
##Get the clean.MCMC outputs and add modelName columns to each element for ID purposes:
readyListP <- mapply(cbind,"modelType" = ModListNames,lapply(Pmod.all.clean, clean.MCMC),SIMPLIFY = F)

##Turn the list of data.frames into one big data.frame:
#P.mcmcOutputs <- as.data.frame(do.call(rbind, readyListP), stringsAsFactors = FALSE)

##--------------------------
#Haemoproteus
ModListNames <- list("Hmod1.full","Hmod1.red1","Hmod1.red2","Hmod1.red3","Hmod1.red3b",
                     "Hmod1.red4","H.SppInd1","H.SppInd2","H.SppInd3")
##Get the clean.MCMC outputs and add modelName columns to each element for ID purposes:
readyListH <- mapply(cbind,"modelType" = ModListNames, lapply(Hmod.all.clean, clean.MCMC), SIMPLIFY = F)

#------------------------------------------------------------------------
##Turn the list of data.frames into one big data.frame including P and H
mcmcOutputs <- bind_rows(readyListP, readyListH, .id="ParSpp")
write.xlsx(mcmcOutputs, file="MCMCoutput_all_final.xlsx")

#Getting an alternate format of MCMC results obtained by tidy
#MCMCoutput_altFormat <- bind_rows(Pmod.all.res.df, Hmod.all.res.df, .id="ParSpp")
#MCMCoutput_altFormat$ParSpp[1:59] <- "Plasm"
#MCMCoutput_altFormat$ParSpp[60:118] <- "Haem"
#write.csv(MCMCoutput_altFormat, file="MCMCoutput_all_AltFormat.csv")

##-----------------------------------------------------------------------
##Merging and writing diagnostics output for P and H
MCMCdiag <- list(P=Pdiag.List, H=Hdiag.List)
capture.output(MCMCdiag, file = "Spp_MCMCdiagnostics.txt")
##-------------------------------------------------------------------------

###########################################################################
##MCMCglmm Analysis COMPLETE
###########################################################################
###########################################################################
#Estimating phylogenetic variance
###########################################################################
## Posterior probability of the phylogenetic signal lambda(just like estimating heritability)
## this is basically estimating variation explained by each variance components
#-------------------------------------------------------------------------------
#---------------------
#Plasmodium full model
#---------------------
lambda.sp <- (Pmod1.full$VCV[,1]/(Pmod1.full$VCV[,1]+Pmod1.full$VCV[,2] + Pmod1.full$VCV[,3]+1))
lambda.loc <- (Pmod1.full$VCV[,2]/(Pmod1.full$VCV[,1]+Pmod1.full$VCV[,2] + Pmod1.full$VCV[,3]+1))

P.full.lambda.sp <- as.data.frame(cbind(posterior.mode(lambda.sp),HPDinterval(lambda.sp)))
P.full.lambda.loc <- as.data.frame(cbind(posterior.mode(lambda.loc),HPDinterval(lambda.loc)))

##-----------------------------
#Plasmodium reduced final model
#------------------------------
lambda.sp <- (Pmod1.red4$VCV[,1]/(Pmod1.red4$VCV[,1]+Pmod1.red4$VCV[,2] + Pmod1.red4$VCV[,3]+1))
lambda.loc <- (Pmod1.red4$VCV[,2]/(Pmod1.red4$VCV[,1]+Pmod1.red4$VCV[,2] + Pmod1.red4$VCV[,3]+1))

P.red.lambda.sp <- as.data.frame(cbind(posterior.mode(lambda.sp),HPDinterval(lambda.sp)))
P.red.lambda.loc <- as.data.frame(cbind(posterior.mode(lambda.loc),HPDinterval(lambda.loc)))

##-----------------------------------------------------------------------------
##Haemoproteus
##-----------------------------------------------------------------------------
lambda.sp <- (Hmod1.full$VCV[,1]/(Hmod1.full$VCV[,1]+Hmod1.full$VCV[,2] + Hmod1.full$VCV[,3]+1))
lambda.loc <- (Hmod1.full$VCV[,2]/(Hmod1.full$VCV[,1]+Hmod1.full$VCV[,2] + Hmod1.full$VCV[,3]+1))

H.full.lambda.sp <- as.data.frame(cbind(posterior.mode(lambda.sp),HPDinterval(lambda.sp)))
H.full.lambda.loc <- as.data.frame(cbind(posterior.mode(lambda.loc),HPDinterval(lambda.loc)))


lambda.sp <- (H.SppInd1$VCV[,1]/(H.SppInd1$VCV[,1]+H.SppInd1$VCV[,2] + H.SppInd1$VCV[,3]+1))
lambda.loc <- (H.SppInd1$VCV[,2]/(H.SppInd1$VCV[,1]+H.SppInd1$VCV[,2] + H.SppInd1$VCV[,3]+1))
  
H.red.lambda.sp <- as.data.frame(cbind(posterior.mode(lambda.sp),HPDinterval(lambda.sp)))
H.red.lambda.loc <- as.data.frame(cbind(posterior.mode(lambda.loc),HPDinterval(lambda.loc)))
 

#Summarizing phylogenetic variance output for P and H

postmode.var <- bind_rows(P.full.lambda.sp, P.full.lambda.loc,P.red.lambda.sp, P.red.lambda.loc,
                          H.full.lambda.sp, H.full.lambda.loc,H.red.lambda.sp, H.red.lambda.loc, .id="ID")
rownames(postmode.var) <- c("P.full.lambda.sp", "P.full.lambda.loc","P.red.lambda.sp", "P.red.lambda.loc",
                            "H.full.lambda.sp", "H.full.lambda.loc", "H.red.lambda.sp", "H.red.lambda.loc")
write.csv(postmode.var, file="PhyloVariance_estimates.csv")

#####################################################################
#####################################################################
## PLOT MCMCGLMM ANALYSIS RESULTS
#####################################################################
#####################################################################
#Coeffecients Plot
library(ggplot2)
library(gridExtra); library(grid);library(ggpubr)
library(ggthemes); library(scales)

##################################################################
##Publication quality Plots
##################################################################
##Function------------------------------------------------------
##---------------------------------------------------------------
#Code Source: https://rpubs.com/Koundy/71792
##----------------------------------------------------------------
## Axis title could be made bold [axis.title = element_text(size = rel(1.1))]
## margin(t, r, b, l);left margin 1.5 for Haem plot
################################################################
theme_Publication <- function(base_size=13, base_family="Helvetica") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold.italic", size = rel(1.1)),
           text = element_text(),
           #panel.background = element_rect(colour = NA),
           panel.background =element_rect(fill = "grey92", colour = NA),
           plot.background = element_rect(colour = NA),
           #plot.margin=unit(c(3,3,3,3),"mm"),
           plot.margin=unit(c(3,1,2,3),"mm"), #Plas
           #plot.margin=unit(c(3,3,2,1),"mm"),  #Haem
           #panel.border = element_rect(colour = NA),
           #axis.title = element_text(size = rel(1.1)),
           #axis.title = element_text(size = 16), #used for MCMCglmm plots
           axis.title = element_text(size = 14), #used for sp prediction prob.
           axis.title.y = element_text(angle=90),
           axis.title.x = element_text(vjust = -0.3),
           axis.text = element_text(), 
           #axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           #panel.grid.major = element_blank(),
           panel.grid.major = element_line(colour='white', size=0.5),
           panel.grid.minor = element_line(colour='white', size=0.2),
           #panel.grid.minor = element_blank(),
           panel.border = element_rect(colour = "black", fill=NA, size=rel(1.1)),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.direction = "horizontal",
           legend.key.size= unit(0.2, "cm"),
           #egend.margin = unit(0, "cm"),
           legend.title = element_blank(),
           strip.background = element_blank(),
           strip.text = element_blank()
           #strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           #strip.text = element_text(face="bold")
           
   ))
  
}

##Function to get number format for x axis
#--------------------------------------------------------------------------
plain <- function(x,...) {
  format(x, ..., scientific = FALSE, drop0trailing = TRUE)
}
#--------------------------------------------------------------------------
##---------Additional options----------------------------------------------
# Format y axis tick mark labels to show exponents
p + scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(2^.x)))

## Flip a plot on it’s side using (coord_flip())
##-------------------------------------------------------------------------
## Remove axis.label (theme(axis.title.y = element_blank()))
###########################################################################
###########################################################################

##MAKE PLOTS---------------------------------------------------------------

#Rename variables for the plot
vars <- c("(Intercept)" = expression("Intercept"),
          "ECO_feedStrata02.high"=expression(paste("Feed.Strata"[HIGH],sep="")),
          "ECO_socialRoost02.social"=expression(paste(Roosting[SOCIAL],sep="")),
          "ECO_sdm02.yes"=expression(paste(SDM[PRESENT],sep="")),
          "ECO_minElev.div100"=expression("Min.Elev"),
          "IND_Cen.BodyCI2.div10"=expression(paste("Body.Cond",sep="")))
###########################################################################
#Getting input data
#--------------------------------------------------------------------------
P.exp.estimates <- read.csv("Plasm_expEstimates.csv", sep=",")

H.Ind.exp.estimates <- read.csv("Haem_expEstimates.csv", sep=",")
H.Ind.exp.estimates2 <- read.csv("SppInd_exp2.csv", sep=",")
H.cont.exp <- read.csv("H_continuous_vars.csv", sep=",")
H.cat.exp <- read.csv("H_categorical_vars.csv", sep=",")

##-------------------------------------------------------------------------
#Order factor levels
##-------------------------------------------------------------------------
##https://janhove.github.io/analysis/2016/08/18/ordering-factor-levels
source("https://raw.githubusercontent.com/janhove/janhove.github.io/master/RCode/sortLvls.R")
##-------------------------------------------------------------------------
P.exp.estimates$vn <- sortLvls.fnc(P.exp.estimates$vn, c(1,4,3,2))
levels(P.exp.estimates$vn)

H.Ind.exp.estimates$vn <- sortLvls.fnc(H.Ind.exp.estimates$vn, c(1,3, 2,4))
levels(H.Ind.exp.estimates$vn)

H.Ind.exp.estimates2$vn
H.Ind.exp.estimates2$vn <- sortLvls.fnc(H.Ind.exp.estimates2$vn, c(1,3, 2,4))
levels(H.Ind.exp.estimates2$vn)
##-------------------------------------------------------------------------
###########################################################################

#Vertical plots next to each other.
#Used for final MS - 7-6-2020

##########################################################################
##Plasmodium Plot
###########################################################################

Plasm1 <- ggplot(P.exp.estimates, aes(vn, Posterior.Mode)) + theme_few()+
  scale_y_log10(labels = plain)+
  geom_point(color= "#0072B2", fill="#0072B2", size=3)+
  geom_errorbar(aes(ymin = lower, ymax = upper, width=0.15), color="#0072B2")+ scale_x_discrete(labels =vars)+
  geom_hline(yintercept=1,lty='dashed', color= "grey50")+
  labs(title="Plasmodium",  y= "Infection Risk (Odds Ratio)", x= "")
#plot.margin=unit(c(3,3,1,3),"mm")
Plasm2 <- Plasm1 + (theme_Publication())

##########################################################################
##Haemoproteus Plot
###########################################################################
##------------------------------------------------------------------------------
Haem3 <- ggplot(H.Ind.exp.estimates, aes(vn, Posterior.Mode))+
  geom_point(color= "#D55E00", fill="#D55E00", size=3, shape=c(19,19,15,15))+ 
  geom_errorbar(aes(ymin = lower, ymax = upper, width=0.15), color="#D55E00")+ scale_x_discrete(labels =vars)+
  geom_hline(yintercept=1,lty="dashed", color= "grey50")+
  labs(title="Haemoproteus",  y= "", x= "")

#plot.margin=unit(c(1,3,3,3),"mm")
Haem4 <- Haem3+ scale_y_log10(labels = plain)+
  facet_wrap(grp ~ ., scales='free') + (theme_Publication())

###########################################################################
#Save the final plot (individual P and H plot)
##This works
ggsave(filename = "Plasm_Exp.png", plot = Plasm2, width = 6, height = 4, dpi = 400, units = "in")
ggsave(filename = "Haemm_Exp.png", plot = Haem2, width = 6, height = 4, dpi = 400, units = "in")

##-------------------------------------------------------------------------
# Put two plots one below other
##-------------------------------------------------------------------------

multiplot <- grid.arrange(Plasm2, Haem2, ncol=1)
ggsave(filename = "Plasm_Haem_Plot_v1_5-7.png", plot = multiplot, width = 5, height = 7, dpi = 400, units = "in")

##-------------------------------------------------------------------------
# Put two plots side by side (Final format used)- 4-27-20
##-------------------------------------------------------------------------

#Used for side by side plots on 
multiplot3 <- ggarrange(Plasm2, Haem4, labels = c("a", "b"), ncol = 2, align = "hv")
multiplot3 <- annotate_figure(multiplot3, bottom = text_grob("Ecological predictors", color = "black",size = 16, vjust=-0.5))

#Make 10x7..extra to fit vars labels properly
ggsave(filename = "Plasm_Haem_Plot_test4.png", plot = multiplot3, width = 10, height = 5, dpi = 400, units = "in")


####################################################################
#Plotting Phylovariance estimates
###################################################################
theme_Publication2 <- function(base_size=16, base_family="helvetica") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold.italic", size = rel(1.2)),
           text = element_text(),
           #panel.background = element_rect(colour = NA),
           panel.background =element_rect(fill = "grey92", colour = NA),
           plot.background = element_rect(colour = NA),
           #panel.border = element_rect(colour = NA),
           axis.title = element_text(size = 14),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           #axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           #panel.grid.major = element_blank(),
           #panel.grid.minor = element_blank(),
           panel.grid.major = element_line(colour='white', size=0.5),
           panel.grid.minor = element_line(colour='white', size=0.2),
           #panel.border = element_rect(colour = NA),
           panel.border = element_rect(colour = "black", fill=NA, size=rel(0.8)),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.direction = "horizontal",
           legend.key.size= unit(0.5, "cm"),
           legend.margin = unit(0, "cm"),
           legend.title = element_blank(),
           legend.text = element_text(face = "italic"),
           plot.margin=unit(c(4,4,4,4),"mm") #Plas
           #plot.margin=unit(c(3,2.5,3,1.5),"mm")  #Haem
           #strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           #strip.text = element_text(face="bold")
   ))
  
}

#getting data
post.var <- read.csv("PH_PhyloVar_estimates_plot.csv", sep=",")

post.var$ParSpp <- sortLvls.fnc(post.var$ParSpp, c(2,1))
levels(post.var$ParSpp)

p <- ggplot(post.var, aes(x=variance, y=V1*100, fill=ParSpp)) + 
  geom_bar(stat="identity", color="grey10", position="dodge")+
  geom_errorbar(aes(ymin=lower*100, ymax=upper*100), width=0.2,
                position=position_dodge(0.9)) 
 
# Finished bar plot
p2 <- p+labs(x="", y = "% Phylogenetic Variance")+
   (theme_Publication2()) +
  scale_fill_manual(values=alpha(c('#0072B2', "#D55E00"), 0.8))

#ppi <- 500
#png("Phylovar_plot_3-4.png", width = 3.3, height = 4.4, res=ppi, units = "in")

ggsave(filename = "Phylovar_plot_10-4-20.tiff", plot = p2, width = 4, height = 5, dpi = 800, units = "in")


#################################################################################
#END CODE
