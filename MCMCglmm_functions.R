#Function for MCMCglmm
#data = whatever data needs to be used
#vars = whatever predictors to be used for the model
#--------------------------------------------------------------------
#Plasmodium
P.mcmcrun <- function(data,vars) {
  Call <- as.formula(paste("Inf_Plas~", paste(vars, collapse="+")))
  do.call("MCMCglmm", args = list(Call,data=data,random=~Species+Minor.Loc,  
                                  family = "categorical",ginverse=list(Species=Ainv.m1),
                                  slice=T,prior=prior1,nitt=NITT,burnin=BITT,thin=THIN))
}

##-----------------------------------------------------------------------
#Haemoproteus
H.mcmcrun <- function(data,vars) {
  Call <- as.formula(paste("Inf_Haem~", paste(vars, collapse="+")))
  do.call("MCMCglmm", args = list(Call,data=data,random=~Species+Minor.Loc,  
                                  family = "categorical",ginverse=list(Species=Ainv.m1),
                                  slice=T,prior=prior1,nitt=NITT,burnin=BITT,thin=THIN))
}

###################################################################
#Getting MCMC model full outputs for P and H together
##################################################################
##Source code: https://gkhajduk.github.io/2017-10-25-cleanMCMCglmm/

clean.MCMC <- function(x) {
  sols <- summary(x)$solutions  ## pull out relevant info from model summary
  Gcovs <- summary(x)$Gcovariances
  #Rcovs <- summary(x)$Rcovariances
  fixed <- data.frame(row.names(sols), sols, row.names = NULL)  ## convert to dataframes with the row.names as the first col
  random <- data.frame(row.names(Gcovs), Gcovs, row.names = NULL)
  #residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
  names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  ## change the columns names to variable, so they all match
  names(random)[names(random) == "row.names.Gcovs."] <- "variable"
  #names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
  fixed$effect <- "fixed"  ## add ID column for type of effect (fixed, random, residual)
  random$effect <- "random"
  #residual$effect <- "residual"
  modelTerms <- as.data.frame(bind_rows(fixed, random))  # merge it all together
}

########################################################################3
#We write a small utility function to automate the process of combining the fixed-effect (Sol) and variance-covariance (VCV) parameter chains and making the trace plot:
#Only works with a single model
tfun <- function(mod) {
        plotTrace(as.mcmc(cbind(mod$Sol,mod$VCV)))
}

#########################################################################
Pmod.all.res <- lapply(resMCMC.List[["P"]],function(x)summary(x)$solutions)
############################################################################
plotTrace(as.mcmc(cbind(H.Ind.mod2$Sol, H.Ind.mod2$VCV)), axes=TRUE,las=1)
