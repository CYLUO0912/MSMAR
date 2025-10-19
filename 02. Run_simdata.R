##############################################################################-#
## The following codes replicate our results in simulation study               #
## Author: Caiying Luo                                                         #
##############################################################################-#


################################################################################
## Model the simulation datasets for each scenario                             #
## This step will cost several days due to the large number of simulation times#
## You may need a super computer to run the codes                              #
################################################################################
### Load data and function -----
path <- getwd()
`%+%` <- function(x,y) paste0(x,y)
source(path%+%"/fun_MCMAR.R")
source(path%+%"/fun_MSEMAR.R")
source(path%+%"/fun_MSLMAR.R")
load(path%+%"/data/simdata.Rdata")
library(mvmeta)

### Fit MMR, MCMAR, and MSMAR with or without covariates -----
method <- "ml" 
smvcontrol <- list(maxiter = 200,factr = 1e7,hessian = F,opt.iter = 1,opt.iter.show = F)

scennames <- names(simdata)[1:length(simdata)]
for (scen in scennames ) {
  data <- simdata[[scen]]
  Sall <- data$trueparameter$Sall
  covariate <- data$trueparameter$covariate
  Wmatrix <- data$trueparameter$Wmatrix
  Cmatrix <- data$trueparameter$Cmatrix
  
  nsim <- 1000
  res_MMR_intercept <- res_MMR_covariate <- vector(mode = "list",length = nsim) 
  res_MCMAR_intercept <- res_MCMAR_covariate <- vector(mode = "list",length = nsim)
  res_MSLMAR_intercept <- res_MSLMAR_covariate <- vector(mode = "list",length = nsim)
  res_MSEMAR_intercept <- res_MSEMAR_covariate <- vector(mode = "list",length = nsim)
  
  for (i in 1:nsim) {
    yall <- data$Simy[[i]]
    #MMR、MCMAR、MSLMAR、MSEMAR
    fit0_intercept <- mvmeta(yall,Sall,method = method)
    fit0_covariate <- mvmeta(yall ~ covariate[,-1],Sall,method = method)
    
    fit1_intercept <- smvmeta(yall,S=Sall,Cmatrix = Cmatrix,
                              method=method,control = smvcontrol)
    fit1_covariate <- smvmeta(yall ~ covariate[,-1],S=Sall,Cmatrix = Cmatrix,
                              method=method,control = smvcontrol)
    
    fit2_intercept <- Lmvmeta(yall,S=Sall,Wmatrix = Wmatrix,
                                   method=method,control = smvcontrol)
    fit2_covariate <- Lmvmeta(yall ~ covariate[,-1],S=Sall,Wmatrix = Wmatrix,
                                   method=method,control = smvcontrol)
    
    fit3_intercept <- Emvmeta(yall,S=Sall,Wmatrix = Wmatrix,
                              method=method,control = smvcontrol)
    fit3_covariate <- Emvmeta(yall ~ covariate[,-1],S=Sall,Wmatrix = Wmatrix,
                              method=method,control = smvcontrol)
    
    res_MMR_intercept[[i]] <- fit0_intercept
    res_MMR_covariate[[i]] <- fit0_covariate
    res_MCMAR_intercept[[i]] <- fit1_intercept
    res_MCMAR_covariate[[i]] <- fit1_covariate
    res_MSLMAR_intercept[[i]] <- fit2_intercept[-c(15:17,32)]
    res_MSLMAR_covariate[[i]] <- fit2_covariate[-c(15:17,32)]
    res_MSEMAR_intercept[[i]] <- fit3_intercept[-c(15:17,32)]
    res_MSEMAR_covariate[[i]] <- fit3_covariate[-c(15:17,32)]
    
    if(i%%100 == 0) save(res_MMR_intercept,res_MMR_covariate,
                         res_MCMAR_intercept,res_MCMAR_covariate,
                         res_MSLMAR_intercept,res_MSLMAR_covariate,
                         res_MSEMAR_intercept,res_MSEMAR_covariate,
                         file = "simres_"%+%scen%+%"_"%+%method%+%".Rdata")
    cat(i," ")
  }
  
}

