##############################################################################-#
## The following codes replicate our results of the computation time in  #     #
## simulation study                                                            #
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
source(path%+%"/fun_MSEMAR.R")
source(path%+%"/fun_MSLMAR.R")
source(path%+%"/fun_MCMAR.R")
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
  res_time <- as.data.frame(matrix(NA,nsim,8))
  colnames(res_time) <- c("MMR_intercept","MMR_covariate","MCMAR_intercept","MCMAR_covariate",
                          "MSLMAR_intercept","MSLMAR_covariate","MSEMAR_intercept","MSEMAR_covariate")
  
  for (i in 1:nsim) {
    yall <- data$Simy[[i]]
    #MMR、MCMAR、MSLMAR、MSEMAR
    fit0_intercept <- system.time(mvmeta(yall,Sall,method = method))["elapsed"]
    fit0_covariate <- system.time(mvmeta(yall ~ covariate[,-1],Sall,method = method))["elapsed"]
    
    fit1_intercept <- system.time(smvmeta(yall,S=Sall,Cmatrix = Cmatrix,
                                          method=method,control = smvcontrol))["elapsed"]
    fit1_covariate <- system.time(smvmeta(yall ~ covariate[,-1],S=Sall,Cmatrix = Cmatrix,
                                          method=method,control = smvcontrol))["elapsed"]
    
    fit2_intercept <- system.time(Lmvmeta(yall,S=Sall,Wmatrix = Wmatrix,
                                          method=method,control = smvcontrol))["elapsed"]
    fit2_covariate <- system.time(Lmvmeta(yall ~ covariate[,-1],S=Sall,Wmatrix = Wmatrix,
                                          method=method,control = smvcontrol))["elapsed"]
    
    fit3_intercept <- system.time(Emvmeta(yall,S=Sall,Wmatrix = Wmatrix,
                                          method=method,control = smvcontrol))["elapsed"]
    fit3_covariate <- system.time(Emvmeta(yall ~ covariate[,-1],S=Sall,Wmatrix = Wmatrix,
                                          method=method,control = smvcontrol))["elapsed"]
    
    res_time[i,1] <- fit0_intercept
    res_time[i,2] <- fit0_covariate
    res_time[i,3] <- fit1_intercept
    res_time[i,4] <- fit1_covariate
    res_time[i,5] <- fit2_intercept
    res_time[i,6] <- fit2_covariate
    res_time[i,7] <- fit3_intercept
    res_time[i,8] <- fit3_covariate
    
    if(i%%100 == 0) save(res_time,
                         file = "simrestime_"%+%scen%+%"_"%+%method%+%".Rdata")
    cat(i," ")
  }
  
}

