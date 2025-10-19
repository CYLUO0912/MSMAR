##############################################################################-#
## The following codes replicate our results in simulation study               #
## Author: Caiying Luo                                                         #
##############################################################################-#


################################################################################
##### Get performance for simulation study######################################
################################################################################

########## Load function and data -----
library(mvmeta)
`%+%` <- function(x,y) paste0(x,y)
`%>%` <- magrittr::`%>%`
load("data\\simdata.Rdata")


GetErrorBeta0 <- function(model,truevalues){
  sqrt(sum((model$coefficients[1,]/truevalues - 1)^2))
}

GetErrorBeta1 <- function(model,truevalues){
  
  if (nrow(model$coefficients)==1) {
    coef <- model$coefficients
    mat <- matrix(0,3,4)
    model$coefficients <- rbind(coef,mat)
  }
  
  if (any(truevalues == 0))
    sqrt(sum((model$coefficients[2:4,]-truevalues)^2))
  else 
    sqrt(sum((model$coefficients[2:4,]/truevalues - 1)^2))
}

GetAIC <- function(model){
  class(model) <- "mvmeta"
  df <- attr(logLik(model),"df")
  -logLik(model) * 2 + 2 * df
}

waldtest <- function(model,var){
  class(model) <- "mvmeta"
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  pvalue <- 1-pchisq(waldstat,df)
  wald <- c(waldstat,df,pvalue)
  wald
}

Getcoveragebeta <- function(model,mu,var){
  class(model) <- "mvmeta"
  ind <- grep(var,names(coef(model)))
  s2 <- vcov(model)[ind,ind]
  est <- coef(model)[ind] - mu
  b <- est %*% solve(s2) %*% est
  pchisq(q = b,df = length(est)) < 0.95
}


method <- "ml"
scennames <- names(simdata)
########## get performance for specific scenario ------
ARMAE_beta <- aic <- power <- coverage95 <- NULL

for (scen in scennames) {
    cat(scen,"")
    load(getwd()%+%"/data/Simdataresults/simres_"%+%scen%+%"_"%+%method%+%".Rdata")
    truedata <- simdata[[scen]]
    m <- 108
    k <- 4
    if (grepl("MCMAR-int",scen)) {truedata$trueparameter$beta <- rbind(truedata$trueparameter$beta,matrix(0,3,4))}
    
    ## beta -----------
    beta0_true <- truedata$trueparameter$beta[1,]
    beta1_true <- truedata$trueparameter$beta[2:4,]
    
    beta0_RMAE_MMR_intercept <- sapply(res_MMR_intercept, GetErrorBeta0,
                                       truevalues = beta0_true) %>% mean %>% round(3) 
    beta0_RMAE_MMR_covariate <- sapply(res_MMR_covariate, GetErrorBeta0,
                                       truevalues = beta0_true) %>% mean %>% round(3) 
    beta0_RMAE_MCMAR_intercept <- sapply(res_MCMAR_intercept, GetErrorBeta0,
                                         truevalues = beta0_true) %>% mean %>% round(3) 
    beta0_RMAE_MCMAR_covariate <- sapply(res_MCMAR_covariate, GetErrorBeta0,
                                         truevalues = beta0_true) %>% mean %>% round(3) 
    beta0_RMAE_MSLMAR_intercept <- sapply(res_MSLMAR_intercept, GetErrorBeta0,
                                          truevalues = beta0_true) %>% mean %>% round(3) 
    beta0_RMAE_MSLMAR_covariate <- sapply(res_MSLMAR_covariate, GetErrorBeta0,
                                          truevalues = beta0_true) %>% mean %>% round(3) 
    beta0_RMAE_MSEMAR_intercept <- sapply(res_MSEMAR_intercept, GetErrorBeta0,
                                          truevalues = beta0_true) %>% mean %>% round(3) 
    beta0_RMAE_MSEMAR_covariate <- sapply(res_MSEMAR_covariate, GetErrorBeta0,
                                          truevalues = beta0_true) %>% mean %>% round(3) 
    
    
    beta1_RMAE_MMR_intercept <- sapply(res_MMR_intercept, GetErrorBeta1,
                                       truevalues = beta1_true) %>% mean %>% round(3) 
    beta1_RMAE_MMR_covariate <- sapply(res_MMR_covariate, GetErrorBeta1,
                                       truevalues = beta1_true) %>% mean %>% round(3) 
    beta1_RMAE_MCMAR_intercept <- sapply(res_MCMAR_intercept, GetErrorBeta1,
                                         truevalues = beta1_true) %>% mean %>% round(3) 
    beta1_RMAE_MCMAR_covariate <- sapply(res_MCMAR_covariate, GetErrorBeta1,
                                         truevalues = beta1_true) %>% mean %>% round(3) 
    beta1_RMAE_MSLMAR_intercept <- sapply(res_MSLMAR_intercept, GetErrorBeta1,
                                          truevalues = beta1_true) %>% mean %>% round(3) 
    beta1_RMAE_MSLMAR_covariate <- sapply(res_MSLMAR_covariate, GetErrorBeta1,
                                          truevalues = beta1_true) %>% mean %>% round(3) 
    beta1_RMAE_MSEMAR_intercept <- sapply(res_MSEMAR_intercept, GetErrorBeta1,
                                          truevalues = beta1_true) %>% mean %>% round(3) 
    beta1_RMAE_MSEMAR_covariate <- sapply(res_MSEMAR_covariate, GetErrorBeta1,
                                          truevalues = beta1_true) %>% mean %>% round(3) 
    
    
    beta_RMAE <- c(beta0_RMAE_MMR_intercept,beta0_RMAE_MCMAR_intercept,beta0_RMAE_MSLMAR_intercept,beta0_RMAE_MSEMAR_intercept,
                   beta0_RMAE_MMR_covariate,beta0_RMAE_MCMAR_covariate,beta0_RMAE_MSLMAR_covariate,beta0_RMAE_MSEMAR_covariate,
                   beta1_RMAE_MMR_intercept,beta1_RMAE_MCMAR_intercept,beta1_RMAE_MSLMAR_intercept,beta1_RMAE_MSEMAR_intercept,
                   beta1_RMAE_MMR_covariate,beta1_RMAE_MCMAR_covariate,beta1_RMAE_MSLMAR_covariate,beta1_RMAE_MSEMAR_covariate)
    ARMAE_beta <- rbind(ARMAE_beta,beta_RMAE)
    colnames(ARMAE_beta) <- c("beta0_MMR_intercept","beta0_MCMAR_intercept","beta0_MSLMAR_intercept","beta0_MSEMAR_intercept",
                              "beta0_MMR_covariate","beta0_MCMAR_covariate","beta0_MSLMAR_covariate","beta0_MSEMAR_covariate",
                              "beta1_MMR_intercept","beta1_MCMAR_intercept","beta1_MSLMAR_intercept","beta1_MSEMAR_intercept",
                              "beta1_MMR_covariate","beta1_MCMAR_covariate","beta1_MSLMAR_covariate","beta1_MSEMAR_covariate")
    
    ## AIC----------
    aic_MMR_intercept <- sapply(res_MMR_intercept, GetAIC) %>% mean %>% round(3) 
    aic_MMR_covariate <- sapply(res_MMR_covariate, GetAIC) %>% mean %>% round(3) 
    aic_MCMAR_intercept <- sapply(res_MCMAR_intercept, GetAIC) %>% mean %>% round(3) 
    aic_MCMAR_covariate <- sapply(res_MCMAR_covariate, GetAIC) %>% mean %>% round(3) 
    aic_MSLMAR_intercept <- sapply(res_MSLMAR_intercept, GetAIC) %>% mean %>% round(3) 
    aic_MSLMAR_covariate <- sapply(res_MSLMAR_covariate, GetAIC) %>% mean %>% round(3) 
    aic_MSEMAR_intercept <- sapply(res_MSEMAR_intercept, GetAIC) %>% mean %>% round(3) 
    aic_MSEMAR_covariate <- sapply(res_MSEMAR_covariate, GetAIC) %>% mean %>% round(3) 
    
    aici <- c(aic_MMR_intercept,aic_MCMAR_intercept,aic_MSLMAR_intercept,aic_MSEMAR_intercept,
              aic_MMR_covariate,aic_MCMAR_covariate,aic_MSLMAR_covariate,aic_MSEMAR_covariate)
    aic <- rbind(aic,aici)
    colnames(aic) <- c("aic_MMR_intercept","aic_MCMAR_intercept","aic_MSLMAR_intercept","aic_MSEMAR_intercept",
                       "aic_MMR_covariate","aic_MCMAR_covariate","aic_MSLMAR_covariate","aic_MSEMAR_covariate")
    
    ## identification of the covariate -----
    p_MMR_pop100 <- sapply(res_MMR_covariate, function(x) waldtest(x,"pop100")[3])
    p_MMR_Phigh <- sapply(res_MMR_covariate, function(x) waldtest(x,"Phigh")[3])
    p_MMR_Punem <- sapply(res_MMR_covariate, function(x) waldtest(x,"Punem")[3])
    
    p_MCMAR_pop100 <- sapply(res_MCMAR_covariate, function(x) waldtest(x,"pop100")[3])
    p_MCMAR_Phigh <- sapply(res_MCMAR_covariate, function(x) waldtest(x,"Phigh")[3])
    p_MCMAR_Punem <- sapply(res_MCMAR_covariate, function(x) waldtest(x,"Punem")[3])
    
    p_MSLMAR_pop100 <- sapply(res_MSLMAR_covariate, function(x) waldtest(x,"pop100")[3])
    p_MSLMAR_Phigh <- sapply(res_MSLMAR_covariate, function(x) waldtest(x,"Phigh")[3])
    p_MSLMAR_Punem <- sapply(res_MSLMAR_covariate, function(x) waldtest(x,"Punem")[3])
    
    p_MSEMAR_pop100 <- sapply(res_MSEMAR_covariate, function(x) waldtest(x,"pop100")[3])
    p_MSEMAR_Phigh <- sapply(res_MSEMAR_covariate, function(x) waldtest(x,"Phigh")[3])
    p_MSEMAR_Punem <- sapply(res_MSEMAR_covariate, function(x) waldtest(x,"Punem")[3])
    
    
    poweri <- c(sum(p_MMR_pop100<0.05),sum(p_MMR_Phigh<0.05),sum(p_MMR_Punem<0.05),
                sum(p_MCMAR_pop100<0.05),sum(p_MCMAR_Phigh<0.05),sum(p_MCMAR_Punem<0.05),
                sum(p_MSLMAR_pop100<0.05),sum(p_MSLMAR_Phigh<0.05),sum(p_MSLMAR_Punem<0.05),
                sum(p_MSEMAR_pop100<0.05),sum(p_MSEMAR_Phigh<0.05),sum(p_MSEMAR_Punem<0.05))
    
    power <- rbind(power,poweri) %>% round(3)
    colnames(power) <- c("MMR_pop100","MMR_Phigh","MMR_Punem",
                         "MCMAR_pop100","MCMAR_Phigh","MCMAR_Punem",
                         "MSLMAR_pop100","MSLMAR_Phigh","MSLMAR_Punem",
                         "MSEMAR_pop100","MSEMAR_Phigh","MSEMAR_Punem")
    
    ## 95% coverage rates for beta -----
    beta0_true <- truedata$trueparameter$beta[1,]
    beta1pop100_true <- truedata$trueparameter$beta[2,]
    beta1Phigh_true <- truedata$trueparameter$beta[3,]
    beta1Punem_true <- truedata$trueparameter$beta[4,]
    
    va <- "Intercept"
    beta0_coverate_MMR_intercept <- sapply(res_MMR_intercept, Getcoveragebeta,
                                           mu = beta0_true,var=va) %>% mean %>% round(3) 
    beta0_coverate_MMR_covariate <- sapply(res_MMR_covariate, Getcoveragebeta,
                                           mu = beta0_true,var=va) %>% mean %>% round(3) 
    beta0_coverate_MCMAR_intercept <- sapply(res_MCMAR_intercept, Getcoveragebeta,
                                             mu = beta0_true,var=va) %>% mean %>% round(3) 
    beta0_coverate_MCMAR_covariate <- sapply(res_MCMAR_covariate, Getcoveragebeta,
                                             mu = beta0_true,var=va) %>% mean %>% round(3) 
    beta0_coverate_MSLMAR_intercept <- sapply(res_MSLMAR_intercept, Getcoveragebeta,
                                              mu = beta0_true,var=va) %>% mean %>% round(3) 
    beta0_coverate_MSLMAR_covariate <- sapply(res_MSLMAR_covariate, Getcoveragebeta,
                                              mu = beta0_true,var=va) %>% mean %>% round(3) 
    beta0_coverate_MSEMAR_intercept <- sapply(res_MSEMAR_intercept, Getcoveragebeta,
                                              mu = beta0_true,var=va) %>% mean %>% round(3) 
    beta0_coverate_MSEMAR_covariate <- sapply(res_MSEMAR_covariate, Getcoveragebeta,
                                              mu = beta0_true,var=va) %>% mean %>% round(3) 
    #beta1
    va <- "pop100"
    beta1pop100_coverate_MMR_covariate <- sapply(res_MMR_covariate, Getcoveragebeta,
                                                 mu = beta1pop100_true,var=va) %>% mean %>% round(3) 
    beta1pop100_coverate_MCMAR_covariate <- sapply(res_MCMAR_covariate, Getcoveragebeta,
                                                   mu = beta1pop100_true,var=va) %>% mean %>% round(3) 
    beta1pop100_coverate_MSLMAR_covariate <- sapply(res_MSLMAR_covariate, Getcoveragebeta,
                                                    mu = beta1pop100_true,var=va) %>% mean %>% round(3) 
    beta1pop100_coverate_MSEMAR_covariate <- sapply(res_MSEMAR_covariate, Getcoveragebeta,
                                                    mu = beta1pop100_true,var=va) %>% mean %>% round(3) 
    
    va <- "Phigh"
    beta1Phigh_coverate_MMR_covariate <- sapply(res_MMR_covariate, Getcoveragebeta,
                                                mu = beta1Phigh_true,var=va) %>% mean %>% round(3) 
    beta1Phigh_coverate_MCMAR_covariate <- sapply(res_MCMAR_covariate, Getcoveragebeta,
                                                  mu = beta1Phigh_true,var=va) %>% mean %>% round(3) 
    beta1Phigh_coverate_MSLMAR_covariate <- sapply(res_MSLMAR_covariate, Getcoveragebeta,
                                                   mu = beta1Phigh_true,var=va) %>% mean %>% round(3) 
    beta1Phigh_coverate_MSEMAR_covariate <- sapply(res_MSEMAR_covariate, Getcoveragebeta,
                                                   mu = beta1Phigh_true,var=va) %>% mean %>% round(3) 
    
    va <- "Punem"
    beta1Punem_coverate_MMR_covariate <- sapply(res_MMR_covariate, Getcoveragebeta,
                                                mu = beta1Punem_true,var=va) %>% mean %>% round(3) 
    beta1Punem_coverate_MCMAR_covariate <- sapply(res_MCMAR_covariate, Getcoveragebeta,
                                                  mu = beta1Punem_true,var=va) %>% mean %>% round(3) 
    beta1Punem_coverate_MSLMAR_covariate <- sapply(res_MSLMAR_covariate, Getcoveragebeta,
                                                   mu = beta1Punem_true,var=va) %>% mean %>% round(3) 
    beta1Punem_coverate_MSEMAR_covariate <- sapply(res_MSEMAR_covariate, Getcoveragebeta,
                                                   mu = beta1Punem_true,var=va) %>% mean %>% round(3) 
    
    a <- c(beta0_coverate_MMR_intercept,beta0_coverate_MCMAR_intercept,beta0_coverate_MSLMAR_intercept,beta0_coverate_MSEMAR_intercept,
           beta0_coverate_MMR_covariate,beta0_coverate_MCMAR_covariate,beta0_coverate_MSLMAR_covariate,beta0_coverate_MSEMAR_covariate,
           beta1pop100_coverate_MMR_covariate,beta1pop100_coverate_MCMAR_covariate,beta1pop100_coverate_MSLMAR_covariate,beta1pop100_coverate_MSEMAR_covariate,
           beta1Phigh_coverate_MMR_covariate,beta1Phigh_coverate_MCMAR_covariate,beta1Phigh_coverate_MSLMAR_covariate,beta1Phigh_coverate_MSEMAR_covariate,
           beta1Punem_coverate_MMR_covariate,beta1Punem_coverate_MCMAR_covariate,beta1Punem_coverate_MSLMAR_covariate,beta1Punem_coverate_MSEMAR_covariate)
    coverage95 <- rbind(coverage95,a)
    colnames(coverage95) <- c("beta0_MMR_intercept","beta0_MCMAR_intercept","beta0_MSLMAR_intercept","beta0_MSEMAR_intercept",
                              "beta0_MMR_covariate","beta0_MCMAR_covariate","beta0_MSLMAR_covariate","beta0_MSEMAR_covariate", 
                              "beta1pop100_MMR_covariate","beta1pop100_MCMAR_covariate","beta1pop100_MSLMAR_covariate","beta1pop100_MSEMAR_covariate",
                              "beta1Phigh_MMR_covariate","beta1Phigh_MCMAR_covariate","beta1Phigh_MSLMAR_covariate","beta1Phigh_MSEMAR_covariate",
                              "beta1Punem_MMR_covariate","beta1Punem_MCMAR_covariate","beta1Punem_MSLMAR_covariate","beta1Punem_MSEMAR_covariate")

}

rm(list=ls()[which(!ls() %in% c('ARMAE_beta','ARMAE_scc','aic','power','coverage95','scennames','%+%','method','%>%'))])
rownames(ARMAE_beta) <- rownames(aic) <- rownames(power) <- rownames(coverage95) <- scennames
powerbind <- sapply(1:4, function(x) rowMeans(power[,(3*x-2):(3*x)])) %>% round(3)
coverage95bind <- sapply(1:4, function(x) rowMeans(coverage95[,c(x+8,x+12,x+16)])) %>% round(3)
colnames(powerbind) <- colnames(coverage95bind) <- c("MMR","MCMAR","MSLMAR","MSEMAR")


xlsx::write.xlsx(ARMAE_beta,file = "performance-results-"%+%method%+%".xlsx",sheetName = "ARMAE_beta")
xlsx::write.xlsx(aic,file = "performance-results-"%+%method%+%".xlsx",sheetName = "aic",append = T)
xlsx::write.xlsx(power,file = "performance-results-"%+%method%+%".xlsx",sheetName = "power",append = T)
xlsx::write.xlsx(powerbind,file = "performance-results-"%+%method%+%".xlsx",sheetName = "powerbind",append = T)
xlsx::write.xlsx(coverage95,file = "performance-results-"%+%method%+%".xlsx",sheetName = "coverage95",append = T)
xlsx::write.xlsx(coverage95bind,file = "performance-results-"%+%method%+%".xlsx",sheetName = "coverage95bind",append = T)

