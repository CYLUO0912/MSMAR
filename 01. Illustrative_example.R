##############################################################################-#
## The following codes replicate our results in illustrative examples               #
## Author: Caiying Luo                                                         #
##############################################################################-#


################################################################################
# R code for the analysis in:
#
##########Main case example of temperature-mortality relationships##############
################################################################################
`%+%` <- function(x,y) paste0(x,y)
path <- getwd()
source(path%+%"/fun_MCMAR.R")
source(path%+%"/fun_MSEMAR.R")
source(path%+%"/fun_MSLMAR.R")

# Wald test for MMR, MCMAR and MSMAR
waldtest <- function(model){
  m <- model
  coef <- m$coefficients[2,]
  vcov <- m$vcov[seq(2,8,2),seq(2,8,2)]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  pvalue <- 1-pchisq(waldstat,df)
  wald <- c(waldstat,df,pvalue)
  wald
}

# Load packages
library(mvmeta) ; library(dlnm) ; library(scales);
library(sp) ; library(spdep); library(maps);library(mixmeta)

# Load coef/vcov from first-stage models
tmeanpar <- read.csv(file="data/tmeanpar.csv")
yall <- as.matrix(tmeanpar[,grep("coef", names(tmeanpar))])
Sall <- as.matrix(tmeanpar[,grep("vcov", names(tmeanpar))])
# Link with census data
cityind <- tmeanpar[,1:4,]
latlong <- read.csv("data/latlong.csv")
cityind <- merge(cityind, latlong, by="city")
citycensus <- read.csv("data/citycensus.csv")
cityind <- merge(cityind, citycensus, by="city")

################################################################################
##### Get the spatially adjacent matrix            #############################
################################################################################
## k-nearest neighbor method-----
k <- 3
temp <- cityind
coordinates(temp) <- ~ long + lat 
nb <- knearneigh(temp, k=k) 
nb <- knn2nb(nb)
plot(nb,coordinates(temp),pch = 16)
W <- nb2mat(nb)
Wmatrix <- W
n <- ncol(W)
C <- W * apply(W, 1, function(x) sum(x!=0))
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    if(C[i,j]!= C[j,i]){
      C[i,j] <- 1
      C[j,i] <- 1
    } 
  }
}
R <- diag(rowSums(C)) - C
Cmatrix <- diag(n) - R

################################################################################
############# Compare MMR with MCMAR and MSMAR with only intercept##############
################################################################################
dfvar <- 4

##MMR
fit0 <- mvmeta(yall,Sall,method = "ml")
aic0 <- (fit0$logLik * (-2) + 2 * (dfvar) + dfvar*(dfvar+1)) %>% round(1)

##EMMR
fit1 <- mixmeta(as.matrix(tmeanpar[,5:8]), tmeanpar[,9:18], data=tmeanpar, random=~1|state/city,method="ml")
aic1 <- summary(fit1)$AIC  %>% round(1)

##MCMAR
system.time(
  fit2 <- smvmeta(yall,S=Sall,Cmatrix = Cmatrix,method= "ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic2 <- (fit2$logLik * (-2) + 2 * (dfvar) + dfvar*(dfvar+1) + 2)  %>% round(1)

##MSLMAR
system.time(
  fit3 <- Lmvmeta(yall,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic3 <- (fit3$logLik * (-2) + 2 * (dfvar) + dfvar*(dfvar+1) + 2)  %>% round(1)
fit3coe <- Outoerr_L(fit3)

##MSEMAR
system.time(
  fit4 <- Emvmeta(yall,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic4 <- (fit4$logLik * (-2) + 2 * (dfvar) + dfvar*(dfvar+1) + 2)  %>% round(1)
fit4coe <- Outoerr_E(fit4)

aic0; aic1; aic2; aic3; aic4
fit2$rho; fit3$Rho; fit4$Lambda
1-pchisq((fit2$logLik - fit0$logLik) * 2,df = 1)
1-pchisq((fit3$logLik - fit0$logLik) * 2,df = 1)
1-pchisq((fit4$logLik - fit0$logLik) * 2,df = 1)


# Average temperature-mortality association
# Load Average temperature distribution across cities
avgtmeansum <- read.csv("data/avgtmeansum.csv")
tmean <- avgtmeansum$tmean
# Define spline transformation originally used in first-stage models
knots <- tmean[avgtmeansum$perc %in% paste0(c(50,90), ".0%")]
bvar <- onebasis(tmean, fun="bs", degree=2, knots=knots)
# Define the centering point (at point of minimum risk)
cen <- tmean[avgtmeansum$perc %in% "50.0%"]
# plotting labels
xperc <- c(0,1,5,25,50,75,95,99,100)
xval <- tmean[avgtmeansum$perc %in% paste0(xperc, ".0%")]
# Predict the association
cp <- crosspred(bvar, coef=coef(fit0), vcov=vcov(fit0), model.link="log",
                at=tmean, cen=cen)
cp1 <- crosspred(bvar, coef= as.numeric(fit1$coefficients), vcov=fit1$vcov, model.link="log",
                 at=tmean, cen=cen)
cp2 <- crosspred(bvar, coef= as.numeric(fit3coe$coefsim), vcov=fit3coe$vcovsim, model.link="log",
                 at=tmean, cen=cen)
cp3 <- crosspred(bvar, coef= as.numeric(fit4coe$coefsim), vcov=fit4coe$vcovsim, model.link="log",
                 at=tmean, cen=cen)

##MMR vs MCMAR
bluetrans <- rgb(0, 0,255, 120, maxColorValue=255) 
par(mfrow=c(1,3))
plot(cp, ylim=c(0.9,1.4), xlab="Temperature percentile", ylab="RR",
      xaxt="n", cex.axis=0.8, col="red",
     lwd=2)
axis(1, at=xval, labels=paste0(xperc, "%"))
lines(cp1,col="blue",lty=4,lwd=2,ci="area",ci.arg=list(density=20,col=bluetrans))
legend(x="top",inset =0, legend=c("MMR      (AIC = "%+%aic0%+%")", "MCMAR (AIC = "%+%aic2%+%")"),
       lwd=1.5, lty=1, col=c("red", "blue"), bty="n",ncol=1, cex=1.0)

##MMR vs MSLMAR
plot(cp, ylim=c(0.9,1.4), xlab="Temperature percentile", ylab="RR",
     xaxt="n",  cex.axis=0.8, col="red",
     lwd=2)
axis(1, at=xval, labels=paste0(xperc, "%"))
lines(cp2,col="blue",lty=4,lwd=2,ci="area",ci.arg=list(density=20,col=bluetrans))
legend(x="top",inset =0, legend=c("MMR      (AIC = "%+%aic0%+%")", "MSLMAR (AIC = "%+%aic3%+%")"),
       lwd=1.5, lty=1, col=c("red", "blue"), bty="n",ncol=1, cex=1.0)

##MMR vs MSEMAR
plot(cp, ylim=c(0.9,1.4), xlab="Temperature percentile", ylab="RR",
     xaxt="n", cex.axis=0.8, col="red",
     lwd=2)
axis(1, at=xval, labels=paste0(xperc, "%"))
lines(cp3,col="blue",lty=4,lwd=2,ci="area",ci.arg=list(density=20,col=bluetrans))
legend(x="top",inset =0, legend=c("MMR      (AIC = "%+%aic0%+%")", "MSEMAR (AIC = "%+%aic4%+%")"),
       lwd=1.5, lty=1, col=c("red", "blue"), bty="n",ncol=1, cex=1.0)


################################################################################
######### Compare models of MMR with MCMAR and MSMAR with covariates ###########
################################################################################
varnum <- 1+1
dfvarm <- dfvar*varnum
#~cityind$pop100
##MMR
fit5 <- mvmeta(yall~cityind$pop100,Sall,method = "ml")
aic5 <- (fit5$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) )  %>% round(1)

##EMMR
fit6 <- mixmeta(as.matrix(tmeanpar[,5:8])~cityind$pop100, tmeanpar[,9:18], data=tmeanpar, random=~1|state/city,method="ml")
aic6 <- summary(fit6)$AIC  %>% round(1)

#MCMAR
system.time(
  fit7 <- smvmeta(yall~cityind$pop100,S=Sall,Cmatrix = Cmatrix,method= "ml",
                control = list(maxiter = 200,factr = 1e7,hessian = T,
                               opt.iter = 1,opt.iter.show = T))
)
aic7 <- (fit7$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

##MSLMAR
system.time(
  fit8 <- Lmvmeta(yall~cityind$pop100,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic8 <- (fit8$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

##MSEMAR
system.time(
  fit9 <- Emvmeta(yall~cityind$pop100,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic9 <- (fit9$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

aic5; aic6; aic7; aic8; aic9
fit7$rho; fit8$Rho; fit9$Lambda
1-pchisq((fit7$logLik - fit5$logLik) * 2,df = 1)
1-pchisq((fit8$logLik - fit5$logLik) * 2,df = 1)
1-pchisq((fit9$logLik - fit5$logLik) * 2,df = 1)
waldtest(fit5); waldtest(fit7); waldtest(fit8); waldtest(fit9);

#~cityind$Phigh
##MMR
fit5 <- mvmeta(yall~cityind$Phigh,Sall,method = "ml")
aic5 <- (fit5$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) )  %>% round(1)

##EMMR
fit6 <- mixmeta(as.matrix(tmeanpar[,5:8])~cityind$Phigh, tmeanpar[,9:18], data=tmeanpar, random=~1|state/city,method="ml")
aic6 <- summary(fit6)$AIC  %>% round(1)

#MCMAR
system.time(
  fit7 <- smvmeta(yall~cityind$Phigh,S=Sall,Cmatrix = Cmatrix,method= "ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic7 <- (fit7$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

##MSLMAR
system.time(
  fit8 <- Lmvmeta(yall~cityind$Phigh,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic8 <- (fit8$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

##MSEMAR
system.time(
  fit9 <- Emvmeta(yall~cityind$Phigh,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic9 <- (fit9$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

aic5; aic6; aic7; aic8; aic9
fit7$rho; fit8$Rho; fit9$Lambda
1-pchisq((fit7$logLik - fit5$logLik) * 2,df = 1)
1-pchisq((fit8$logLik - fit5$logLik) * 2,df = 1)
1-pchisq((fit9$logLik - fit5$logLik) * 2,df = 1)
waldtest(fit5); waldtest(fit7); waldtest(fit8); waldtest(fit9);


#~cityind$Punem
##MMR
fit5 <- mvmeta(yall~cityind$Punem,Sall,method = "ml")
aic5 <- (fit5$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) )  %>% round(1)

##EMMR
fit6 <- mixmeta(as.matrix(tmeanpar[,5:8])~cityind$Punem, tmeanpar[,9:18], data=tmeanpar, random=~1|state/city,method="ml")
aic6 <- summary(fit6)$AIC  %>% round(1)

#MCMAR
system.time(
  fit7 <- smvmeta(yall~cityind$Punem,S=Sall,Cmatrix = Cmatrix,method= "ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic7 <- (fit7$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

##MSLMAR
system.time(
  fit8 <- Lmvmeta(yall~cityind$Punem,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic8 <- (fit8$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

##MSEMAR
system.time(
  fit9 <- Emvmeta(yall~cityind$Punem,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic9 <- (fit9$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

aic5; aic6; aic7; aic8; aic9
fit7$rho; fit8$Rho; fit9$Lambda
1-pchisq((fit7$logLik - fit5$logLik) * 2,df = 1)
1-pchisq((fit8$logLik - fit5$logLik) * 2,df = 1)
1-pchisq((fit9$logLik - fit5$logLik) * 2,df = 1)
waldtest(fit5); waldtest(fit7); waldtest(fit8); waldtest(fit9);

#~cityind$pop100+cityind$Phigh+cityind$Punem
varnum <- 1+3
dfvarm <- dfvar*varnum

##MMR
fit5 <- mvmeta(yall~cityind$pop100+cityind$Phigh+cityind$Punem,Sall,method = "ml")
aic5 <- (fit5$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) )  %>% round(1)

##EMMR
fit6 <- mixmeta(as.matrix(tmeanpar[,5:8])~cityind$pop100+cityind$Phigh+cityind$Punem, tmeanpar[,9:18], data=tmeanpar, random=~1|state/city,method="ml")
aic6 <- summary(fit6)$AIC  %>% round(1)

#MCMAR
system.time(
  fit7 <- smvmeta(yall~cityind$pop100+cityind$Phigh+cityind$Punem,S=Sall,Cmatrix = Cmatrix,method= "ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic7 <- (fit7$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

##MSLMAR
system.time(
  fit8 <- Lmvmeta(yall~cityind$pop100+cityind$Phigh+cityind$Punem,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic8 <- (fit8$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

##MSEMAR
system.time(
  fit9 <- Emvmeta(yall~cityind$pop100+cityind$Phigh+cityind$Punem,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic9 <- (fit9$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

aic5; aic6; aic7; aic8; aic9
fit7$rho; fit8$Rho; fit9$Lambda
1-pchisq((fit7$logLik - fit5$logLik) * 2,df = 1)
1-pchisq((fit8$logLik - fit5$logLik) * 2,df = 1)
1-pchisq((fit9$logLik - fit5$logLik) * 2,df = 1)

################################################################################
###############Prediction for city-specific association#########################
################################################################################
# plot the map for RR referring to 50
objtem <- c(10,25,75,90)
predbeta0 <- blup(fit0)
predbeta1 <- fit3$fitted.values.spatial
objbasis <- onebasis(tmean[avgtmeansum$perc %in% paste0(objtem, ".0%")], fun="bs", degree=2, knots=knots)
refbasis <- onebasis(tmean[avgtmeansum$perc %in% paste0(50, ".0%")], fun="bs", degree=2, knots=knots)
difx <- t(t(objbasis) - as.numeric(refbasis))
allRRfit0 <- exp(predbeta0 %*% t(difx)) 
allRRfit1 <- exp(predbeta1 %*% t(difx))
allRRfit <- cbind(cityind$city,allRRfit0,allRRfit1) %>% as.data.frame()
names(allRRfit) <- c("city","MMR_"%+%objtem,"MSLMAR_"%+%objtem)

citypoints <- cityind[,1:6]
citypoints <- merge(citypoints,allRRfit,by = "city")
lat <- as.numeric(as.character(cityind$lat))
long <- -as.numeric(as.character(cityind$long))

##10%
par(mfrow=c(1,2))
par(mar=c(1,1,1,1))
cutoff <- pretty(as.matrix(citypoints[,c(7,11)]), n=7)
labmap <- paste(format(cutoff)[-length(cutoff)], format(cutoff)[-1], sep="-")
colmap <- rev(c("#8B0000","#BE0000","#FF1516","#FFC0CB","#ADD8E6","#94B9D9"))
rrcat <- cut(as.numeric(citypoints[,7]), cutoff,  labels=labmap)
map("state", interior=F)
map("state", lty=2, add=T)
points(long, lat, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
rrcat <- cut(as.numeric(citypoints[,11]), cutoff,  labels=labmap)
map("state", interior=F)
map("state", lty=2, add=T)
points(long, lat, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
legend(-80.5962 ,37.88514, labmap, pch=19, col=colmap, bty="n", pt.cex=1.5, cex=0.8,
       title="RR", inset=0.02, xpd=T)

##25%
par(mfrow=c(1,2))
par(mar=c(1,1,1,1))
cutoff <- pretty(as.matrix(citypoints[,c(8,12)]), n=7)
labmap <- paste(format(cutoff)[-length(cutoff)], format(cutoff)[-1], sep="-")
colmap <- rev(c("#FF6A70","#ADD8E6","#94B9D9","#7B9ACC","#627BBF","#4A5CB2","#00008B"))
rrcat <- cut(as.numeric(citypoints[,8]), cutoff,  labels=labmap)
map("state", interior=F)
map("state", lty=2, add=T)
points(long, lat, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
rrcat <- cut(as.numeric(citypoints[,12]), cutoff,  labels=labmap)
map("state", interior=F)
map("state", lty=2, add=T)
points(long, lat, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
legend(-80.5962 ,37.88514, labmap, pch=19, col=colmap, bty="n", pt.cex=1.5, cex=0.8,
       title="RR", inset=0.02, xpd=T)

##75%
par(mfrow=c(1,2))
par(mar=c(1,1,1,1))
cutoff <- pretty(as.matrix(citypoints[,c(9,13)]), n=7)
labmap <- paste(format(cutoff)[-length(cutoff)], format(cutoff)[-1], sep="-")
colmap <- rev(c("#8B0000","#FF1516","#FFC0CB","#ADD8E6","#00008B"))
rrcat <- cut(as.numeric(citypoints[,9]), cutoff,  labels=labmap)
map("state", interior=F)
map("state", lty=2, add=T)
points(long, lat, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
rrcat <- cut(as.numeric(citypoints[,13]), cutoff,  labels=labmap)
map("state", interior=F)
map("state", lty=2, add=T)
points(long, lat, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
legend(-80.5962 ,37.88514, labmap, pch=19, col=colmap, bty="n", pt.cex=1.5, cex=0.8,
       title="RR", inset=0.02, xpd=T)

##90%
par(mfrow=c(1,2))
par(mar=c(1,1,1,1))
cutoff <- pretty(as.matrix(citypoints[,c(10,14)]), n=6)
labmap <- paste(format(cutoff)[-length(cutoff)], format(cutoff)[-1], sep="-")
colmap <- rev(c("#8B0000","#BE0000","#FF1516","#FFC0CB","#4A5CB2"))
rrcat <- cut(as.numeric(citypoints[,10]), cutoff,  labels=labmap)
map("state", interior=F)
map("state", lty=2, add=T)
points(long, lat, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
rrcat <- cut(as.numeric(citypoints[,14]), cutoff,  labels=labmap)
map("state", interior=F)
map("state", lty=2, add=T)
points(long, lat, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
legend(-80.5962 ,37.88514, labmap, pch=19, col=colmap, bty="n", pt.cex=1.5, cex=0.8,
       title="RR", inset=0.02, xpd=T)



################################################################################
# R code for the analysis in:
#
#########An added case example of temperature-HFMD relationships################
################################################################################
`%+%` <- function(x,y) paste0(x,y)
path <- getwd()
source(path%+%"/fun_MCMAR.R")
source(path%+%"/fun_SEM-MVMR.R")
source(path%+%"/fun_SLM-MVMR.R")
source(path%+%"/fun_common.R")

# Wald test for MMR, MCMAR and MSMAR
waldtest <- function(model){
  m <- model
  coef <- m$coefficients[2,]
  vcov <- m$vcov[seq(2,10,2),seq(2,10,2)]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  pvalue <- 1-pchisq(waldstat,df)
  wald <- c(waldstat,df,pvalue)
  wald
}

# Load packages
library(mvmeta);library(sp);library(mixmeta)
# Load data
load(path%+%"/data/stage1data.Rdata")
load(path%+%"/data/mapdata_China_province_project.Rdata")

################################################################################
##### Get the spatially adjacent matrix            #############################
################################################################################
## k-nearest neighbor method-----
k <- 4
Rmethod <- "_"%+%k%+%"n"
if(k == 4) Rmethod <- ""
temp <- CINF
coordinates(temp) <- ~ POINT_X + POINT_Y 
nb <- knearneigh(temp, k=k) 
nb <- knn2nb(nb)
plot(nb,coordinates(temp),pch = 16)
W <- nb2mat(nb)
n <- ncol(W)
C <- W * apply(W, 1, function(x) sum(x!=0))
# get symatric C
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    if(C[i,j]!= C[j,i]){
      C[i,j] <- 1
      C[j,i] <- 1
    } 
  }
}
R <- diag(rowSums(C)) - C
Cmatrix <- diag(n) - R
Wmatrix <- W

################################################################################
############# Compare MMR with MCMAR and MSMAR with only intercept##############
################################################################################
method <- "ml"
xvar <- seq(0,100,by=5)
bvar <- do.call("onebasis",c(list(x=xvar),attr(M.CB,"argvar")))
dfvar <- 5

## MMR
fit0 <- mvmeta(yall,Sall,method = method)
aic0 <- (fit0$logLik * (-2) + 2 * dfvar + dfvar*(dfvar+1)) %>% round(1)

## EMMR
Sallmat <- do.call("rbind",lapply(Sall, function(x) x[lower.tri(x, diag = TRUE)]))
tmeanpar <- cbind(CINF$`Climate region`,CINF$City,yall,Sallmat) %>% as.data.frame
colnames(tmeanpar)[1:2] <- c("Climateregion","City")
tmeanpar[,3:22] <- apply(tmeanpar[,3:22], 2, as.numeric)
fit1 <- mixmeta(as.matrix(tmeanpar[,3:7]), tmeanpar[,8:22], data=tmeanpar, random=~1|Climateregion/City,method="ml")
aic1 <- summary(fit1)$AIC %>% round(1)

## MCMAR
system.time(
  fit2 <- smvmeta(yall,S=Sall,Cmatrix = Cmatrix,method=method,
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 10,opt.iter.show = T))
)
aic2 <- (fit2$logLik * (-2) + 2 * dfvar + dfvar*(dfvar+1) + 2)  %>% round(1)

## MSLMAR
system.time(
  fit3 <- Lmvmeta(yall,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic3 <- (fit3$logLik * (-2) + 2 * (dfvar) + dfvar*(dfvar+1) + 2)  %>% round(1)
fit3coe <- Outoerr_L(fit3)

## MSEMAR
system.time(
  fit4 <- Emvmeta(yall,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic4 <- (fit4$logLik * (-2) + 2 * (dfvar) + dfvar*(dfvar+1) + 2)  %>% round(1)
fit4coe <- Outoerr_E(fit4)

aic0; aic1; aic2; aic3; aic4
fit2$rho; fit3$Rho; fit4$Lambda
1-pchisq((fit2$logLik - fit0$logLik) * 2,df = 1)
1-pchisq((fit3$logLik - fit0$logLik) * 2,df = 1)
1-pchisq((fit4$logLik - fit0$logLik) * 2,df = 1)

# Average temperature-mortality association
# Predict the association
cp <- crosspred(bvar, coef= as.numeric(fit0$coefficients), vcov=fit0$vcov, model.link="log",by=5)
cp1 <- crosspred(bvar, coef= as.numeric(fit2$coefficients), vcov=fit2$vcov, model.link="log",by=5)
cp2 <- crosspred(bvar, coef= as.numeric(fit3coe$coefsim), vcov=fit3coe$vcovsim, model.link="log",by=5)
cp3 <- crosspred(bvar, coef= as.numeric(fit4coe$coefsim), vcov=fit4coe$vcovsim, model.link="log",by=5)

bluetrans <- rgb(0, 0,255, 120, maxColorValue=255) 
par(mfrow=c(1,3))

# MMR vs MCMAR
plot(cp,ylim=c(0.6,1.6), xlab="Temperature percentile",
     ylab="RR",xaxt="n", cex.axis=0.8, col="red",lwd=2)
xperc <- seq(0,100,by=20)
axis(1,at=xperc, labels=paste0(xperc, "%"))
lines(cp1,col="blue",lty=4,lwd=2,ci="area",ci.arg=list(density=20,col=bluetrans))
legend(x="top",inset =0, legend=c("MMR      (AIC = "%+%aic0%+%")", "MCMAR (AIC = "%+%aic2%+%")"),
       lwd=1.5, lty=1, col=c("red", "blue"), bty="n",ncol=1, cex=1.0)

# MMR vs MSLMAR
plot(cp,ylim=c(0.6,1.6), xlab="Temperature percentile",
     ylab="RR",xaxt="n", cex.axis=0.8, col="red",lwd=2)
xperc <- seq(0,100,by=20)
axis(1,at=xperc, labels=paste0(xperc, "%"))
lines(cp2,col="blue",lty=4,lwd=2,ci="area",ci.arg=list(density=20,col=bluetrans))
legend(x="top",inset =0, legend=c("MMR      (AIC = "%+%aic0%+%")", "MSLMAR (AIC = "%+%aic3%+%")"),
       lwd=1.5, lty=1, col=c("red", "blue"), bty="n",ncol=1, cex=1.0)

# MMR vs MSEMAR
plot(cp,ylim=c(0.6,1.6), xlab="Temperature percentile",
     ylab="RR",xaxt="n", cex.axis=0.8, col="red",lwd=2)
xperc <- seq(0,100,by=20)
axis(1,at=xperc, labels=paste0(xperc, "%"))
lines(cp3,col="blue",lty=4,lwd=2,ci="area",ci.arg=list(density=20,col=bluetrans))
legend(x="top",inset =0, legend=c("MMR      (AIC = "%+%aic0%+%")", "MSEMAR (AIC = "%+%aic4%+%")"),
       lwd=1.5, lty=1, col=c("red", "blue"), bty="n",ncol=1, cex=1.0)


################################################################################
######### Compare models of MMR with MCMAR and MSMAR with covariates ###########
################################################################################
varnum <- 1+1
dfvarm <- dfvar*varnum
#~CINF$Altitude
## MMR
fit5 <- mvmeta(yall~CINF$Altitude,Sall,method = "ml")
aic5 <- (fit5$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) ) %>% round(1)

## EMMR
fit6 <- mixmeta(as.matrix(tmeanpar[,3:7])~CINF$Altitude, tmeanpar[,8:22], data=tmeanpar, random=~1|Climateregion/City,method="ml")
aic6 <- summary(fit6)$AIC %>% round(1)

## MCMAR
system.time(
  fit7 <- smvmeta(yall~CINF$Altitude,S=Sall,Cmatrix = Cmatrix,method= "ml",
                control = list(maxiter = 200,factr = 1e7,hessian = T,
                               opt.iter = 1,opt.iter.show = T))
)
aic7 <- (fit7$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

## MSLMAR
system.time(
  fit8 <- Lmvmeta(yall~CINF$Altitude,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic8 <- (fit8$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

## MSEMAR
system.time(
  fit9 <- Emvmeta(yall~CINF$Altitude,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic9 <- (fit9$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

aic5; aic6; aic7; aic8; aic9
fit7$rho; fit8$Rho; fit9$Lambda
1-pchisq((fit7$logLik - fit5$logLik) * 2,df = 1)
1-pchisq((fit8$logLik - fit5$logLik) * 2,df = 1)
1-pchisq((fit9$logLik - fit5$logLik) * 2,df = 1)
waldtest(fit5); waldtest(fit7); waldtest(fit8); waldtest(fit9);

#~CINF$Rainfall
## MMR
fit5 <- mvmeta(yall~CINF$Rainfall,Sall,method = "ml")
aic5 <- (fit5$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) ) %>% round(1)

## EMMR
fit6 <- mixmeta(as.matrix(tmeanpar[,3:7])~CINF$Rainfall, tmeanpar[,8:22], data=tmeanpar, random=~1|Climateregion/City,method="ml")
aic6 <- summary(fit6)$AIC %>% round(1)

## MCMAR
system.time(
  fit7 <- smvmeta(yall~CINF$Rainfall,S=Sall,Cmatrix = Cmatrix,method= "ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic7 <- (fit7$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

## MSLMAR
system.time(
  fit8 <- Lmvmeta(yall~CINF$Rainfall,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic8 <- (fit8$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

## MSEMAR
system.time(
  fit9 <- Emvmeta(yall~CINF$Rainfall,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic9 <- (fit9$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

aic5; aic6; aic7; aic8; aic9
fit7$rho; fit8$Rho; fit9$Lambda
1-pchisq((fit7$logLik - fit5$logLik) * 2,df = 1)
1-pchisq((fit8$logLik - fit5$logLik) * 2,df = 1)
1-pchisq((fit9$logLik - fit5$logLik) * 2,df = 1)
waldtest(fit5); waldtest(fit7); waldtest(fit8); waldtest(fit9);

#~CINF$Altitude+CINF$Rainfall
varnum <- 1+2
dfvarm <- dfvar*varnum
## MMR
fit5 <- mvmeta(yall~CINF$Altitude+CINF$Rainfall,Sall,method = "ml")
aic5 <- (fit5$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) ) %>% round(1)

## EMMR
fit6 <- mixmeta(as.matrix(tmeanpar[,3:7])~CINF$Altitude+CINF$Rainfall, tmeanpar[,8:22], data=tmeanpar, random=~1|Climateregion/City,method="ml")
aic6 <- summary(fit6)$AIC %>% round(1)

## MCMAR
system.time(
  fit7 <- smvmeta(yall~CINF$Altitude+CINF$Rainfall,S=Sall,Cmatrix = Cmatrix,method= "ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic7 <- (fit7$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

## MSLMAR
system.time(
  fit8 <- Lmvmeta(yall~CINF$Altitude+CINF$Rainfall,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic8 <- (fit8$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

## MSEMAR
system.time(
  fit9 <- Emvmeta(yall~CINF$Altitude+CINF$Rainfall,S=Sall,Wmatrix = Wmatrix,method="ml",
                  control = list(maxiter = 200,factr = 1e7,hessian = T,
                                 opt.iter = 1,opt.iter.show = T))
)
aic9 <- (fit9$logLik * (-2) + 2 * (dfvarm) + dfvar*(dfvar+1) + 2)  %>% round(1)

aic5; aic6; aic7; aic8; aic9
fit7$rho; fit8$Rho; fit9$Lambda
1-pchisq((fit7$logLik - fit5$logLik) * 2,df = 1)
1-pchisq((fit8$logLik - fit5$logLik) * 2,df = 1)
1-pchisq((fit9$logLik - fit5$logLik) * 2,df = 1)


################################################################################
###############Prediction for city-specific association#########################
################################################################################
# plot the map for RR referring to 50
objtem <- c(10,25,75,90)
predbeta0 <- blup(fit0)
predbeta1 <- fit3$fitted.values.spatial
objbasis <- do.call("onebasis",c(list(x=objtem),attr(M.CB,"argvar"))) 
refbasis <- do.call("onebasis",c(list(x=50),attr(M.CB,"argvar")))
difx <- t(t(objbasis) - as.numeric(refbasis))
allRRfit0 <- exp(predbeta0 %*% t(difx)) 
allRRfit1 <- exp(predbeta1 %*% t(difx))
allRRfit <- cbind(CINF$citycd,allRRfit0,allRRfit1) %>% as.data.frame()
names(allRRfit) <- c("citycd","MMR_"%+%objtem,"MSLMAR_"%+%objtem)

citypoints <- CINF[,1:6]
citypoints <- merge(citypoints,allRRfit,by = "citycd")

##10%
par(mar = c(1, 1, 1, 1))
par(mfrow=c(1,2))
cutoff <- pretty(as.matrix(citypoints[,c(7,11)]), n=7)
labmap <- paste(format(cutoff)[-length(cutoff)], format(cutoff)[-1], sep="-")
colmap <- rev(c("#8B0000","#FF1516","#FFC0CB","#ADD8E6","#94B9D9","#7B9ACC","#627BBF","#4A5CB2","#313DA4","#00008B"))
rrcat <- cut(as.numeric(citypoints[,7]), cutoff,  labels=labmap)
plot(mapdata,lty=2)
points(citypoints$POINT_X, citypoints$POINT_Y, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
rrcat <- cut(as.numeric(citypoints[,11]), cutoff,  labels=labmap)
plot(mapdata,lty=2)
points(citypoints$POINT_X, citypoints$POINT_Y, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
legend(-3300000, 3589877, labmap, pch=19, col=colmap, bty="n", pt.cex=1.5, cex=0.8,
       title="RR", inset=0.02, xpd=T)

##25%
par(mar = c(1, 1, 1, 1))
par(mfrow=c(1,2))
cutoff <- pretty(as.matrix(citypoints[,c(8,12)]), n=7)
labmap <- paste(format(cutoff)[-length(cutoff)], format(cutoff)[-1], sep="-")
colmap <- rev(c("#8B0000","#FF1516","#FF6A70","#FFE4E1","#ADD8E6","#94B9D9","#7B9ACC","#4A5CB2","#313DA4","#00008B"))
rrcat <- cut(as.numeric(citypoints[,8]), cutoff,  labels=labmap)
plot(mapdata,lty=2)
points(citypoints$POINT_X, citypoints$POINT_Y, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
rrcat <- cut(as.numeric(citypoints[,12]), cutoff,  labels=labmap)
plot(mapdata,lty=2)
points(citypoints$POINT_X, citypoints$POINT_Y, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
legend(-3300000, 3589877, labmap, pch=19, col=colmap, bty="n", pt.cex=1.5, cex=0.8,
       title="RR", inset=0.02, xpd=T)

##75%
par(mar = c(1, 1, 1, 1))
par(mfrow=c(1,2))
cutoff <- pretty(as.matrix(citypoints[,c(9,13)]), n=7)
labmap <- paste(format(cutoff)[-length(cutoff)], format(cutoff)[-1], sep="-")
colmap <- rev(c("#8B0000","#FF1516","#FF6A70","#FFC0CB","#FFE4E1","#ADD8E6","#7B9ACC","#4A5CB2","#00008B"))
rrcat <- cut(as.numeric(citypoints[,9]), cutoff,  labels=labmap)
plot(mapdata,lty=2)
points(citypoints$POINT_X, citypoints$POINT_Y, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
rrcat <- cut(as.numeric(citypoints[,13]), cutoff,  labels=labmap)
plot(mapdata,lty=2)
points(citypoints$POINT_X, citypoints$POINT_Y, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
legend(-3300000, 3589877, labmap, pch=19, col=colmap, bty="n", pt.cex=1.5, cex=0.8,
       title="RR", inset=0.02, xpd=T)

##90%
par(mar = c(1, 1, 1, 1))
par(mfrow=c(1,2))
cutoff <- pretty(as.matrix(citypoints[,c(10,14)]), n=7)
labmap <- paste(format(cutoff)[-length(cutoff)], format(cutoff)[-1], sep="-")
colmap <- rev(c("#8B0000","#BE0000","#FF1516","#FF6A70","#FFC0CB","#FFE4E1","#ADD8E6","#00008B"))
rrcat <- cut(as.numeric(citypoints[,10]), cutoff,  labels=labmap)
plot(mapdata,lty=2)
points(citypoints$POINT_X, citypoints$POINT_Y, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
rrcat <- cut(as.numeric(citypoints[,14]), cutoff,  labels=labmap)
plot(mapdata,lty=2)
points(citypoints$POINT_X, citypoints$POINT_Y, col=alpha(colmap[rrcat], 0.6), pch=19, cex=0.8)
legend(-3300000, 3589877, labmap, pch=19, col=colmap, bty="n", pt.cex=1.5, cex=0.8,
       title="RR", inset=0.02, xpd=T)
