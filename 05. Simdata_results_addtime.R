##############################################################################-#
## The following codes replicate our results of the computation time in       #
## simulation study                                                            #
## Author: Caiying Luo                                                         #
##############################################################################-#


################################################################################
##### Get performance for simulation study######################################
################################################################################

########## Load function and data -----
`%+%` <- function(x,y) paste0(x,y)
load("data\\simdata.Rdata")
method <- "ml"
scennames <- names(simdata)
########## get performance for specific scenario ------
finaltime <- matrix(NA,length(scennames),8)
n=0
for (scen in scennames) {
    cat(scen,"")
    n=n+1
    load(getwd()%+%"/data/Simdataresults_time/simrestime_"%+%scen%+%"_"%+%method%+%".Rdata")
    finaltime[n,] <- apply(res_time,2,mean)
}
rownames(finaltime) <- scennames
colnames(finaltime) <- colnames(res_time)
finaltime <- finaltime[,c(1,3,5,7,2,4,6,8)]
finaltime <- round(finaltime,3)

xlsx::write.xlsx(finaltime,file = "performancetime-results-"%+%method%+%".xlsx",sheetName = "time")