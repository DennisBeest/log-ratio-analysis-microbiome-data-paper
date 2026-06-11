

rm(list=ls(all=TRUE))

gitdir <- "E:/1Documenten/Werk/gits/gitlabwur/projects/3microbiome/LRApaper/results def"
setwd(gitdir)

#---load file with functions
source("functions def.R")



################################################################################################
#--load data
################################################################################################



load("dataMidges.RData")
load("dataRice.RData")
ls()



################################################################################################
#---Generate data
################################################################################################


set.seed(321)
outx3 <- dataSim(sdSite=1, sdSpecies=2, weightResp=3, drawDistr="nb", Max=1, strengthCorrResp=0, zeroThres=5)



################################################################################################



#savepics <- TRUE
savepics <- FALSE
saveloc <- paste(gitdit,"/CombinedPics", sep="")
saveInfo <- list()
saveInfo[[1]] <- savepics
saveInfo[[2]] <- saveloc



################################################################################################


ss1 <- colSums(dataRice$x!=0) >= 150 & colSums(dataRice$x) >= 2000
newxRice <- dataRice$x[, ss1]

ss2 <- colSums(dataMidges$x!=0) >= 75 & colSums(dataMidges$x) >= 2000
newxKnut <- dataMidges$x[, ss2]

index3 <- colSums(outx3$xs > 0) >= 30 & colSums(outx3$xs) >= 500
newx3 <- outx3$xs[, index3]


################################################################################################


#dim(newxRice)
#sum(newxRice==0)/length(newxRice)
 
#dim(newxKnut)
#sum(newxKnut ==0)/length(newxKnut)

#dim(newx3)
#sum(newx3==0)/length(newx3)

#dim(dataRice$x)
#dim(dataMidges$x)

#1-sum(newx3)/sum(outx3$xs)
#[1] 0.1366326
#1-sum(newxKnut)/sum(dataMidges$x)
#[1] 0.04328441
#1-sum(newxRice)/sum(dataRice$x)
#[1] 0.004759742



################################################################################################
#---LRA filtered
################################################################################################


saveInfo[[3]] <- "LRAFilt.pdf"

outRice <- LRA.ade4(x=newxRice, response=dataRice$treatment, pseudocount=1)
outKnut <- LRA.ade4(x=newxKnut, response=dataMidges$response, pseudocount=1)
outSim <- LRA.ade4(x=newx3, response=outx3$response)

bothPics2(lraRice=outRice, lraKnut=outKnut, lraSim=outSim, saveInfo=saveInfo, OrientOrdi = c(1, 1, -1, 1, -1, -1), textAnalysis="Log-ratio PCA & filtered", hjustTitle2=1)


rm(outRice, outKnut, outSim)



###########################################################################################
#----Zero imputation
###########################################################################################



saveInfo[[3]] <- "LRAImpFrac.pdf"


riceImp <- as.matrix(cmultRepl(dataRice$x/rowSums(dataRice$x)))
knutImp <- as.matrix(cmultRepl(dataMidges$x/rowSums(dataMidges$x)))
simImp <- as.matrix(cmultRepl(outx3$xs/rowSums(outx3$xs)))

outRice <- LRA.ade4(x=riceImp, response=dataRice$treatment, pseudocount=0, xReal=dataRice$x)
outKnut <- LRA.ade4(x=knutImp, response=dataMidges$response, pseudocount=0, xReal=dataMidges$x)
outSim <- LRA.ade4(simImp, outx3$response, pseudocount=0, xReal=outx3$xs)


bothPics2(outRice, outKnut, outSim , saveInfo=saveInfo, OrientOrdi = c(1, 1, -1, -1, 1, -1), textAnalysis="Log-ratio PCA & GBM", hjustTitle2=1.2)


rm(outRice, outKnut, outSim)



################################################################################################
#---Weighted LRA
################################################################################################



saveInfo[[3]] <- "WLRA.pdf"


outRice <- wLRA.ade4(x=dataRice$x, response=dataRice$treatment)
outKnut <- wLRA.ade4(x=dataMidges$x, response=dataMidges$response)
outSim <- wLRA.ade4(outx3$xs, outx3$response)


bothPics2(outRice, outKnut, outSim, saveInfo=saveInfo, OrientOrdi = c(1, 1, 1, 1, 1, 1), textAnalysis="Weigthed log-ratio PCA", hjustTitle2=1)


rm(outRice, outKnut, outSim)



################################################################################################
#---Weighted LRA filtered
################################################################################################



saveInfo[[3]] <- "WLRAFilt.pdf"


outRice <- wLRA.ade4(x=newxRice, response=dataRice$treatment)
outKnut <- wLRA.ade4(x=newxKnut, response=dataMidges$response)
outSim <- wLRA.ade4(newx3, outx3$response)


bothPics2(outRice, outKnut, outSim, saveInfo=saveInfo, OrientOrdi = c(-1, 1, -1, 1, -1, 1), textAnalysis="Weigthed log-ratio PCA & filtered", hjustTitle2=1)


rm(outRice, outKnut, outSim)



################################################################################################
#---LogFrac
################################################################################################



saveInfo[[3]] <- "LogFrac.pdf"


outRice <- LogFrac.ade4(x=dataRice$x, response=dataRice$treatment, pseudoFrac=0.001)
outKnut <- LogFrac.ade4(x=dataMidges$x, response=dataMidges$response, pseudoFrac=0.001)
outSim <- LogFrac.ade4(x=outx3$xs, response=outx3$response, pseudoFrac=0.001)


bothPics1(outRice, outKnut, outSim, saveInfo=saveInfo, OrientOrdi = c(-1, 1, -1, 1, 1, 1), textAnalysis="Log proportions PCA", hjustTitle2=0.9)


rm(outRice, outKnut, outSim)



################################################################################################
#---LogFrac filtered
################################################################################################



saveInfo[[3]] <- "LogFracFilt.pdf"


outRice <- LogFrac.ade4(x=newxRice, response=dataRice$treatment, pseudoFrac=0.001)
outKnut <- LogFrac.ade4(x=newxKnut, response=dataMidges$response, pseudoFrac=0.001)
outSim <- LogFrac.ade4(newx3, outx3$response, pseudoFrac=0.001)


bothPics1(outRice, outKnut, outSim, saveInfo=saveInfo, OrientOrdi = c(-1, 1, -1, -1, 1, -1), textAnalysis="Log prop. PCA & filtered", hjustTitle2=0.7)


rm(outRice, outKnut, outSim)



################################################################################################
#---CA untransformed
################################################################################################



saveInfo[[3]] <- "CAUntrans.pdf"


outRice <- CA.ade4(x=dataRice$x, response=dataRice$treatment, xReal=dataRice$x)
outKnut <- CA.ade4(x=dataMidges$x, dataMidges$response, xReal=dataMidges$x)
outSim <- CA.ade4(x=outx3$xs, response=outx3$response, xReal=outx3$xs)


bothPics2(outRice, outKnut, outSim, saveInfo=saveInfo, OrientOrdi = c(1, -1, 1, -1, -1, -1), textAnalysis="CA-counts", hjustTitle2=0.7)


rm(outRice, outKnut, outSim)



################################################################################################
#---CA untransformed
################################################################################################


saveInfo[[3]] <- "CAUntransFilt.pdf"


outRice <- CA.ade4(x=newxRice, response=dataRice$treatment, xReal=newxRice)
outKnut <- CA.ade4(x=newxKnut, dataMidges$response, xReal=newxKnut)
outSim <- CA.ade4(x=newx3, response=outx3$response, xReal=newx3)


bothPics2(outRice, outKnut, outSim, saveInfo=saveInfo, OrientOrdi = c(-1, -1, 1, -1, -1, -1), textAnalysis="CA-counts & filtered", hjustTitle2=1.2)


################################################################################################
#---CA Sqrt
################################################################################################


saveInfo[[3]] <- "CASqrt.pdf"


outRice <- CA.ade4(x=sqrt(dataRice$x), response=dataRice$treatment, xReal=dataRice$x)
outKnut <- CA.ade4(x=sqrt(dataMidges$x), dataMidges$response, xReal=dataMidges$x)
outSim <- CA.ade4(x=sqrt(outx3$xs), response=outx3$response, xReal=outx3$xs)


bothPics2(outRice, outKnut, outSim, saveInfo=saveInfo, OrientOrdi = c(1, -1, -1, -1, 1, -1), textAnalysis="CA-sqrt", hjustTitle2=0.9)


rm(outRice, outKnut, outSim)



################################################################################################
#---CA Sqrt filtered
################################################################################################



saveInfo[[3]] <- "CASqrtFilt.pdf"


outRice <- CA.ade4(x=sqrt(newxRice), response=dataRice$treatment, xReal=newxRice)
outKnut <- CA.ade4(x=sqrt(newxKnut), dataMidges$response, xReal=newxKnut)
outSim <- CA.ade4(x=sqrt(newx3), response=outx3$response, xReal=newx3)


bothPics2(outRice, outKnut, outSim, saveInfo=saveInfo, OrientOrdi = c(-1, 1, -1, -1, 1, -1), textAnalysis="CA-sqrt & filtered", hjustTitle2=0.5)


rm(outRice, outKnut, outSim)


################################################################################################
#---CA log
################################################################################################


saveInfo[[3]] <- "CALog.pdf"


outRice <- CA.ade4(x=log(dataRice$x+1), response=dataRice$treatment, xReal=dataRice$x)
outKnut <- CA.ade4(x=log(dataMidges$x+1), dataMidges$response, xReal=dataMidges$x)
outSim <- CA.ade4(x=log(outx3$xs+1), response=outx3$response, xReal=outx3$xs)

bothPics2(outRice, outKnut, outSim, saveInfo=saveInfo, OrientOrdi = c(-1, 1, -1, -1, 1, -1), textAnalysis="CA-log", hjustTitle2=0.9)

rm(outRice, outKnut, outSim)



################################################################################################
#---CA log filtered
################################################################################################


saveInfo[[3]] <- "CALogFilt.pdf"


outRice <- CA.ade4(x=log(newxRice + 1), response=dataRice$treatment, xReal=newxRice)
outKnut <- CA.ade4(x=log(newxKnut + 1), dataMidges$response, xReal=newxKnut)
outSim <- CA.ade4(x=log(newx3 + 1), response=outx3$response, xReal=newx3)


bothPics2(outRice, outKnut, outSim, saveInfo=saveInfo, OrientOrdi = c(1, 1, -1, -1, 1, -1), textAnalysis="CA-log & filtered", hjustTitle2=0.3)


rm(outRice, outKnut, outSim)



################################################################################################
#---RCM
################################################################################################


saveInfo[[3]] <- "RCMfig.pdf"
runRCM <- FALSE

if(runRCM)
{
	#----------------------------------------------------

	simRCM <- RCM(outx3$xs, k=2, prevCutOff=0, minFraction=0)
	simRCM2 <- RCM(outx3$xs, k=2)
	simRCM3 <- RCM(newx3, k=2)

	#----------------------------------------------------

	knutRCM <- RCM(dataMidges$x, k=2, prevCutOff=0, minFraction=0)
	knutRCM2 <- RCM(dataMidges$x, k=2)
	knutRCM3 <- RCM(newxKnut, k=2)

	#----------------------------------------------------

	riceRCM <- RCM(dataRice$x, k=2, prevCutOff=0, minFraction=0)
	riceRCM2 <- RCM(dataRice$x, k=2)
	riceRCM3 <- RCM(newxRice, k=2)

	setwd(gitdit)
	save(simRCM, simRCM2, simRCM3, knutRCM, knutRCM2, knutRCM3, riceRCM, riceRCM2, riceRCM3, file="CombinedPics/rcmObjects.RData" )
}

if(!runRCM)
{
	setwd(gitdit)
	load("CombinedPics/rcmObjects.RData" )
}


#dim(simRCM$cMat)
#dim(simRCM2$cMat)
#dim(simRCM3$cMat)
#dim(knutRCM$cMat)
#dim(knutRCM2$cMat)
#dim(knutRCM3$cMat)
#dim(riceRCM$cMat)
#dim(riceRCM2$cMat)
#dim(riceRCM3$cMat)


outRice1=outRice2=outRice3 <- c()
outRice1$Sc <- cbind(riceRCM$rMat[, 1], riceRCM$rMat[, 2])
outRice2$Sc <- cbind(riceRCM2$rMat[, 1], riceRCM2$rMat[, 2])
outRice3$Sc <- cbind(riceRCM3$rMat[, 1], riceRCM3$rMat[, 2])
outRice1$resp=outRice2$resp=outRice3$resp <- dataRice$treatment

outKnut1=outKnut2=outKnut3 <- c()
outKnut1$Sc <- cbind(knutRCM$rMat[, 1], knutRCM$rMat[, 2])
outKnut2$Sc <- cbind(knutRCM2$rMat[, 1], knutRCM2$rMat[, 2])
outKnut3$Sc <- cbind(knutRCM3$rMat[, 1], knutRCM3$rMat[, 2])
outKnut1$resp=outKnut2$resp=outKnut3$resp <- dataMidges$response

outSim1=outSim2=outSim3 <- c()
outSim1$Sc <- cbind(simRCM$rMat[, 1], simRCM$rMat[, 2])
outSim2$Sc <- cbind(simRCM2$rMat[, 1], simRCM2$rMat[, 2])
outSim3$Sc <- cbind(simRCM3$rMat[, 1], simRCM3$rMat[, 2])
outSim1$resp=outSim2$resp=outSim3$resp <- outx3$response


picRCM(outSim1, outSim3, outRice1, outRice3, outKnut1, outKnut3, saveInfo=saveInfo, OrientOrdi = c(-1, 1, -1, 1, 1, 1, -1, 1, -1, 1, -1, 1))




