
rm(list=ls(all=TRUE))
#---load functions


gitdir <- "E:/1Documenten/Werk/gits/gitlabwur/projects/3microbiome/LRApaper/results def"
#gitdir <- "C:/1Files/2Git/smallProjects/microbiome/LRApaper/results def"
#gitdir <- "/home/dennis/Documents/gitwur/smallProjects/microbiome/LRApaper/results def"
setwd(gitdir)


#---load file with functions
source("functions def.R")


################################################################################################
#---load data
################################################################################################


load("dataMidges.RData")
load("dataRice.RData")
ls()


################################################################################################


#dim(dataMidges$x)
#dim(dataRice$x)

#range(colSums(dataMidges$x>0))
#range(colSums(dataRice$x>0))

#sum(dataMidges$x==0)/length(dataMidges$x)
#sum(dataRice$x==0)/length(dataRice$x)

#range(rowSums(dataMidges$x))
#range(rowSums(dataRice$x))

#Xtr1 <- DoubleCenter(log(dataRice$x + 1))
#cor(as.numeric(dataRice$treatment), attributes(Xtr1)$rowmeans)
# -0.39765
#Xtr2 <- DoubleCenter(log(dataMidges$x + 1))
#cor(as.numeric(dataMidges$response), attributes(Xtr2)$rowmeans)
#0.5440085



################################################################################################
#---Generate data
################################################################################################


sdSite <- c(0, 0.5, 1)
setseed <- 321
set.seed(setseed)
outx1 <- dataSim(sdSite=sdSite[1], sdSpecies=2, weightResp=3, drawDistr="nb", Max=1, strengthCorrResp=0, zeroThres=5)
set.seed(setseed)
outx2 <- dataSim(sdSite=sdSite[2], sdSpecies=2, weightResp=3, drawDistr="nb", Max=1, strengthCorrResp=0, zeroThres=5)
set.seed(setseed)
outx3 <- dataSim(sdSite=sdSite[3], sdSpecies=2, weightResp=3, drawDistr="nb", Max=1, strengthCorrResp=0, zeroThres=5)


saveInfo <- list()
#saveInfo[[1]] <- TRUE 		#sets whether or not pics are saved
saveInfo[[1]] <- FALSE
saveInfo[[2]] <- paste(gitdir,"/mainFigs", sep="")



#dim(outx1$xs)
#dim(outx2$xs)
#dim(outx3$xs)

#range(rowSums(outx1$xs))
#range(rowSums(outx2$xs))
#range(rowSums(outx3$xs))

#sum(outx1$xs==0)/length(outx1$xs)
#sum(outx2$xs==0)/length(outx2$xs)
#sum(outx3$xs==0)/length(outx3$xs)



########################################################################
#---LRA picture, ordination, contribution, correlation
########################################################################



saveInfo[[3]] <- "dataLRACombi.pdf"


lraRice <- LRA.ade4(x=dataRice$x, response=dataRice$treatment, pseudocount=1)
lraKnut <- LRA.ade4(x=dataMidges$x, response=dataMidges$response, pseudocount=1)

dataPics4(lraRice, lraKnut, saveInfo=saveInfo)



################################################################
#----------Log-ratio PCA (LRA)
################################################################
  
  
saveInfo[[3]] <- "SimLRA.pdf"

res1 <- LRA.ade4(outx1$xs, outx1$response)
res2 <- LRA.ade4(outx2$xs, outx2$response)
res3 <- LRA.ade4(outx3$xs, outx3$response)
  
  
simPic1(res1, res2, res3, saveInfo, OrientOrdi = c(1, 1, 1, -1, -1, 1))
  



