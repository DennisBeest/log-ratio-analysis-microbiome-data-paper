
rm(list=ls(all=TRUE))
gitdit <- "E:/1Documenten/Werk/gits/gitlabwur/smallProjects/microbiome/LRApaper/results def"
setwd(gitdit)
#---load file with functions
source("functions def.R")


###############################################################################


nsim <- 2000
#nsim <- 10

sdSiteVec <- c(0.25, 0.5, 1)
sdLength <- length(sdSiteVec)
nScen <- 3


###############################################################################


#--this sets the sd. Repeat this for j = 1, 2, 3.
j <- 1


###############################################################################


savecorr <- matrix(nrow=nsim, ncol=3)
outxsAll <- list()
outxsAll[[1]] <- list()
outxsAll[[2]] <- list()
outxsAll[[3]] <- list()

set.seed(1122)
for(i in 1:nsim)
{
	#i <- 10
	if(i %% 10 == 0) print(i)

	################################################
	
	outxsAll[[1]][[i]] <- dataSim(sdSite=sdSiteVec[j], sdSpecies=2, weightResp=1, drawDistr="nb", Max=1, strengthCorrResp=0, zeroThres=5)
	outxsAll[[2]][[i]] <- dataSim(sdSite=sdSiteVec[j], sdSpecies=2, weightResp=1, drawDistr="nb", Max=1, strengthCorrResp=2, zeroThres=5)
	outxsAll[[3]][[i]] <- dataSim(sdSite=sdSiteVec[j], sdSpecies=2, weightResp=1.5, drawDistr="nb", Max=1, strengthCorrResp=0, zeroThres=5)

	for(k in 1:3)
	{
		y <- outxsAll[[k]][[i]]$resp
		x <- outxsAll[[k]][[i]]$x
		savecorr[i, k] <- cor(rowMeans(log(x + 1)), y)
	}
}




###############################################################################



pvalsLRA <- array(dim = c(sdLength, nScen, nsim))
pvalsLogFrac <- array(dim = c(sdLength, nScen, nsim))
pvalsWLRA <- array(dim = c(sdLength, nScen, nsim))
pvalsCA <- array(dim = c(sdLength, nScen, nsim))
pvalsCASqrt <- array(dim = c(sdLength, nScen, nsim))
pvalsCALog <- array(dim = c(sdLength, nScen, nsim))
pvalsAnosim <- array(dim = c(sdLength, nScen, nsim))
pvalsPermanova <- array(dim = c(sdLength, nScen, nsim))
pvalsLRA.ZI <- array(dim = c(sdLength, nScen, nsim))
pvalsLRA.ZI2 <- array(dim = c(sdLength, nScen, nsim))



resLRA <- matrix(nrow=sdLength, ncol=nScen)
resLogFrac <- matrix(nrow=sdLength, ncol=nScen)
resWLRA <- matrix(nrow=sdLength, ncol=nScen)
resCA <- matrix(nrow=sdLength, ncol=nScen)
resCASqrt <- matrix(nrow=sdLength, ncol=nScen)
resCALog <- matrix(nrow=sdLength, ncol=nScen)
resAnosim <- matrix(nrow=sdLength, ncol=nScen)
resPermanova <- matrix(nrow=sdLength, ncol=nScen)
resLRA.ZI <- matrix(nrow=sdLength, ncol=nScen)
resLRA.ZI2 <- matrix(nrow=sdLength, ncol=nScen)


choices <- 1:10


print("j")
print(sdSiteVec[j])



for(i in 1:nsim)
{
	#i <- 1

	if(i %% 2 == 0) print(i)

	for(k in 1:3)  #--3 scenarios
	{
		#k <- 1

		X <- outxsAll[[k]][[i]]$x
		resp <- outxsAll[[k]][[i]]$response


		################################################
		#---LRA
		################################################


		if(any(choices == 1))
		{
			Xtr <- DoubleCenter(log(X + 1))

			pca1 <- dudi.pca(as.data.frame(Xtr), scannf = FALSE, scale = FALSE, center = FALSE)
			rda1 <- pcaiv(dudi=pca1, df=resp, scannf = FALSE)
			rt <- randtest(rda1, nrepet = 999)
			pvalsLRA[j, k, i] <- rt$pvalue
		}


		################################################
		#---LogFrac
		################################################


		if(any(choices == 2))
		{
			Xtr <- log(X/rowSums(X) + 0.001)
			pca1 <- dudi.pca(as.data.frame(Xtr), scannf = FALSE, scale = FALSE, center = TRUE)
			rda1 <- pcaiv(dudi=pca1, df=resp, scannf = FALSE)
			rt <- randtest(rda1, nrepet = 999)

			pvalsLogFrac[j, k, i] <- rt$pvalue
		}


		################################################
		#---WLRA
		################################################


		if(any(choices == 3))
		{
			xlog <- log(X + 1)
			RowWeights <- rowSums(xlog)
			ColWeights <- colSums(xlog)
			Xtr <- DoubleCenter(xlog, RowWeights=RowWeights, ColWeights=ColWeights)
			pca1 <- dudi.pca(as.data.frame(Xtr), scannf = FALSE, scale = FALSE, center = FALSE, row.w = RowWeights, col.w = ColWeights)
			rda1 <- pcaiv(dudi=pca1, df=resp, scannf = FALSE)

			test <- randtest(rda1, nrepet = 999)
			pvalsWLRA[j, k, i] <- test$pvalue
		}


		################################################
		#---CA Untransformed
		################################################


		if(any(choices == 4))
		{
			Xtr <- X
			res.ca <- dudi.coa(Xtr, scannf = FALSE)
			cca1 <- pcaiv(dudi=res.ca, df=resp, scannf = FALSE)
			an <- randtest(cca1, nrepet = 999)

			pvalsCA[j, k, i] <- an$pvalue
		}


		################################################
		#---CA Sqrt
		################################################


		if(any(choices == 5))
		{
			Xtr <- sqrt(X)
			res.ca <- dudi.coa(Xtr, scannf = FALSE)
			cca1 <- pcaiv(dudi=res.ca, df=resp, scannf = FALSE)
			an <- randtest(cca1)
			pvalsCASqrt[j, k, i] <- an$pvalue
		}


		################################################
		#---CA log
		################################################


		if(any(choices == 6))
		{
			Xtr <- log(X + 1)
			res.ca <- dudi.coa(Xtr, scannf = FALSE)
			cca1 <- pcaiv(dudi=res.ca, df=resp, scannf = FALSE)
			an <- randtest(cca1)

			pvalsCALog[j, k, i] <- an$pvalue
		}


		################################################
		#---Anosim
		################################################


		if(any(choices == 7))
		{
			ans <- anosim(X/rowSums(X), resp)
			pvalsAnosim[j, k, i] <- ans$signif
		}


		################################################
		#---Permanova
		################################################


		if(any(choices == 8))
		{
			ans <- adonis(X/rowSums(X) ~ resp)
			pvalsPermanova[j, k, i] <- ans$aov.tab[1, 6]
		}


		################################################
		#---Zeroreplacement
		################################################


		if(any(choices == 9))
		{
			X.Imp <- as.matrix(cmultRepl(X))
			Xtr <- DoubleCenter(log(X.Imp))
			pca1 <- dudi.pca(as.data.frame(Xtr), scannf = FALSE, scale = FALSE, center = FALSE)
			rda1 <- pcaiv(dudi=pca1, df=resp, scannf = FALSE)
			rt <- randtest(rda1, nrepet = 999)
			pvalsLRA.ZI[j, k, i] <- rt$pvalue
		}

		if(any(choices == 10))
		{
			fracs <- X/rowSums(X)
			X.Imp <- as.matrix(cmultRepl(fracs))
			Xtr <- DoubleCenter(log(X.Imp))
			pca1 <- dudi.pca(as.data.frame(Xtr), scannf = FALSE, scale = FALSE, center = FALSE)
			rda1 <- pcaiv(dudi=pca1, df=resp, scannf = FALSE)
			rt <- randtest(rda1, nrepet = 999)
			pvalsLRA.ZI2[j, k, i] <- rt$pvalue
		}
	}
}


for(k in 1:nScen)
{
	resLRA[j, k] <- mean(pvalsLRA[j, k, ] < 0.05)
	resLogFrac[j, k] <- mean(pvalsLogFrac[j, k, ] < 0.05)
	resWLRA[j, k] <- mean(pvalsWLRA[j, k, ] < 0.05)

	resCA[j, k] <- mean(pvalsCA[j, k, ] < 0.05)
	resCASqrt[j, k] <- mean(pvalsCASqrt[j, k, ] < 0.05)
	resCALog[j, k] <- mean(pvalsCALog[j, k, ] < 0.05)

	resAnosim[j, k] <- mean(pvalsAnosim[j, k, ] < 0.05)
	resPermanova[j, k] <- mean(pvalsPermanova[j, k, ] < 0.05)

	resLRA.ZI[j, k] <- mean(pvalsLRA.ZI[j, k, ] < 0.05)
	resLRA.ZI2[j, k] <- mean(pvalsLRA.ZI2[j, k, ] < 0.05)
}


outres <- rbind(resLRA[j, ], resLogFrac[j, ], resWLRA[j, ], resCA[j, ], resCASqrt[j, ], resCALog[j, ], resAnosim[j, ], resPermanova[j, ], resLRA.ZI[j, ], resLRA.ZI2[j, ])
rownames(outres) <- paste(rep(c("resLRA", "resLogFrac", "resWLRA", "resCA", "resCASqrt", "resCALog", "Anosim", "Permanova", "LRA.ZI", "LRA.ZI2"), each=1), rownames(outres))
colnames(outres) <- c("type1Base", "type1Cor", "power")


outres



###########################################################
#---Result 25 september
###########################################################


#----------------------------------------
2000 runs
#----------------------------------------
j <- 1
            type1Base type1Cor  power
resLRA         0.0525   0.0890 0.7955
resLogFrac     0.0475   0.0525 0.6920
resWLRA        0.0440   0.0545 0.5665
resCA          0.0405   0.0480 0.3475
resCASqrt      0.0390   0.0340 0.7825
resCALog       0.0370   0.0510 0.7285
Anosim         0.0475   0.0510 0.2900
Permanova      0.0460   0.0540 0.2995
LRA.ZI         0.0420   0.0530 0.7835
LRA.ZI2        0.0460   0.0435 0.7910


#----------------------------------------
j <- 2
            type1Base type1Cor  power
resLRA         0.0540   0.6685 0.6755
resLogFrac     0.0530   0.0445 0.6785
resWLRA        0.0530   0.1250 0.5315
resCA          0.0630   0.0260 0.3025
resCASqrt      0.0430   0.0595 0.7500
resCALog       0.0405   0.3530 0.6725
Anosim         0.0535   0.0505 0.2915
Permanova      0.0515   0.0460 0.2975
LRA.ZI         0.0520   0.1560 0.7415
LRA.ZI2        0.0515   0.0480 0.7615


#----------------------------------------
j <- 3
            type1Base type1Cor  power
resLRA         0.0510   0.9905 0.1650
resLogFrac     0.0525   0.1425 0.6365
resWLRA        0.0525   0.7420 0.4560
resCA          0.0495   0.0095 0.1930
resCASqrt      0.0425   0.2035 0.6875
resCALog       0.0405   0.9330 0.3055
Anosim         0.0505   0.1325 0.2950
Permanova      0.0485   0.1145 0.2945
LRA.ZI         0.0520    0.906 0.6060
LRA.ZI2        0.0505    0.237 0.6615




















