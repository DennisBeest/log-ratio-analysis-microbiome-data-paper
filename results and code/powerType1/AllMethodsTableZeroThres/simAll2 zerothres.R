


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


#--this sets the sd. Repeat this for j = 1, 2, 3.
j <- 1



###############################################################################


outxsAll <- list()
outxsAll[[1]] <- list()
outxsAll[[2]] <- list()
outxsAll[[3]] <- list()


set.seed(1122)
for(i in 1:nsim)
{
	#i <- 1
	if(i %% 10 == 0) print(i)

	outxsAll[[1]][[i]] <- dataSim(sdSite=sdSiteVec[j], sdSpecies=2, weightResp=1, drawDistr="nb", Max=1, strengthCorrResp=0, zeroThres=5)
	outxsAll[[2]][[i]] <- dataSim(sdSite=sdSiteVec[j], sdSpecies=2, weightResp=1, drawDistr="nb", Max=1, strengthCorrResp=2, zeroThres=5)
	outxsAll[[3]][[i]] <- dataSim(sdSite=sdSiteVec[j], sdSpecies=2, weightResp=1.5, drawDistr="nb", Max=1, strengthCorrResp=0, zeroThres=5)
}


for(i in 1:nsim)
{
	ss <- colSums(outxsAll[[1]][[i]]$xs!=0) >= 30 & colSums(outxsAll[[1]][[i]]$xs) >= 500
	outxsAll[[1]][[i]]$xs <- outxsAll[[1]][[i]]$xs[, ss]
	ss <- colSums(outxsAll[[2]][[i]]$xs!=0) >= 30 & colSums(outxsAll[[2]][[i]]$xs) >= 500
	outxsAll[[2]][[i]]$xs <- outxsAll[[2]][[i]]$xs[, ss]
	ss <- colSums(outxsAll[[3]][[i]]$xs!=0) >= 30 & colSums(outxsAll[[3]][[i]]$xs) >= 500
	outxsAll[[3]][[i]]$xs <- outxsAll[[3]][[i]]$xs[, ss]
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
	if(i %% 2 == 0) print(i)

	for(k in 1:3)
	{
		#k <- 2

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



########################################################################################################################
#----26 sept
########################################################################################################################



j <- 1
            type1Base type1Cor  power
resLRA         0.0450   0.0455 0.4010
resLogFrac     0.0480   0.0495 0.4050
resWLRA        0.0440   0.0465 0.3610
resCA          0.0445   0.0500 0.2680
resCASqrt      0.0410   0.0410 0.4130
resCALog       0.0365   0.0425 0.3675
Anosim         0.0480   0.0490 0.2380
Permanova      0.0470   0.0545 0.2425
LRA.ZI         0.0460   0.0500 0.3880
LRA.ZI2        0.0460   0.0500 0.3845


j <- 2
            type1Base type1Cor  power
resLRA         0.0480   0.0460 0.4185
resLogFrac     0.0530   0.0425 0.4095
resWLRA        0.0570   0.0390 0.3595
resCA          0.0640   0.0270 0.2420
resCASqrt      0.0450   0.0280 0.4100
resCALog       0.0340   0.0940 0.3705
Anosim         0.0515   0.0445 0.2375
Permanova      0.0505   0.0420 0.2410
LRA.ZI         0.0430   0.0420 0.3710
LRA.ZI2        0.0465   0.0450 0.3760


j <- 3
            type1Base type1Cor  power
resLRA         0.0535   0.0925 0.4740
resLogFrac     0.0550   0.0930 0.4650
resWLRA        0.0525   0.0240 0.3785
resCA          0.0515   0.0105 0.1720
resCASqrt      0.0435   0.0125 0.3965
resCALog       0.0430   0.7655 0.3820
Anosim         0.0535   0.0565 0.2430
Permanova      0.0495   0.0540 0.2550
LRA.ZI         0.0520   0.0755 0.4100
LRA.ZI2        0.0510   0.0775 0.4025










