
rm(list=ls(all=TRUE))
gitdit <- "E:/1Documenten/Werk/gits/gitlabwur/projects/3microbiome/LRApaper/results def"
setwd(gitdit)
#---load file with functions
source("functions def.R")


respvec <- c(1, 1.25, 1.5, 1.75, 2, 2.25, 2.5)
nresp <- length(respvec)
nsim <- 2000
#nsim <- 4
sdSite <- c(0.25, 0.5, 1)
sdLength <- length(sdSite)
pvals=nzeros <- array(dim = c(nresp, sdLength, nsim))
pvals2=nzeros2 <- matrix(nrow=nresp, ncol=sdLength)


set.seed(2000)
for(j in 1:nresp)
{
	#j <- 1
	print("j")
	print(j)

	for(k in 1:sdLength)
	{
		#k <- 2

		for(i in 1:nsim)
		{
			#i <- 1
			outx <- dataSim(sdSite=sdSite[k], sdSpecies=2, weightResp=respvec[j], drawDistr="nb", Max=1, strengthCorrResp=0)

			x <- outx$xs
			Xtr <- DoubleCenter(log(x + 1))
			pca1 <- dudi.pca(as.data.frame(Xtr), scannf = FALSE, scale = FALSE, center = FALSE)
			rda1 <- pcaiv(dudi=pca1, df=outx$resp, scannf = FALSE)
			rt <- randtest(rda1, nrepet = 999)
			pvals[j, k, i] <- rt$pvalue
			nzeros[j, k, i] <- sum(x==0)/length(x)

		}

		pvals2[j, k] <- mean(pvals[j, k, ] < 0.05)
		nzeros2[j, k] <- mean(nzeros[j, k, ])
	}
}



rownames(pvals2) <- paste("resp", respvec)
colnames(pvals2) <- paste("sd", sdSite)
pvals2


#> pvals2
#          sd 0.25 sd 0.5   sd 1
#resp 1     0.0520 0.0500 0.0555
#resp 1.25  0.2000 0.1555 0.0675
#resp 1.5   0.7935 0.6840 0.1705
#resp 1.75  0.9960 0.9815 0.5245
#resp 2     1.0000 1.0000 0.8670
#resp 2.25  1.0000 1.0000 0.9840
#resp 2.5   1.0000 1.0000 1.0000


plot(respvec, pvals2[,1], cex=0.1)
lines(respvec, pvals2[,1], col="red")

points(respvec, pvals2[,2], pch=20)
lines(respvec, pvals2[,2], pch=20)

points(respvec, pvals2[,3], col="red")
lines(respvec, pvals2[,3], col="red")


save(respvec, pvals, pvals2, file="powerType1/LRAmain/power/resPowerade4.RData")



