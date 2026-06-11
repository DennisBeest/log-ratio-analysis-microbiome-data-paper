


rm(list=ls(all=TRUE))
gitdit <- "E:/1Documenten/Werk/gits/gitlabwur/projects/3microbiome/LRApaper/results def"
setwd(gitdit)
#---load file with functions
source("functions def.R")


respvec <- c(1)
nresp <- length(respvec)
nsim <- 2000
sdSiteVec <- c(0.25, 0.5, 1)
sdLength <- length(sdSiteVec)
strengthsCorr <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)
strenghLength <- length(strengthsCorr)


pvals <- array(dim = c(sdLength, strenghLength, nsim))
pvals2 <- matrix(nrow=sdLength, ncol=strenghLength)


set.seed(2000)
for(j in 1:sdLength)
{
	#j <- 1
	print(j)
	print("j")

	for(k in 1:strenghLength)
	{
		#k <- 1
		#[1] 1660

		for(i in 1:nsim)
		{
			#i <- 39

			outx <- dataSim(sdSite=sdSiteVec[j], sdSpecies=2, weightResp=respvec, drawDistr="nb", Max=1, strengthCorrResp=strengthsCorr[k])

			x <- outx$xs
			Xtr <- DoubleCenter(log(x + 1))
			pca1 <- dudi.pca(as.data.frame(Xtr), scannf = FALSE, scale = FALSE, center = FALSE)
			rda1 <- pcaiv(dudi=pca1, df=outx$resp, scannf = FALSE)
			rt <- randtest(rda1, nrepet = 999)
			pvals[j, k, i] <- rt$pvalue
		}

		pvals2[j, k] <- mean(pvals[j, k, ] < 0.05)
	}
}


colnames(pvals2) <- paste("strengthsCorr", strengthsCorr)
rownames(pvals2) <- paste("sdSiteVec", sdSiteVec)
t(pvals2)
 

#resType1ade4
#                  sdSiteVec 0.25 sdSiteVec 0.5 sdSiteVec 1
#strengthsCorr 0           0.0520        0.0515      0.0435
#strengthsCorr 0.5         0.0585        0.0945      0.3840
#strengthsCorr 1           0.0595        0.2560      0.8415
#strengthsCorr 1.5         0.0710        0.4780      0.9695
#strengthsCorr 2           0.0790        0.6685      0.9920
#strengthsCorr 2.5         0.1180        0.8175      0.9990
#strengthsCorr 3           0.1560        0.8795      0.9985

save(strengthsCorr, pvals, pvals2, file="powerType1/LRAmain/type1/resType1ade4.RData")


plot(strengthsCorr, pvals2[3,])
lines(strengthsCorr, pvals2[3, ])
points(strengthsCorr, pvals2[2,], col="red")
lines(strengthsCorr, pvals2[2,], col="red")
points(strengthsCorr, pvals2[1,], col="blue")
lines(strengthsCorr, pvals2[1,], col="blue")

