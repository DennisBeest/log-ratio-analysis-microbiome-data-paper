


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


corrs1 <- array(dim = c(sdLength, strenghLength, nsim))
corrs2 <- array(dim = c(sdLength, strenghLength, nsim))


#####################################################################


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
			#i <- 1

			outx <- dataSim(sdSite=sdSiteVec[j], sdSpecies=2, weightResp=respvec, drawDistr="nb", Max=1, strengthCorrResp=strengthsCorr[k])
			Xtr <- DoubleCenter(log(outx$xs + 1))

			corrs1[j, k, i] <- cor(outx$siteEffect, outx$resp)
			corrs2[j, k, i] <- cor(attr(Xtr, "rowmeans"), outx$resp)
		}
	}
}


#####################################################################



means1 <- matrix(nrow = sdLength, ncol = strenghLength)
means2 <- matrix(nrow = sdLength, ncol = strenghLength)

for(j in 1:sdLength)
{
	#j <- 1
	print(j)
	print("j")

	for(k in 1:strenghLength) 
	{
		means1[j, k] <- mean(corrs1[j, k, ])
		means2[j, k] <- mean(corrs2[j, k, ])
	}
}



#save(means, file = "picCorGamma/meansCorrAX.RData")

attr(means1, "what") <- "mean cor(outx$siteEffect, outx$resp)"
attr(means2, "what") <- "mean cor(attr(Xtr, rowmeans), outx$resp"
save(means1, means2, file = "picCorGamma/meansCorrRhoAX.RData")




#####################################################################



#load("picCorGamma/meansCorrAX.RData")
load("picCorGamma/meansCorrRhoAX.RData")


pdf("picCorGamma/picCorGamma.pdf", width=8, height=8)


par(mfrow = c(2, 1))
par(oma = c(1, 1, 0, 0))
par(mar = c(3.3, 3.3, 1, 1))
par(mgp = c(2.3, 1, 0))


ii <- 1
plot(strengthsCorr, means1[ii, ], pch = 20, xlim = range(strengthsCorr) * 1.02, ylim = c(0, 0.7), xlab = expression(gamma), ylab = expression(paste(rho[ax])), cex.lab = 1.2)
lines(strengthsCorr, means1[ii, ])
ii <- 2
points(strengthsCorr, means1[ii, ], pch = 20)
lines(strengthsCorr, means1[ii, ], lty = 2)
ii <- 3
points(strengthsCorr, means1[ii, ], pch = 20)
lines(strengthsCorr, means1[ii, ], lty = 3)

shift <- 0.05
legend(0.01, 0.63 + shift,  expression(paste(sigma[a], " = 1")), lty = 3, lwd = 2, bty = "n", cex = 1.1)
legend(0.01, 0.57 + shift, expression(paste(sigma[a], " = 0.5")), lty = 2, lwd = 2, bty = "n", cex = 1.1)
legend(0.01, 0.51 + shift, expression(paste(sigma[a], " = 0.25")), lty = 1, lwd = 2, bty = "n", cex = 1.1)


ii <- 1
plot(strengthsCorr, means2[ii, ], pch = 20, xlim = range(strengthsCorr) * 1.02, ylim = c(0, 0.7), xlab = expression(gamma), ylab = expression(paste(rho[rx])), cex.lab = 1.2)
lines(strengthsCorr, means2[ii, ])
ii <- 2
points(strengthsCorr, means2[ii, ], pch = 20)
lines(strengthsCorr, means2[ii, ], lty = 2)
ii <- 3
points(strengthsCorr, means2[ii, ], pch = 20)
lines(strengthsCorr, means2[ii, ], lty = 3)



dev.off()


