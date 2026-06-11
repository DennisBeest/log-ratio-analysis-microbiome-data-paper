

rm(list=ls(all=TRUE))

gitdit <- "E:/1Documenten/Werk/gits/gitlabwur/smallProjects/microbiome/LRApaper/results def"
#gitdit <- "C:/1Files/2Git/smallProjects/microbiome/LRApaper/results def"
setwd(gitdit)

#---load file with functions
source("functions def.R")



###############################################################################


nsim <- 500
#nsim <- 100
#nsim <- 3


sdSiteVec <- c(0.25, 0.5, 1)
sdLength <- length(sdSiteVec)
nScen <- 3



###############################################################################


tryFunc <- function(physeq1) {
    out <- tryCatch(
        {
		suppressWarnings(RCM(physeq1, covariates = "Y", k=1, maxItOut=5))
        },
        error=function(cond) {
            return(NULL)
        },
        warning=function(cond) {
            return(NULL)
        },
        finally={}
    )    
    return(out)
}




outxsAll <- list()
outxsAll[[1]] <- list()
outxsAll[[2]] <- list()
outxsAll[[3]] <- list()


#---sets SD vec
j <- 1


#---each fit takes 1.2 sec
#---A simulation with 200 permutations takes 4 min
#---500 sims take 2000 min or 33 hours. Rougly 15 hours accounting for convergence.
#---200 sims takes 13 hours
#---if 50% of the sims converge this goes down to about 7 hours


set.seed(1122)
for(i in 1:nsim)
{
	#i <- 10
	if(i %% 10 == 0) print(i)

	################################################
	
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


pvalsRCM <- array(dim = c(sdLength, nScen, nsim))
resRCM <- matrix(nrow=sdLength, ncol=nScen)
pvalsRCMConv <- array(dim = c(sdLength, nScen, nsim))
resRCMConv <- matrix(nrow=sdLength, ncol=nScen)


#set the 3 Scenarios
k <- 3



simConv <- logical(nsim)
npermConv <- numeric(nsim) + NA

for(i in 1:nsim)
{	
	print("i")
	print(i)
	#i <- 1

	X <- outxsAll[[k]][[i]]$x
	resp <- outxsAll[[k]][[i]]$response

	rownames(X) <- NULL
	OTU <- otu_table(X, taxa_are_rows = FALSE)

	respDF <- data.frame(Y = factor(resp))
	samp <- sample_data(respDF)
	physeq1 <- phyloseq(OTU, samp)

	simRCM <- tryFunc(physeq1)

#simRCM$runtimeInMins
#simRCM$converged
#simRCM2$iter

	if(!is.null(simRCM)) if(simRCM$converged)
	{
		testVar <- liks(simRCM)[2, 2]
		simConv[i] <- simRCM$converged

		nperm <- 199
		permVar <- numeric(nperm) + NA
		checkConPerm <- logical(nperm)

		for(p in 1:nperm)
		{
			#p <- 148

			#set.seed(p)
			respDF <- data.frame(Y = factor(resp[sample(50)]))
			samp <- sample_data(respDF)
			physeq1 <- phyloseq(OTU, samp)

			simRCMPerm <- tryFunc(physeq1)

			if(is.null(simRCMPerm) == FALSE)
			{
				permVar[p] <- liks(simRCMPerm)[2, 2]
				checkConPerm[p] <- simRCMPerm$converged
			}

			#simRCMPerm <- suppressWarnings(try(RCM(physeq1, covariates = "Y", k=1, maxItOut=5), silent=TRUE))
			#if(class(simRCMPerm) != "try-error")
			#print(simRCMPerm$converged)
			#print(simRCMPerm$iter)
			#liks(simRCMPerm)
		}

		if(any(permVar[checkConPerm] < 0)) print("negative var explained")

		pvalsRCM[j, k, i] <- sum(c(permVar[checkConPerm], testVar) >= testVar) / (sum(checkConPerm) + 1)
		#--note: if not converged variance explained can be negative, non-converged should be removed

		npermConv[i] <- sum(checkConPerm)

	} else {
		pvalsRCM[j, k, i] <- NA
	}
}


pvalsRCM[j, k, ]


paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
print(sum(simConv, na.rm=TRUE))
range(npermConv, na.rm=TRUE)



##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=1 & k=1"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0.018
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.052
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.068
> print(sum(simConv, na.rm=TRUE))
[1] 500
> range(npermConv, na.rm=TRUE)
[1] 199 199
 

##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=2 & k=1"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0.01
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.05
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.064
> print(sum(simConv, na.rm=TRUE))
[1] 500
> range(npermConv, na.rm=TRUE)
[1] 199 199


##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=3 & k=1"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0.002
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.034
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.076
> print(sum(simConv, na.rm=TRUE))
[1] 500
> range(npermConv, na.rm=TRUE)
[1] 199 199


##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=1 & k=2"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0.012
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.04
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.048
> print(sum(simConv, na.rm=TRUE))
[1] 500
> range(npermConv, na.rm=TRUE)
[1] 199 199


##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=2 & k=2"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0.012
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.052
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.066
> print(sum(simConv, na.rm=TRUE))
[1] 500
> range(npermConv, na.rm=TRUE)
[1] 199 199



##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=3 & k=2"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0.014
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.032
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.092
> print(sum(simConv, na.rm=TRUE))
[1] 500
> range(npermConv, na.rm=TRUE)
[1] 199 199


##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=1 & k=3"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0.236
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.544
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.624
> print(sum(simConv, na.rm=TRUE))
[1] 500
> range(npermConv, na.rm=TRUE)
[1] 199 199


##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=2 & k=3"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0.19
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.508
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.612
> print(sum(simConv, na.rm=TRUE))
[1] 500
> range(npermConv, na.rm=TRUE)
[1] 199 199


##################################################


 paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=3 & k=3"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0.102
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.41
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.732
> print(sum(simConv, na.rm=TRUE))
[1] 500
> range(npermConv, na.rm=TRUE)
[1] 199 199




