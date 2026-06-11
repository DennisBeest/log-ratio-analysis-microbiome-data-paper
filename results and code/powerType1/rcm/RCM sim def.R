

rm(list=ls(all=TRUE))
gitdit <- "E:/1Documenten/Werk/gits/gitlabwur/smallProjects/microbiome/LRApaper/results def"
#gitdit <- "C:/1Files/2Git/smallProjects/microbiome/LRApaper/results def"
setwd(gitdit)

#---load file with functions
source("functions def.R")


###############################################################################


nsim <- 500


sdSiteVec <- c(0.25, 0.5, 1)
sdLength <- length(sdSiteVec)
nScen <- 3



###############################################################################
#---convergence, j = 1 & k = 1:  46/100
#---convergence, j = 1 & k = 2:  56/100
#---convergence, j = 1 & k = 3:  52/100
#---convergence, j = 2 & k = 1:  52/100
#---convergence, j = 2 & k = 2:  34/100
#---convergence, j = 2 & k = 3:  45/100
#---convergence, j = 3 & k = 2:  5/100
#---convergence, j = 3 & k = 3:  44/100
#---convergence, j = 3 & k = 1:  37/100
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
j <- 3


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


###############################################################################


pvalsRCM <- array(dim = c(sdLength, nScen, nsim))
resRCM <- matrix(nrow=sdLength, ncol=nScen)
pvalsRCMConv <- array(dim = c(sdLength, nScen, nsim))
resRCMConv <- matrix(nrow=sdLength, ncol=nScen)


#the 3 Scenarios
k <- 2

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
[1] 0.004115226
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.004115226
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.02057613
> print(sum(simConv))
[1] 243
> range(npermConv, na.rm=TRUE)
[1]  63 146


##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=2 & k=1"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0.008032129
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.008032129
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.02008032
> print(sum(simConv, na.rm=TRUE))
[1] 249
> range(npermConv, na.rm=TRUE)
[1]  60 142


##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=3 & k=1"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0.00877193
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.02192982
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.04824561
> print(sum(simConv, na.rm=TRUE))
[1] 228
> range(npermConv, na.rm=TRUE)
[1]  60 153

##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=1 & k=2"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0.007692308
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.007692308
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.03076923
> print(sum(simConv, na.rm=TRUE))
[1] 260
> range(npermConv, na.rm=TRUE)
[1]  57 148


##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=2 & k=2"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.01863354
> print(sum(simConv, na.rm=TRUE))
[1] 161
> range(npermConv, na.rm=TRUE)
[1]  60 162

##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=3 & k=2"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.05
> print(sum(simConv, na.rm=TRUE))
[1] 20
> range(npermConv, na.rm=TRUE)
[1]  83 129


##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=1 & k=3"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0.3417722
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.3586498
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.4303797
> print(sum(simConv))
[1] 237
> range(npermConv, na.rm=TRUE)
[1]  67 147


##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=2 & k=3"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0.3914894
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.4255319
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.4851064
> print(sum(simConv))
[1] 235
> range(npermConv, na.rm=TRUE)
[1]  67 145


##################################################


> paste(paste("j=",j, sep=""), paste("k=",k, sep=""), sep=" & ")
[1] "j=3 & k=3"
> mean(pvalsRCM[j, k, ] < 0.05, na.rm=TRUE)
[1] 0.27897
> mean(pvalsRCM[j, k, ] < 0.10, na.rm=TRUE)
[1] 0.3304721
> mean(pvalsRCM[j, k, ] < 0.20, na.rm=TRUE)
[1] 0.4506438
> print(sum(simConv, na.rm=TRUE))
[1] 233
> range(npermConv, na.rm=TRUE)
[1]  64 142



