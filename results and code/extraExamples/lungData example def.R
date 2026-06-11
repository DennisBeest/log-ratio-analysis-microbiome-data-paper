

rm(list=ls(all=TRUE))


gitdit <- "E:/1Documenten/Werk/gits/gitlabwur/smallProjects/microbiome/LRApaper/results def"
setwd(gitdit)

#---load file with functions
source("functions def.R")


library(vegan)
library(ade4)
library(ggplot2)
library(gridExtra)
library(grid)
library(vegan)


################################################################################################
#---Load data
################################################################################################



library(metagenomeSeq)
data(lungData)
Yin <- t(MRcounts(lungData))
Yin <- as.matrix(Yin)
dim(Yin)


smokingStatus <- pData(lungData)$SmokingStatus
subsetNA <- !is.na(smokingStatus)
Yin2 <- Yin[subsetNA, ]
response <- smokingStatus[subsetNA]


subset1 <- colSums(Yin2 > 0) >= 10
subset2 <- colSums(Yin2 > 0) >= 20 & colSums(Yin2) >= 250

X1 <- Yin2[, subset1]
X2 <- Yin2[, subset2]


#dim(X1)
#dim(X2)
#sum(Yin2[, subset1])
#sum(Yin2[, subset2])
#1-sum(X2)/sum(X1)
#[1] 0.523432
#range(rowSums(X1))
#sum(X1==0)/length(X1)
#sum(X2==0)/length(X2)


################################################################################################


out1 <- LRA.ade4(X1)
out2 <- LRA.ade4(X2)


################################################################################################
#----Call the general pictures


p1 <- FigCor(out1)

p2 <- FigContr(plotInfo1 = out1, whichContrib=1)

p3 <- FigCor(plotInfo1 = out2)

p4 <- FigContr(plotInfo1 = out2, whichContrib=1)



##############################################
#----Some modifying of the general pictures


yloc <- 0.95
grobA <- grobTree(textGrob("(A)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))
grobB <- grobTree(textGrob("(B)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))
grobC <- grobTree(textGrob("(C)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))
grobD <- grobTree(textGrob("(D)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))
p1 <- p1 + annotation_custom(grobA)
p2 <- p2 + annotation_custom(grobB)
p3 <- p3 + annotation_custom(grobC)
p4 <- p4 + annotation_custom(grobD)


p1 <- p1 + ggtitle("Lung")
p1 <- p1 + theme(plot.title = element_text(size = 15, face = "plain", hjust = 0))


p3 <- p3 + ggtitle("Lung - extra filtering")
p3 <- p3 + theme(plot.title = element_text(size = 15, face = "plain", hjust = 0))



##############################################



dev.new(width=7, height=7)
grid.arrange(
	p1, p2,
	p3, p4,
	widths = c(3.3, 3.3),
nrow=2, ncol=2)



doSave <- TRUE
locwd <- paste(gitdit,"/extraExamples", sep="")
nameFile <- "lungdata.pdf"
if(doSave)
{
	setwd(locwd)
	pdf(nameFile,width=7, height=7)
		grid.arrange(
			p1, p2,
			p3, p4,
		nrow=2, ncol=2)
	dev.off() 
}




