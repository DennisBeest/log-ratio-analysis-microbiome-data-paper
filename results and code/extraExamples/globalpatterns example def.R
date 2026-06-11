
rm(list=ls(all=TRUE))
library(phyloseq)

#gitdit <- "E:/1Documenten/Werk/gits/gitlabwur/smallProjects/microbiome/LRApaper/results def"
gitdit <- "/home/dennis/Documents/gitwur/smallProjects/microbiome/LRApaper/results def"
setwd(gitdit)
#---load file with functions
source("functions def.R")


##############################################



data(GlobalPatterns)
Yin <- t(otu_table(GlobalPatterns))
sample_data(GlobalPatterns)


Yin2 <- as.matrix(Yin)
subset1 <- colSums(Yin2 > 0) >= 5
subset2 <- colSums(Yin2 > 0) >= 20 & colSums(Yin2) >= 1000
X1 <- Yin2[, subset1]
X2 <- Yin2[, subset2]
dim(X1)
dim(X2)


#range(rowSums(X1))
#sum(X1==0)/length(X1)
#sum(X2==0)/length(X2)
#1-sum(X2)/sum(X1)
#[1] 0.2104025


##############################################


out1 <- LRA.ade4(X1)
out2 <- LRA.ade4(X2)


##############################################


p1 <- FigCor(plotInfo1 = out1)


p2 <- FigContr(plotInfo1 = out1, whichContrib=1)


p3 <- FigCor(plotInfo1 = out2)


p4 <- FigContr(plotInfo1 = out2, whichContrib=1)



##############################################



yloc <- 0.95
grobA <- grobTree(textGrob("(A)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))
grobB <- grobTree(textGrob("(B)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))
grobC <- grobTree(textGrob("(C)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))
grobD <- grobTree(textGrob("(D)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))
p1 <- p1 + annotation_custom(grobA)
p2 <- p2 + annotation_custom(grobB)
p3 <- p3 + annotation_custom(grobC)
p4 <- p4 + annotation_custom(grobD)


p1 <- p1 + ggtitle("Global patterns")
p1 <- p1 + theme(plot.title = element_text(size = 15, face = "plain", hjust = 0))


p3 <- p3 + ggtitle("Global patterns-extra filtering")
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

nameFile <- "globalpatterns.pdf"
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



