

rm(list=ls(all=TRUE))

gitdit <- "E:/1Documenten/Werk/gits/gitlabwur/smallProjects/microbiome/LRApaper/results def"
setwd(gitdit)
#---load file with functions
source("functions def.R")


##############################################


Yin <- read.table("extraExamples/yatsunenko2012.txt",h=T,sep="\t",quote="",row.names = 1) 
head(Yin)
dim(Yin)
#downloaded from http://bioinformatics.math.chalmers.se/metabayes/


##############################################



Yin2 <- as.matrix(t(Yin))
subset1 <- colSums(Yin2 > 0) >= 10
subset2 <- colSums(Yin2 > 0) >= 50 & colSums(Yin2) >= 2e3


X1 <- Yin2[, subset1]
X2 <- Yin2[, subset2]
dim(X1)
dim(X2)


#range(rowSums(X1))
#sum(X1==0)/length(X1)
#sum(X2==0)/length(X2)
#1-sum(X2)/sum(X1)
#[1] 0.1907341


##############################################


out1 <- LRA.ade4(X1)
out2 <- LRA.ade4(X2)


##############################################
#----Call the general pictures

p1 <- FigCor(plotInfo1 = out1)

p2 <- FigContr(plotInfo1 = out1, whichContrib=2)
p2 <- p2 + ylab(expression(paste("Log contrib ", 2^nd, " axis")))

p3 <- FigCor(plotInfo1 = out2)

p4 <- FigContr(plotInfo1 = out2, whichContrib=2)
p4 <- p4 + ylab(expression(paste("Log contrib ", 2^nd, " axis")))



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


p1 <- p1 + ggtitle("Human gut")
p1 <- p1 + theme(plot.title = element_text(size = 15, face = "plain", hjust = 0))


p3 <- p3 + ggtitle("Human gut - extra filtering")
p3 <- p3 + theme(plot.title = element_text(size = 15, face = "plain", hjust = 0))



##############################################


dev.new(width=7, height=7)
grid.arrange(
	p1, p2,
	p3, p4,
nrow=2, ncol=2)


doSave <- TRUE
locwd <- paste(gitdit,"/extraExamples", sep="")
nameFile <- "humangut.pdf"
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




