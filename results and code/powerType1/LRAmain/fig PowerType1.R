
rm(list=ls(all=TRUE))


####################################################################################
library(gridExtra)
library(grid)
library(ggplot2)
library(vegan)


gitdit <- "E:/1Documenten/Werk/gits/gitlabwur/projects/3microbiome/LRApaper/results def"
#gitdit <- "C:/1Files/2Git/smallProjects/microbiome/LRApaper/results def"
setwd(gitdit)


load("powerType1/LRAmain/type1/resType1ade4.RData")
ls()
type1 <- t(pvals2)
rm(pvals)
rm(pvals2)




####################################################################################


load("powerType1/LRAmain/power/resPowerade4.RData")
ls()
power <- pvals2[-c(6, 7), ]
respvec <- respvec[-c(6, 7)]
rm(pvals)
rm(pvals2)


####################################################################################



plotData1 <- data.frame(x=rep(respvec, times=dim(power)[2]), y=as.numeric(power), g=factor(rep(1:dim(power)[2], each=dim(power)[1])))
p1 <- ggplot(data=plotData1, aes(x=x, y=y, group=g))
p1 <- p1 + geom_line()
p1 <- p1 + theme_classic()
p1 <- p1 + geom_point(aes(shape=g), size=2.5)
p1 <- p1 + ylab("Fraction rejected hypothesis")
p1 <- p1 + theme(axis.title.y=element_text(size=13, face="plain"))
p1 <- p1 + theme(axis.title.x=element_text(size=13, face="plain"))
p1 <- p1 + theme(legend.position="none")
name1 <- expression("Fold change")
p1 <- p1 + xlab(name1)
p1 <- p1 + ylim(c(0,1))


plotData1 <- data.frame(x=rep(strengthsCorr, times=dim(type1)[2]), y=as.numeric(type1), g=factor(rep(1:dim(type1)[2], each=dim(type1)[1])))
p2 <- ggplot(data=plotData1, aes(x=x, y=y, group=g))
p2 <- p2 + geom_line()
p2 <- p2 + theme_classic()
p2 <- p2 + geom_point(aes(shape=g), size=2.5)
p2 <- p2 + ylab("Fraction rejected hypothesis")
p2 <- p2 + theme(axis.title.y=element_text(size=13, face="plain"))
p2 <- p2 + theme(axis.title.x=element_text(size=17, face="plain"))
p2 <- p2 + theme(legend.position="none")
name2 <- expression(gamma)
p2 <- p2 + xlab(name2)
p2 <- p2 + ylim(c(0,1))


name1 <- expression(paste(sigma[a], " = 0.25"))
name2 <- levels(plotData1$g)[3] <- expression(paste(sigma[a], " = 0.5"))
name3 <- levels(plotData1$g)[4] <- expression(paste(sigma[a], " = 1"))
p1 <- p1 + scale_shape_discrete(name = "", labels = c(name1, name2, name3))
p1 <- p1 + theme(legend.position=c(0.8, 0.24))
p1 <- p1 + theme(legend.key.size = unit(1.15, "cm"))
p1 <- p1 + theme(legend.text = element_text(size = 12))
p1 <- p1 + theme(legend.text = element_text(size = 12))
p1 <- p1 + theme(legend.text.align = 0)
p1 <- p1 + theme(axis.title.y=element_text(size=13, face="plain"))
p1 <- p1 + theme(legend.title=element_blank())
p1 <- p1 + theme(legend.title=element_blank(), legend.margin=margin(c(1, 1, 1, 1)))

yloc <- 0.95
grobA <- grobTree(textGrob("(A)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))
grobB <- grobTree(textGrob("(B)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))
p1 <- p1 + annotation_custom(grobA)
p2 <- p2 + annotation_custom(grobB)



dev.new(width=10, height=5)
grid.arrange(p1, p2, 
	nrow=1, ncol=2)


Save <- FALSE
#Save <- TRUE

if(Save)
{
	pdf("powerType1/LRAmain/PowerType1Main.pdf", width=10, height=5)
	grid.arrange(p1, p2, 
		nrow=1, ncol=2)
	dev.off()

}







