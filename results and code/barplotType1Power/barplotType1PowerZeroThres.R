

rm(list=ls(all=TRUE))

library(ggplot2)
library(gridExtra)
library(grid)

gitdit <- "E:/1Documenten/Werk/gits/gitlabwur/projects/3microbiome/LRApaper/results def"
#gitdit <- "C:/1Files/2Git/smallProjects/microbiome/LRApaper/results def"
setwd(gitdit)


data <- read.table("barplotType1Power/resSuppl.csv", sep = ",", header = TRUE, stringsAsFactor = FALSE)
head(data)


#########################################################################################################################


data$Method <- replace(data$Method, data$Method == "Log-ratio RDA GBM/counts", "Log-ratio RDA GBM/cnts")
data$Method <- replace(data$Method, data$Method == "Anosim", "Anosim (BC-prop)")
data$Method <- replace(data$Method, data$Method == "Permanova", "Permanova (BC-prop)")


data$Method <- factor(data$Method, levels = c("Log-ratio RDA", 
	"Weighted log-ratio RDA",
	"Log-ratio RDA GBM/cnts",
	"Log-ratio RDA GBM/prop",
	"Log proportions RDA",
	"CCA - counts", 
	"CCA - sqrt", 
	"CCA - log",
	"Anosim (BC-prop)",
	"Permanova (BC-prop)",
	"RCM"))

#cbind(levels(data$Method))


fPlot <- function(x) 
{
	p <- ggplot(data=dfIn, aes(x=Method, y=probs))
	p <- p + geom_bar(stat="identity", width=0.95, fill="darkgrey")
	p <- p + theme_classic()
	p <- p + xlab("")
	p <- p + ylim(c(0, 1))
	p <- p + theme(axis.text.y=element_text(size=12))
	#p <- p + theme(axis.title.y=element_text(size = 14))
	p <- p + ylab("")
	p <- p + theme(axis.ticks = element_blank())
	return(p)
}


ps <- list()
for(i in 1:6)
{
	#i <- 1

	if(i == 1) index <- data$sd == 0.25 & data$Filt == "Filt"
	if(i == 2) index <- data$sd == 1 & data$Filt == "Filt"

	if(i == 3) index <- data$sd == 0.25 & data$Filt == "Filt"
	if(i == 4) index <- data$sd == 1 & data$Filt == "Filt"

	if(i == 5) index <- data$sd == 0.25 & data$Filt == "Filt"
	if(i == 6) index <- data$sd == 1 & data$Filt == "Filt"

	if(i == 1 | i == 2)
	{
		dfIn <- data[index, c("Method", "type1Base")]
		colnames(dfIn)[2] <- "probs"
	}
	if(i == 3 | i == 4)
	{
		dfIn <- data[index, c("Method", "power")]
		colnames(dfIn)[2] <- "probs"
	}
	if(i == 5 | i == 6)
	{
		dfIn <- data[index, c("Method", "type1Cor")]
		colnames(dfIn)[2] <- "probs"
	}

	ps[[i]] <- fPlot(dfIn)

	if(!(i == 3 | i == 4)) ps[[i]] <- ps[[i]] + geom_abline(slope=0, intercept=0.05,  lty=1, lwd = 1.5)
}



ps[[5]] <- ps[[5]] + theme(axis.text.x=element_text(size = 12, face = "bold"))
ps[[5]] <- ps[[5]] + theme(axis.text.x = element_text(angle = 90))
ps[[6]] <- ps[[6]] + theme(axis.text.x=element_text(size = 12, face = "bold"))
ps[[6]] <- ps[[6]] + theme(axis.text.x = element_text(angle = 90))

for(i in 1:4) ps[[i]] <- ps[[i]] + theme(axis.text.x=element_blank())


#ps[[1]] <- ps[[1]] + ylim(c(0, 1.05))
#ps[[2]] <- ps[[2]] + ylim(c(0, 1.05))


t2 <- textGrob("Fraction rejected hypothesis", 
               gp = gpar(fontsize = 31), 
               rot = 90, vjust = 0.7, hjust = 0.28)



#top, right, bot, left
ps[[1]] <- ps[[1]] + theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm"))
ps[[2]] <- ps[[2]] + theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm"))
for(i in 3:6) ps[[i]] <- ps[[i]] + theme(plot.margin = unit(c(-0.2, 0, 0, 0), "cm"))



fontSize <- 19
xloc1 <- 0.12
xloc2 <- 0.01
text1 <- expression(paste("(A) ","Type 1 ; ", sigma[a], " = 0.25"," ; ", gamma, " = 0"))
grob1 <- grobTree(textGrob(text1, x=xloc2,  y=0.96, gp=gpar(fontsize = fontSize), hjust = 0))

text2 <- expression(paste("(B) ","Type 1 ; ", sigma[a], " = 1"," ; ", gamma, " = 0", ""))
grob2 <- grobTree(textGrob(text2, x=xloc2,  y=0.96, gp=gpar(fontsize = fontSize), hjust = 0))

text3 <- expression(paste("(C) ","Power ; ", sigma[a], " = 0.25"," ; ", gamma, " = 0", " ; FC = 1.5"))
grob3 <- grobTree(textGrob(text3, x=xloc2,  y=0.96, gp=gpar(fontsize = fontSize), hjust = 0))

text4 <- expression(paste("(D) ","Power ; ", sigma[a], " = 1"," ; ", gamma, " = 0", " ; FC = 1.5"))
grob4 <- grobTree(textGrob(text4, x=xloc2,  y=0.96, gp=gpar(fontsize = fontSize), hjust = 0))

text5 <- expression(paste("(E) ","Type 1 ; ", sigma[a], " = 0.25"," ; ", gamma, " = 2"))
grob5 <- grobTree(textGrob(text5, x=xloc2,  y=0.96, gp=gpar(fontsize = fontSize), hjust = 0))

text6 <- expression(paste("(F) ","Type 1 ; ", sigma[a], " = 1"," ; ", gamma, " = 2"))
grob6 <- grobTree(textGrob(text6, x=xloc2,  y=0.96, gp=gpar(fontsize = fontSize), hjust = 0))



ps[[1]] <- ps[[1]] + annotation_custom(grob1)
ps[[2]] <- ps[[2]] + annotation_custom(grob2)
ps[[3]] <- ps[[3]] + annotation_custom(grob3)
ps[[4]] <- ps[[4]] + annotation_custom(grob4)
ps[[5]] <- ps[[5]] + annotation_custom(grob5)
ps[[6]] <- ps[[6]] + annotation_custom(grob6)



dev.new(width=16, height=12)
grid.arrange(
	ps[[1]], ps[[2]],
	ps[[3]], ps[[4]],
	ps[[5]], ps[[6]],
	heights = c(1, 1, 1.7),
	left = t2,
nrow=3, ncol=2)




##############################################################################




nameFile <- "barplotType1Power/barplotType1PowerZeroThres.pdf"
savefile <- TRUE
if(savefile)
{

	pdf(nameFile, width=16, height=12)
		grid.arrange(
			ps[[1]], ps[[2]],
			ps[[3]], ps[[4]],
			ps[[5]], ps[[6]],
			heights = c(1, 1, 1.7),
			left = t2,
		nrow=3, ncol=2)

	dev.off()

}


