
library(vegan)
library(ggplot2)
library(gridExtra)
library(grid)
library(ade4)
library(MCMCpack)
library(vegan)
library(zCompositions)
library(RCM)
library(phyloseq)



###################################################################################
#---Generic picture to plot contributions
###################################################################################



FigContr <- function(plotInfo1, whichContrib=1, contrained=FALSE)
{
	#plotInfo1 <- lraRice

	yxLabels <- 13
	yxTicks <- 10

	if(!contrained) df1 <- data.frame(x=plotInfo1$logColmeans,y=plotInfo1$Contr[, whichContrib])
	if(contrained) df1 <- data.frame(x=plotInfo1$logColmeans,y=plotInfo1$ContrConstr)

	p <- ggplot(df1, aes(x, y))

	#--plot points
	p <- p + geom_point()

	#--removes background
	p <- p + theme_classic()

	#--name xlab
	p <- p + xlab("log of taxon mean")

	#--naam ylab
	if(!contrained & whichContrib==1) p <- p + ylab(expression(paste("Log contrib ", 1^st, " axis")))
	if(!contrained & whichContrib==2) p <- p + ylab(expression(paste("Log contrib ", 2^nd, " axis")))
	if(contrained) p <- p + ylab("Log contribution r-contrained axis") 

	#--options for ticks
	p <- p + theme(axis.text.x = element_text(angle = 0, hjust = 1, size=yxTicks, color="black"))
	#--options voor text of labels
	p <- p + theme(axis.text.y = element_text(angle = 0, hjust = 1, size=yxTicks, color="black")) 

	#--size x label
	p <- p + theme(axis.title.x=element_text(size=yxLabels, face="plain"))
	#--size y label
	p <- p + theme(axis.title.y=element_text(size=yxLabels, face="plain"))

	#--add title
	p <- p + ggtitle("")
	#--size title
	p <- p + theme(title=element_text(size=14,face="plain"))
 
	return(p)
}


########################################################################################
#----Generic picture with ordination
########################################################################################


FigOrd <- function(plotInfo1)
{
	yxLabels <- 13

	df3 <- data.frame(x=plotInfo1$Sc[,1], y=plotInfo1$Sc[,2], group=factor(plotInfo1$resp))

	p <- ggplot(df3, aes(x, y))

	#--remove background
	p <- p + theme_classic()

	#--plot points with color
	p <- p + geom_point(aes(colour=group))
	#--set color of groups
	p <- p + scale_colour_manual(values=c("black", "red"))

	#--plot points with shape
	#p <- p + geom_point(aes(shape=group))
	#--set shape
	#p <- p + scale_shape_manual(values=c(20, 4))

	#--name xlab
	p <- p + xlab(expression(paste(1^st, " axis")))
	#--name ylab
	p <- p + ylab(expression(paste(2^nd, " axis")))

	#--size x label
	p <- p + theme(axis.title.x=element_text(size=yxLabels, face="plain"))
	#--size y label
	p <- p + theme(axis.title.y=element_text(size=yxLabels, face="plain"))

	#--remove x axis labels/ticks
	p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
	#--remove y axis labels/ticks
	p <- p + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

	#--remove legend
	p <- p + theme(legend.position = "none")

	return(p)
}


####################################################################
#---Generic figure with correlation
####################################################################


FigCor <- function(plotInfo1, plotBoth=FALSE)
{
	yxLabels <- 13
	yxTicks <- 10

	df2 <- data.frame(x=plotInfo1$logColmeans, y=plotInfo1$Corrs)

	#--start
	p <- ggplot(df2, aes(x, y))
	#--plot points
	p <- p + geom_point()
	#--remove background
	p <- p + theme_classic()

	#--name xlab
	p <- p + xlab("log of taxon mean")
	#--options for text of ticks
	p <- p + theme(axis.text.x = element_text(angle = 0, hjust = 1, size=yxTicks,color="black"))
	#--size x label
	p <- p + theme(axis.title.x=element_text(size = yxLabels, face= "plain"))

	#--name ylab
	corExpr <- expression(paste("Correlation (", rho[Sr], ")"))
	p <- p + ylab(corExpr) 
	#--options for text of ticks
	p <- p + theme(axis.text.y = element_text(angle = 0, hjust = 1, size=yxTicks,color="black")) 
	#--size y label
	p <- p + theme(axis.title.y=element_text(size = yxLabels, face = "plain"))

	#--no title
	p <- p + ggtitle("")

	#--set limits y-axis
	p <- p + ylim(c(-1, 1))

	return(p)
}


####################################################################
#----picture supplement
####################################################################


bothPics1 <- function(lraRice, lraKnut, lraSim, saveInfo, whichContrib = c(1, 1, 1), OrientOrdi = c(1, 1, 1, 1, 1, 1), textAnalysis=NULL, hjustTitle2=1)
{

	#--adjust orientation Ordination diagram
	lraRice$Sc[, 1] <- lraRice$Sc[, 1] * OrientOrdi[1]
	lraRice$Sc[, 2] <- lraRice$Sc[, 2] * OrientOrdi[2]

	lraKnut$Sc[, 1] <- lraKnut$Sc[, 1] * OrientOrdi[3]
	lraKnut$Sc[, 2] <- lraKnut$Sc[, 2] * OrientOrdi[4]

	lraSim$Sc[, 1] <- lraSim$Sc[, 1] * OrientOrdi[5]
	lraSim$Sc[, 2] <- lraSim$Sc[, 2] * OrientOrdi[6]


	ps1 <- list()
	ps1[[1]] <- FigOrd(lraSim)
	ps1[[2]] <- FigContr(lraSim, whichContrib=whichContrib[1])

	ps2 <- list()
	ps2[[1]] <- FigOrd(lraRice)
	ps2[[2]] <- FigContr(lraRice, whichContrib=whichContrib[2])

	ps3 <- list()
	ps3[[1]] <- FigOrd(lraKnut)
	ps3[[2]] <- FigContr(lraKnut, whichContrib=whichContrib[3])


	#########################################################################


	yloc <- 0.96
	xloc <- 0.1
	fontsiz <- 9
	grobA <- grobTree(textGrob("(A)", x=xloc,  y=yloc, gp=gpar(fontsize=fontsiz)))
	grobB <- grobTree(textGrob("(B)", x=xloc,  y=yloc, gp=gpar(fontsize=fontsiz)))
	grobC <- grobTree(textGrob("(C)", x=xloc,  y=yloc, gp=gpar(fontsize=fontsiz)))
	grobD <- grobTree(textGrob("(D)", x=xloc,  y=yloc, gp=gpar(fontsize=fontsiz)))
	grobE <- grobTree(textGrob("(E)", x=xloc,  y=yloc, gp=gpar(fontsize=fontsiz)))
	grobF <- grobTree(textGrob("(F)", x=xloc,  y=yloc, gp=gpar(fontsize=fontsiz)))


	ps1[[1]] <- ps1[[1]] + annotation_custom(grobA)
	ps1[[2]] <- ps1[[2]] + annotation_custom(grobB)

	ps3[[1]] <- ps3[[1]] + annotation_custom(grobC)
	ps3[[2]] <- ps3[[2]] + annotation_custom(grobD)

	ps2[[1]] <- ps2[[1]] + annotation_custom(grobE)
	ps2[[2]] <- ps2[[2]] + annotation_custom(grobF)


	#########################################################################


	ps1[[1]] <- ps1[[1]] + ggtitle("Simulated example")
	ps1[[1]] <- ps1[[1]] + theme(plot.title = element_text(size = 16, face = "plain", hjust = 0))

	ps2[[1]] <- ps2[[1]] + ggtitle("Rice example")
	ps2[[1]] <- ps2[[1]] + theme(plot.title = element_text(size = 16, face = "plain", hjust = 0))

	ps3[[1]] <- ps3[[1]] + ggtitle("Biting midges example")
	ps3[[1]] <- ps3[[1]] + theme(plot.title = element_text(size = 16, face = "plain", hjust = 0))

	ps1[[2]] <- ps1[[2]] + ggtitle(textAnalysis)
	ps1[[2]] <- ps1[[2]] + theme(plot.title = element_text(size = 15, face = "plain", hjust = hjustTitle2))


	##########################################################################


	doSave <- saveInfo[[1]]
	locwd <- saveInfo[[2]]
	nameFile <- saveInfo[[3]]

	if(doSave)
	{
		setwd(locwd)
		pdf(nameFile,width=6, height=8)
			grid.arrange(
			ps1[[1]], ps1[[2]],
			ps3[[1]], ps3[[2]],
			ps2[[1]], ps2[[2]],
			widths = c(3.3, 3.7),
		nrow=3, ncol=2)
		dev.off() 
	}

	dev.new(width=6, height=8)
	grid.arrange(
		ps1[[1]], ps1[[2]],
		ps3[[1]], ps3[[2]],
		ps2[[1]], ps2[[2]],
		widths = c(3.3, 3.7),
	nrow=3, ncol=2)
}


####################################################################
#----picture supplement
####################################################################


bothPics2 <- function(lraRice, lraKnut, lraSim, saveInfo, whichContrib=c(1, 1, 1), OrientOrdi = c(1, 1, 1, 1, 1, 1), textAnalysis=NULL, hjustTitle2=0.5)
{

	lraRice$Sc[, 1] <- lraRice$Sc[, 1] * OrientOrdi[1]
	lraRice$Sc[, 2] <- lraRice$Sc[, 2] * OrientOrdi[2]
	lraKnut$Sc[, 1] <- lraKnut$Sc[, 1] * OrientOrdi[3]
	lraKnut$Sc[, 2] <- lraKnut$Sc[, 2] * OrientOrdi[4]
	lraSim$Sc[, 1] <- lraSim$Sc[, 1] * OrientOrdi[5]
	lraSim$Sc[, 2] <- lraSim$Sc[, 2] * OrientOrdi[6]

	ps1 <- list()
	ps1[[1]] <- FigOrd(lraSim)
	ps1[[2]] <- FigCor(lraSim)
	ps1[[3]] <- FigContr(lraSim, whichContrib=whichContrib[1])

	ps2 <- list()
	ps2[[1]] <- FigOrd(lraRice)
	ps2[[2]] <- FigCor(lraRice)
	ps2[[3]] <- FigContr(lraRice, whichContrib=whichContrib[2])

	ps3 <- list()
	ps3[[1]] <- FigOrd(lraKnut)
	ps3[[2]] <- FigCor(lraKnut)
	ps3[[3]] <- FigContr(lraKnut, whichContrib=whichContrib[3])


	#########################################################################


	yloc <- 0.96
	xloc <- 0.1
	fontsiz <- 9
	grobA <- grobTree(textGrob("(A)", x=xloc,  y=yloc, gp=gpar(fontsize=fontsiz)))
	grobB <- grobTree(textGrob("(B)", x=xloc,  y=yloc, gp=gpar(fontsize=fontsiz)))
	grobC <- grobTree(textGrob("(C)", x=xloc,  y=yloc, gp=gpar(fontsize=fontsiz)))
	grobD <- grobTree(textGrob("(D)", x=xloc,  y=yloc, gp=gpar(fontsize=fontsiz)))
	grobE <- grobTree(textGrob("(E)", x=xloc,  y=yloc, gp=gpar(fontsize=fontsiz)))
	grobF <- grobTree(textGrob("(F)", x=xloc,  y=yloc, gp=gpar(fontsize=fontsiz)))
	grobG <- grobTree(textGrob("(G)", x=xloc,  y=yloc, gp=gpar(fontsize=fontsiz)))
	grobH <- grobTree(textGrob("(H)", x=xloc,  y=yloc, gp=gpar(fontsize=fontsiz)))
	grobI <- grobTree(textGrob("(I)", x=xloc,  y=yloc, gp=gpar(fontsize=fontsiz)))


	ps1[[1]] <- ps1[[1]] + annotation_custom(grobA)
	ps1[[2]] <- ps1[[2]] + annotation_custom(grobB)
	ps1[[3]] <- ps1[[3]] + annotation_custom(grobC)

	ps3[[1]] <- ps3[[1]] + annotation_custom(grobD)
	ps3[[2]] <- ps3[[2]] + annotation_custom(grobE)
	ps3[[3]] <- ps3[[3]] + annotation_custom(grobF)

	ps2[[1]] <- ps2[[1]] + annotation_custom(grobG)
	ps2[[2]] <- ps2[[2]] + annotation_custom(grobH)
	ps2[[3]] <- ps2[[3]] + annotation_custom(grobI)


	#########################################################################


	ps1[[1]] <- ps1[[1]] + ggtitle("Simulated example")
	ps1[[1]] <- ps1[[1]] + theme(plot.title = element_text(size = 16, face = "plain", hjust = 0))

	ps2[[1]] <- ps2[[1]] + ggtitle("Rice example")
	ps2[[1]] <- ps2[[1]] + theme(plot.title = element_text(size = 16, face = "plain", hjust = 0))

	ps3[[1]] <- ps3[[1]] + ggtitle("Biting midges")
	ps3[[1]] <- ps3[[1]] + theme(plot.title = element_text(size = 16, face = "plain", hjust = 0))


	##########################################################################


	ps1[[3]] <- ps1[[3]] + ggtitle(textAnalysis)
	ps1[[3]] <- ps1[[3]] + theme(plot.title = element_text(size = 15, face = "plain", hjust = hjustTitle2))

	doSave <- saveInfo[[1]]
	locwd <- saveInfo[[2]]
	nameFile <- saveInfo[[3]]

	if(doSave)
	{
		setwd(locwd)
		pdf(nameFile,width=7, height=7)
			grid.arrange(
			ps1[[1]], ps1[[2]], ps1[[3]],
			ps3[[1]], ps3[[2]], ps3[[3]],
			ps2[[1]], ps2[[2]], ps2[[3]],
		nrow=3, ncol=3)
		dev.off() 
	}

	dev.new(width=7, height=7)
	grid.arrange(
		ps1[[1]], ps1[[2]], ps1[[3]],
		ps3[[1]], ps3[[2]], ps3[[3]],
		ps2[[1]], ps2[[2]], ps2[[3]],
	nrow=3, ncol=3)

}


#######################################################################
#---Two data examples
#---Ordination, contributions, correlations
#########################################################################



dataPics4 <- function(lraRice, lraKnut, saveInfo, whichContrib=c(1, 1), OrientOrdi = c(1, 1, 1, 1))
{
	#--adjust orientation Ordination diagram
	lraKnut$Sc[, 1] <- lraKnut$Sc[, 1] * OrientOrdi[1]
	lraKnut$Sc[, 2] <- lraKnut$Sc[, 2] * OrientOrdi[2]
	lraRice$Sc[, 1] <- lraRice$Sc[, 1] * OrientOrdi[3]
	lraRice$Sc[, 2] <- lraRice$Sc[, 2] * OrientOrdi[4]

	ps1 <- list()
	ps1[[1]] <- FigOrd(lraKnut)
	ps1[[2]] <- FigCor(lraKnut)
	ps1[[3]] <- FigContr(lraKnut, whichContrib=whichContrib[1])

	ps2 <- list()
	ps2[[1]] <- FigOrd(lraRice)
	ps2[[2]] <- FigCor(lraRice)
	ps2[[3]] <- FigContr(lraRice, whichContrib=whichContrib[2])


	#########################################################################


	#----Add letters

	yloc <- 0.95
	grobA <- grobTree(textGrob("(A)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))
	grobB <- grobTree(textGrob("(B)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))
	grobC <- grobTree(textGrob("(C)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))
	grobD <- grobTree(textGrob("(D)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))
	grobE <- grobTree(textGrob("(E)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))
	grobF <- grobTree(textGrob("(F)", x=0.05,  y=yloc, gp=gpar(fontsize=11)))

	ps1[[1]] <- ps1[[1]] + annotation_custom(grobA)
	ps1[[2]] <- ps1[[2]] + annotation_custom(grobB)
	ps1[[3]] <- ps1[[3]] + annotation_custom(grobC)
	ps2[[1]] <- ps2[[1]] + annotation_custom(grobD)
	ps2[[2]] <- ps2[[2]] + annotation_custom(grobE)
	ps2[[3]] <- ps2[[3]] + annotation_custom(grobF)


	#----Add titles

	ps1[[1]] <- ps1[[1]] + ggtitle("Biting midges example")
	ps1[[1]] <- ps1[[1]] + theme(plot.title = element_text(size = 18, face = "plain", hjust = 0))

	ps2[[1]] <- ps2[[1]] + ggtitle("Rice example")
	ps2[[1]] <- ps2[[1]] + theme(plot.title = element_text(size = 18, face = "plain", hjust = 0))


	##########################################################################


	doSave <- saveInfo[[1]]
	locwd <- saveInfo[[2]]
	nameFile <- saveInfo[[3]]

	if(doSave)
	{
		setwd(locwd)
		pdf(nameFile,width=12, height=7.5)
			grid.arrange(
			ps1[[1]], ps1[[2]], ps1[[3]],
			ps2[[1]], ps2[[2]], ps2[[3]],
		nrow=2, ncol=3)
		dev.off() 
	}

	dev.new(width=12, height=7.5)
	grid.arrange(
		ps1[[1]], ps1[[2]], ps1[[3]],
		ps2[[1]], ps2[[2]], ps2[[3]],
	nrow=2, ncol=3)
}



################################################################################################
#---Figure for simulation result.
#---Three simulations
#---plots ordination, contributions, and correlation
################################################################################################



simPic1 <- function(res1, res2, res3, saveInfo, OrientOrdi = c(1, 1, 1, 1, 1, 1))
{

	#--adjust orientation Ordination diagram
	res1$Sc[, 1] <- res1$Sc[, 1] * OrientOrdi[1]
	res1$Sc[, 2] <- res1$Sc[, 2] * OrientOrdi[2]
	res2$Sc[, 1] <- res2$Sc[, 1] * OrientOrdi[3]
	res2$Sc[, 2] <- res2$Sc[, 2] * OrientOrdi[4]
	res3$Sc[, 1] <- res3$Sc[, 1] * OrientOrdi[5]
	res3$Sc[, 2] <- res3$Sc[, 2] * OrientOrdi[6]


	pic1 <- list()
	pic1[[1]] <- FigOrd(res1)
	pic1[[2]] <- FigCor(res1)
	pic1[[3]] <- FigContr(res1)


	pic2 <- list()
	pic2[[1]] <- FigOrd(res2)
	pic2[[2]] <- FigCor(res2)
	pic2[[3]] <- FigContr(res2, whichContrib=2)


	pic3 <- list()
	pic3[[1]] <- FigOrd(res3)
	pic3[[2]] <- FigCor(res3)
	pic3[[3]] <- FigContr(res3)


	#--add x to contribution/correlation plot
	pic3[[2]] <- pic3[[2]] + xlab("log of taxon mean")
	pic3[[2]] <- pic3[[2]] + theme(axis.text.x = element_text(size=7))
	pic3[[3]] <- pic3[[3]] + xlab("log of taxon mean")
	pic3[[3]] <- pic3[[3]] + theme(axis.text.x = element_text(size=7))


	#--remove unneeded xlabels
	pic1[[2]] <- pic1[[2]] + xlab("")
	pic1[[3]] <- pic1[[3]] + xlab("")
	pic2[[2]] <- pic2[[2]] + xlab("")
	pic2[[3]] <- pic2[[3]] + xlab("")


	#--add y labels to contribution plot
	#pic1[[3]] <- pic1[[3]] + ylab("Log contrib 1st axis")
	#pic2[[3]] <- pic2[[3]] + ylab("Log contrib 2nd axis")
	#pic3[[3]] <- pic3[[3]] + ylab("Log contrib 1st axis")


	#--remove ticks contributions
	pic1[[2]] <- pic1[[2]] + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
	pic2[[2]] <- pic2[[2]] + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())


	#--remove ticks correlations
	pic1[[3]] <- pic1[[3]] + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
	pic2[[3]] <- pic2[[3]] + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())


	sdtext1 <- bquote(sigma[a] == .(sdSite[1]))
	sdtext2 <- bquote(sigma[a] == .(sdSite[2]))
	sdtext3 <- bquote(sigma[a] == .(sdSite[3]))

	pic1[[3]] <- pic1[[3]] + ggtitle(sdtext1)
	pic1[[3]] <- pic1[[3]] + theme(plot.title = element_text(size = 18, face = "plain", hjust = 0.9))
	pic2[[3]] <- pic2[[3]] + ggtitle(sdtext2)
	pic2[[3]] <- pic2[[3]] + theme(plot.title = element_text(size = 18, face = "plain", hjust = 0.9))
	pic3[[3]] <- pic3[[3]] + ggtitle(sdtext3)
	pic3[[3]] <- pic3[[3]] + theme(plot.title = element_text(size = 18, face = "plain", hjust = 0.9))


	bot <- 0
	right <- 0
	top <- 0
	left <- 0.25
	bot2 <- 0.25

	pic1[[1]] <- pic1[[1]] + theme(plot.margin = unit(c(top,right,bot,left), "cm")) 	#top, right, bottom, left
	pic1[[2]] <- pic1[[2]] + theme(plot.margin = unit(c(top,right,bot,left), "cm")) 	#top, right, bottom, left
	pic1[[3]] <- pic1[[3]] + theme(plot.margin = unit(c(top,right,bot,left), "cm")) 	#top, right, bottom, left
	pic2[[1]] <- pic2[[1]] + theme(plot.margin = unit(c(top,right,bot,left), "cm")) 	#top, right, bottom, left
	pic2[[2]] <- pic2[[2]] + theme(plot.margin = unit(c(top,right,bot,left), "cm")) 	#top, right, bottom, left
	pic2[[3]] <- pic2[[3]] + theme(plot.margin = unit(c(top,right,bot,left), "cm")) 	#top, right, bottom, left
	pic3[[1]] <- pic3[[1]] + theme(plot.margin = unit(c(top,right,bot2,left), "cm")) 	#top, right, bottom, left
	pic3[[2]] <- pic3[[2]] + theme(plot.margin = unit(c(top,right,bot2,left), "cm")) 	#top, right, bottom, left
	pic3[[3]] <- pic3[[3]] + theme(plot.margin = unit(c(top,right,bot2,left), "cm")) 	#top, right, bottom, left


	#--make grobs
	grobA <- grobTree(textGrob("(A)", x=0.05,  y=0.95, gp=gpar(fontsize=11)))
	grobB <- grobTree(textGrob("(B)", x=0.05,  y=1.05, gp=gpar(fontsize=11)))
	grobC <- grobTree(textGrob("(C)", x=0.05,  y=1.05, gp=gpar(fontsize=11)))
	grobD <- grobTree(textGrob("(D)", x=0.05,  y=0.95, gp=gpar(fontsize=11)))
	grobE <- grobTree(textGrob("(E)", x=0.05,  y=1.05, gp=gpar(fontsize=11)))
	grobF <- grobTree(textGrob("(F)", x=0.05,  y=1.05, gp=gpar(fontsize=11)))
	grobG <- grobTree(textGrob("(G)", x=0.05,  y=0.95, gp=gpar(fontsize=11)))
	grobH <- grobTree(textGrob("(H)", x=0.05,  y=1.05, gp=gpar(fontsize=11)))
	grobI <- grobTree(textGrob("(I)", x=0.05,  y=1.05, gp=gpar(fontsize=11)))


	#--add grobs of clipping
	pic1[[1]] <- pic1[[1]] + annotation_custom(grobA)
	pic1[[2]] <- pic1[[2]] + annotation_custom(grobB)
	pic1[[3]] <- pic1[[3]] + annotation_custom(grobC)
	pic2[[1]] <- pic2[[1]] + annotation_custom(grobD)
	pic2[[2]] <- pic2[[2]] + annotation_custom(grobE)
	pic2[[3]] <- pic2[[3]] + annotation_custom(grobF)
	pic3[[1]] <- pic3[[1]] + annotation_custom(grobG)
	pic3[[2]] <- pic3[[2]] + annotation_custom(grobH)
	pic3[[3]] <- pic3[[3]] + annotation_custom(grobI)


	#--turn of clipping
	pic1[[1]] <- ggplot_gtable(ggplot_build(pic1[[1]]))
	pic1[[1]]$layout$clip[pic1[[1]]$layout$name == "panel"] <- "off"
	pic1[[2]] <- ggplot_gtable(ggplot_build(pic1[[2]]))
	pic1[[2]]$layout$clip[pic1[[2]]$layout$name == "panel"] <- "off"
	pic1[[3]] <- ggplot_gtable(ggplot_build(pic1[[3]]))
	pic1[[3]]$layout$clip[pic1[[3]]$layout$name == "panel"] <- "off"

	pic2[[1]] <- ggplot_gtable(ggplot_build(pic2[[1]]))
	pic2[[1]]$layout$clip[pic2[[1]]$layout$name == "panel"] <- "off"
	pic2[[2]] <- ggplot_gtable(ggplot_build(pic2[[2]]))
	pic2[[2]]$layout$clip[pic2[[2]]$layout$name == "panel"] <- "off"
	pic2[[3]] <- ggplot_gtable(ggplot_build(pic2[[3]]))
	pic2[[3]]$layout$clip[pic2[[3]]$layout$name == "panel"] <- "off"

	pic3[[1]] <- ggplot_gtable(ggplot_build(pic3[[1]]))
	pic3[[1]]$layout$clip[pic3[[1]]$layout$name == "panel"] <- "off"
	pic3[[2]] <- ggplot_gtable(ggplot_build(pic3[[2]]))
	pic3[[2]]$layout$clip[pic3[[2]]$layout$name == "panel"] <- "off"
	pic3[[3]] <- ggplot_gtable(ggplot_build(pic3[[3]]))
	pic3[[3]]$layout$clip[pic3[[3]]$layout$name == "panel"] <- "off"



	#######################################################


	doSave <- saveInfo[[1]]
	locwd <- saveInfo[[2]]
	nameFile <- saveInfo[[3]]

	if(doSave)
	{
		setwd(locwd)

		pdf(nameFile, width=12, height=7.5)
		grid.arrange(
			pic1[[1]], pic1[[2]], pic1[[3]],
			pic2[[1]], pic2[[2]], pic2[[3]],
			pic3[[1]], pic3[[2]], pic3[[3]],
			widths = c(2.0, 2.6, 2.2),
		nrow=3, ncol=3)
		dev.off() 
	}

	dev.new(width=12, height=7.5)
	grid.arrange(
		pic1[[1]], pic1[[2]], pic1[[3]],
		pic2[[1]], pic2[[2]], pic2[[3]],
		pic3[[1]], pic3[[2]], pic3[[3]],
		widths = c(2.0, 2.6, 2.2),
	nrow=3, ncol=3)
}



################################################################################################
#---Function to generate simulated data
################################################################################################



dataSim <- function(nSites=50,
	nSpecies=500,
	sdSite=NULL,
	sdSpecies=NULL,
	weightResp=NULL,
	Max=NULL,
	drawDistr=NULL,
	strengthCorrResp=0,
	ZI=0,
	zeroThres=5,
	nSigGenes=100)
{

	#---get LP
	siteEffect <- rnorm(nSites, sd=sdSite)
	speciesEffect <- rnorm(nSpecies, sd=sdSpecies)
	LP <- outer(siteEffect, speciesEffect, FUN="+")


	#---draw a response
	if(strengthCorrResp > 0 & sdSite > 0)
	{
		draw <- runif(nSites)
		plogisDraw <- plogis(strengthCorrResp * siteEffect)

		fTemp <- function(possibleScales)
		{
			nScale <- length(possibleScales)
			drawMat <- matrix(rep(draw, times=nScale), nrow=nSites, ncol=nScale)
			scalesMat <- matrix(rep(possibleScales, each=nSites), nrow=nSites, ncol=nScale)
			drawMatScaled <- drawMat * scalesMat
			isEqual <- plogisDraw > drawMatScaled
			w <- which(colSums(isEqual)==nSites/2)

			if(length(w)==1) resp <- as.numeric(isEqual[, w])
			if(length(w)>1) resp <- as.numeric(isEqual[, w[ceiling(length(w)/2)]])
			if(length(w)==0) resp <- NULL

			return(resp)
		}

		response <- fTemp(possibleScales=seq(from=0.25, to=1.75, by=0.1))
		if(is.null(response)) response <- fTemp(possibleScales=seq(from=0.01, to=2, by=0.01))
		if(is.null(response)) response <- fTemp(possibleScales=seq(from=0.01, to=2, by=0.001))
		if(is.null(response)) response <- fTemp(possibleScales=seq(from=0.01, to=2, by=0.0001))
		if(is.null(response)) response <- fTemp(possibleScales=seq(from=0.01, to=2, by=0.00001))
		if(is.null(response)) response <- fTemp(possibleScales=seq(from=0.01, to=2, by=0.000001))
		if(is.null(response))
		{
			draw <- runif(nSites)
			plogisDraw <- plogis(strengthCorrResp * siteEffect)
			response <- fTemp(possibleScales=seq(from=0.25, to=1.75, by=0.1))
			if(is.null(response)) response <- fTemp(possibleScales=seq(from=0.01, to=2, by=0.01))
			if(is.null(response)) response <- fTemp(possibleScales=seq(from=0.01, to=2, by=0.001))
			if(is.null(response)) response <- fTemp(possibleScales=seq(from=0.01, to=2, by=0.0001))
			if(is.null(response)) response <- fTemp(possibleScales=seq(from=0.01, to=2, by=0.00001))
			if(is.null(response)) response <- fTemp(possibleScales=seq(from=0.01, to=2, by=0.000001))
		}

	} else {
		#response <- rep(c(0,1), each=nSites/2)
		response <- rep(c(0,1), each=nSites/2)[sample(nSites)]
	}

	change <- sample(c(-1,1), nSigGenes, replace=TRUE) * log(weightResp)
	changeMat <- t(matrix(change, nrow=nSigGenes, ncol=(nSites/2)))
	sigGenes <- sort(sample(1:nSpecies, nSigGenes, replace=FALSE))
	LP[response==1, sigGenes] <- LP[response==1, sigGenes] + changeMat
	LPexp <- Max * exp(LP)


	#################################################
	#--draw from the LPexp matrix
	#################################################


	#--draw with poisson
	if(drawDistr=="poiss") Draw <- rpois(length(LPexp), lambda= LPexp)

	#--draw with negative binomial
	if(drawDistr=="nb") Draw <- rnbinom(length(LPexp), size = 1, mu = LPexp)

	#--make a matrix again
	xtmp <- matrix(Draw,nrow = nSites,ncol = nSpecies)

	#--attach numbers
	rownames(xtmp) <- 1:nSites
	colnames(xtmp) <- 1:nSpecies

	#--extra xero inflation
	ZImat <- matrix(runif(nSites*nSpecies), nrow = nSites, ncol = nSpecies)
	xtmp[ZImat<ZI] <- 0


	#################################################
	#--Output
	#################################################


	#--remove rows and cols
	RemoveCols <- colSums(xtmp!=0) >= zeroThres
	xtmp2 <- xtmp[, RemoveCols]

	changeAll <- c(change, rep(0, times=nSpecies-nSigGenes))
	out <- list(xs=xtmp2, response=response, LP=LP, siteEffect=siteEffect, speciesEffect=speciesEffect[RemoveCols], changeAll=changeAll, sigGenes=sigGenes, sigGenes2 = which((1:nSpecies %in% sigGenes)[RemoveCols]))

	return(out)
}


###################################################################################



wLRA.ade4 <- function(x, response, RowWeights=NULL, ColWeights=NULL)
{
	#--take the log
	xlog <- log(x + 1)
	if(is.null(RowWeights)) RowWeights <- rowSums(xlog)
	if(is.null(ColWeights)) ColWeights <- colSums(xlog)

	#--second part LRA transform
	Xtr <- DoubleCenter(xlog, RowWeights=RowWeights, ColWeights=ColWeights)

	#---intial pca, because ade4 wants so. Output not used.
	pca1 <- dudi.pca(as.data.frame(Xtr), scannf = FALSE, scale = FALSE, center = FALSE, row.w = RowWeights, col.w = ColWeights)
	#note: pca1$tab == Xtr

	#--save eigenvalues
	eigs2 <- pca1$eig/sum(pca1$eig)

	SitesScores <- pca1$li
	principalSpecies <- pca1$co
	standardSpecies <- pca1$c1

	#--do svd for contribution as this is not provided by ade4
	Xsvd <- diag(sqrt(attr(Xtr,"RowWeights"))) %*% Xtr %*% diag(sqrt(attr(Xtr,"ColWeights")))
	svdOut <- svd(Xsvd)
	Contr <- log(svdOut$v[, 1:2]^2)
	#--check svd
	#((svdOut$d^2)/sum(svdOut$d^2))[1:5]-eigs2[1:5]


	#---Correlation between matrix that is input for SVD and its rowmeans
	Corrs <- cor(Xtr, attributes(Xtr)$rowmeans)
	#---Correlation with x
	Corrs2 <- cor(x, attributes(Xtr)$rowmeans)
	#--more output
	logColmeans <- log(colMeans(x))


	######################################


	ExplVar <- data.frame(y = attributes(Xtr)$rowmeans)
	rda1 <- pcaiv(dudi=pca1, df=ExplVar, scannf = FALSE)
	#rda1$tab is predMat

	#--calculate contributions r-contrained axis
	svdPred <- svd(rda1$tab)
	ContribConstr <- log(svdPred$v[,1]^2)
	residMat <- Xtr - rda1$tab
	pca2 <- dudi.pca(residMat, scannf = FALSE, scale = FALSE, center = FALSE, nf=ncol(x), row.w = RowWeights, col.w = ColWeights)
	
	#--save eigenvalues
	eigsContrained <- c(rda1$eig, pca2$eig)
	eigsContrained2 <- eigsContrained/sum(eigsContrained)


	######################################


	out <- list(Sc=SitesScores[, 1:2], Contr=Contr, logColmeans=logColmeans, Corrs=Corrs, Corrs2=Corrs2, resp=response, means=attributes(Xtr)$rowmeans, ContrConstr=ContribConstr, eigUncontr=eigs2, rContr=eigsContrained2, PrincipalSpecies=principalSpecies[,1:2])
	return(out)
}


###################################################################################


LogFrac.ade4 <- function(x, response, pseudoFrac)
{
	#x <- outxs[[2]]$x
	
	Xtr <- log((x/rowSums(x)) + pseudoFrac)
	pca1 <- dudi.pca(as.data.frame(Xtr), scannf = FALSE, scale = FALSE, center = TRUE)
	Sc <- pca1$li

	#--compare
	#rdaout <- rda(X=Xtr)
	#pca1$eig/sum(pca1$eig)
	#summary(rdaout)$cont$importance[2, 1:5]
	#plot(rdaout)
	#plot(pca1$li[, 1], -1*pca1$li[, 2])


	#--do svd for contribution as this is not provided by ade4.
	svdOut <- svd(scale(Xtr, scale=FALSE))
	#((svdOut$d^2)/sum(svdOut$d^2))[1:5]
	Contr <- log(svdOut$v[, 1:2]^2)

	logColmeans <- log(colMeans(x))

	#---Correlation with transformed x
	Corrs <- cor(Xtr, rowMeans(Xtr))
	#---Correlation with x
	Corrs2 <- cor(x, rowMeans(Xtr))

	out <- list(Sc=Sc, Contr=Contr, logColmeans=logColmeans, resp=response, Corrs=Corrs, means=rowMeans(Xtr), Corrs2=Corrs2)
	return(out)
}


###################################################################################


CA.ade4 <- function(x, response, xReal=NULL)
{
	res.ca <- dudi.coa(x, scannf = FALSE)
	Sc <- res.ca$li

	Xsvd <- CATransform(x)
	svdOut <- svd(Xsvd)
	Contr <- log(svdOut$v[, 1:2]^2)

	#--compare methods
	#Xsvd <- diag(sqrt(as.numeric(res.ca$lw))) %*% as.matrix(res.ca$tab) %*% diag(sqrt(as.numeric(res.ca$cw)))
	#svdOut <- svd(Xsvd)
	#ccaout <- cca(X=x)
	#(res.ca$eig/sum(res.ca$eig))[1:5]
	#summary(ccaout)$cont$importance[2, 1:5]
	#((svdOut$d^2)/sum(svdOut$d^2))[1:5]
	#plot(ccaout)
	#plot(res.ca$li[, 1], res.ca$li[, 2])
	#scale(Sc[, 1]) - scale(scores(ccaout)$sites[, 1])

	#---correlation between matrix that is input for SVD and its rowmeans
	#Corrs <- cor(Xsvd, rowMeans(x))
	Corrs <- cor(Xsvd, attr(Xsvd, "rowsum"))

	#---Correlation with x
	Corrs2 <- cor(x, rowMeans(x))

	#---output
	if(is.null(xReal)) logColmeans <- log(colMeans(x))
	if(!is.null(xReal)) logColmeans <- log(colMeans(xReal))

	out <- list(Sc=Sc, Contr=Contr, logColmeans=logColmeans, Corrs=Corrs, Corrs2=Corrs2, resp=response, means=rowMeans(x), Corrs2=Corrs2)
	return(out)
}


###################################################################################



LRA.ade4 <- function(x, response=NULL, pseudocount=1, xReal=NULL)
{
	if(is.null(xReal)) xReal <- x

	Xtr <- DoubleCenter(log(x + pseudocount))
	pca1 <- dudi.pca(as.data.frame(Xtr), scannf = FALSE, scale = FALSE, center = FALSE)
	Sc <- pca1$li

	#--do svd for contribution as this is not provided by ade4.
	svdOut <- svd(Xtr)
	#((svdOut$d^2)/sum(svdOut$d^2))[1:5]
	Contr <- log(svdOut$v[, 1:2]^2)

	#--compare
	#rdaout <- rda(X=Xtr)
	#pca1$eig/sum(pca1$eig)
	#summary(rdaout)$cont$importance[2, 1:5]
	#plot(rdaout)
	#plot(pca1$li[, 1], -1*pca1$li[, 2])
	#((svdOut$d^2)/sum(svdOut$d^2))[1:5]


	#---Correlation between matrix that is input for SVD and its rowmeans
	Corrs <- cor(Xtr, attributes(Xtr)$rowmeans)

	#--r contrained
	rda1 <- pcaiv(dudi=pca1, df=attributes(Xtr)$rowmeans, scannf = FALSE)
	residMat <- Xtr - rda1$tab
	pca2 <- dudi.pca(residMat, scannf = FALSE, scale = FALSE, center = FALSE, nf=ncol(x))
	svdOut <- svd(rda1$tab)
	ContrConstr <- log(svdOut$v[, 1]^2)

	rVariance <- rda1$eig/sum(c(rda1$eig, pca2$eig))


	#--compare
	#--to vegan
	#rdaout2 <- rda(Xtr ~ attributes(Xtr)$rowmeans)
	#summary(rdaout2)$cont$importance[2, 1:5]
	#--to plain svd
	#svdOut2 <- svd(residMat)
	#(c(svdOut$d[1]^2, svdOut2$d^2)/sum(c(svdOut$d^2, svdOut2$d^2)))[1:5]


	#---output
	logColmeans <- log(colMeans(xReal))
	out <- list(Sc=Sc, Contr=Contr, logColmeans=logColmeans, Corrs=Corrs, resp=response, means=attributes(Xtr)$rowmeans, ContrConstr=ContrConstr, rVariance=rVariance)
	return(out)
}



#########################################################################
#----pic for RCM
#########################################################################



picRCM <- function(outSim1, outSim2, outRice1, outRice2, outKnut1, outKnut2, saveInfo, OrientOrdi = c(1, 1, 1, 1, 1, 1))
{

	#########################################################################

	#--adjust orientation Ordination diagram
	outSim1$Sc[, 1] <- outSim1$Sc[, 1] * OrientOrdi[1]
	outSim1$Sc[, 2] <- outSim1$Sc[, 2] * OrientOrdi[2]

	outSim2$Sc[, 1] <- outSim2$Sc[, 1] * OrientOrdi[3]
	outSim2$Sc[, 2] <- outSim2$Sc[, 2] * OrientOrdi[4]

	outRice1$Sc[, 1] <- outRice1$Sc[, 1] * OrientOrdi[5]
	outRice1$Sc[, 2] <- outRice1$Sc[, 2] * OrientOrdi[6]

	outRice2$Sc[, 1] <- outRice2$Sc[, 1] * OrientOrdi[7]
	outRice2$Sc[, 2] <- outRice2$Sc[, 2] * OrientOrdi[8]

	outKnut1$Sc[, 1] <- outKnut1$Sc[, 1] * OrientOrdi[9]
	outKnut1$Sc[, 2] <- outKnut1$Sc[, 2] * OrientOrdi[10]

	outKnut2$Sc[, 1] <- outKnut2$Sc[, 1] * OrientOrdi[11]
	outKnut2$Sc[, 2] <- outKnut2$Sc[, 2] * OrientOrdi[12]



	#########################################################################


	ps1 <- list()

	ps1[[1]] <- FigOrd(outSim1)
	ps1[[2]] <- FigOrd(outRice1)
	ps1[[3]] <- FigOrd(outKnut1)
	ps1[[4]] <- FigOrd(outSim2)
	ps1[[5]] <- FigOrd(outRice2)
	ps1[[6]] <- FigOrd(outKnut2)


	#########################################################################

	
	yloc <- 0.95
	xloc <- 0.05
	xloc2 <- 0.05
	grobA <- grobTree(textGrob("(A)", x=xloc,  y=yloc, gp=gpar(fontsize=11)))
	grobB <- grobTree(textGrob("(B)", x=xloc,  y=yloc, gp=gpar(fontsize=11)))
	grobC <- grobTree(textGrob("(C)", x=xloc2,  y=yloc, gp=gpar(fontsize=11)))
	grobD <- grobTree(textGrob("(D)", x=xloc,  y=yloc, gp=gpar(fontsize=11)))
	grobE <- grobTree(textGrob("(E)", x=xloc,  y=yloc, gp=gpar(fontsize=11)))
	grobF <- grobTree(textGrob("(F)", x=xloc2,  y=yloc, gp=gpar(fontsize=11)))


	ps1[[1]] <- ps1[[1]] + annotation_custom(grobA)
	ps1[[4]] <- ps1[[4]] + annotation_custom(grobB)

	ps1[[3]] <- ps1[[3]] + annotation_custom(grobC)
	ps1[[6]] <- ps1[[6]] + annotation_custom(grobD)

	ps1[[2]] <- ps1[[2]] + annotation_custom(grobE)
	ps1[[5]] <- ps1[[5]] + annotation_custom(grobF)

	##########################################################################


	ps1[[1]] <- ps1[[1]] + ggtitle("Simulated example")
	ps1[[2]] <- ps1[[2]] + ggtitle("Rice example")
	ps1[[3]] <- ps1[[3]] + ggtitle("Biting midges example")

	ps1[[4]] <- ps1[[4]] + ggtitle("With additional filtering")
	ps1[[5]] <- ps1[[5]] + ggtitle("With additional filtering")
	ps1[[6]] <- ps1[[6]] + ggtitle("With additional filtering")

	ps1[[1]] <- ps1[[1]] + theme(plot.title = element_text(size = 15, face = "plain", hjust = 0.05))
	ps1[[2]] <- ps1[[2]] + theme(plot.title = element_text(size = 15, face = "plain", hjust = 0.05))
	ps1[[3]] <- ps1[[3]] + theme(plot.title = element_text(size = 15, face = "plain", hjust = 0.05))


	##########################################################################


	doSave <- saveInfo[[1]]
	locwd <- saveInfo[[2]]
	nameFile <- saveInfo[[3]]

	if(doSave)
	{
		setwd(locwd)
		pdf(nameFile,width=6, height=9)
			grid.arrange(
			ps1[[1]], ps1[[4]],
			ps1[[3]], ps1[[6]],
			ps1[[2]], ps1[[5]],
		nrow=3, ncol=2)
		dev.off() 
	}

	dev.new(width=6, height=9)
		grid.arrange(
			ps1[[1]], ps1[[4]],
			ps1[[3]], ps1[[6]],
			ps1[[2]], ps1[[5]],
	nrow=3, ncol=2)
}



##########################################################################
#---transformation for logratio pca
##########################################################################



DoubleCenter <-function(X, RowWeights = NULL, ColWeights = NULL)
{
	if(is.null(RowWeights)) RowWeights <- rep(1,times=dim(X)[1])
	if(is.null(ColWeights)) ColWeights <- rep(1,times=dim(X)[2])

	RowWeightsRel <- RowWeights/sum(RowWeights)
	ColWeightsRel <- ColWeights/sum(ColWeights)
	TWeightsRel <- outer(RowWeightsRel,ColWeightsRel)

	rowmeans <- rowSums(X %*% diag(ColWeightsRel))
	colmeans <- colSums(diag(RowWeightsRel) %*% X)
	globalmean <- sum(TWeightsRel * X)

	Rows <- X*0 + rowmeans
	Cols <- t(t(X)*0 + colmeans)
	X_double_centered = X - Rows - Cols + globalmean

	attr(X_double_centered,"rowmeans") <- rowmeans
	attr(X_double_centered,"colmeans") <- colmeans
	attr(X_double_centered,"globalmean") <- globalmean
	attr(X_double_centered,"RowWeights") <- RowWeights
	attr(X_double_centered,"ColWeights") <- ColWeights
	attr(X_double_centered,"RowWeightsRel") <- RowWeightsRel
	attr(X_double_centered,"ColWeightsRel") <- ColWeightsRel
	attr(X_double_centered,"Analysis") <- "WPCA"

	X_double_centered
}


##########################################################################
#--Transformation for CA
##########################################################################



CATransform <-function(Counts)
{
	globalsum <- sum(Counts)
	SMat <- Counts/globalsum
	rowsum <- rowSums(SMat)
	colsum <- colSums(SMat)
	smat1 <- matrix(rowsum,nrow=nrow(SMat),ncol=ncol(SMat))
	smat2 <- t(matrix(colsum,ncol=nrow(SMat),nrow=ncol(SMat)))
	CAMat <- (SMat-smat1*smat2)/sqrt(smat1*smat2)

	attr(CAMat,"rowsum") <- rowsum
	attr(CAMat,"colsum") <- colsum
	attr(CAMat,"globalsum") <- globalsum
	attr(CAMat,"Analysis") <- "CA"

	return(CAMat)
}



