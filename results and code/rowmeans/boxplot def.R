

gitdit <- "E:/1Documenten/Werk/gits/gitlabwur/projects/3microbiome/LRApaper/results def"
#gitdit <- "C:/1Files/2Git/smallProjects/microbiome/LRApaper/results def"
setwd(gitdit)
library(vegan)

load("dataMidges.RData")
load("dataRice.RData")
ls()




#################################################################################################



set.seed(1)
x1Rar <- rrarefy(dataRice$x, min(rowSums(dataRice$x)))
xRiceRar <- x1Rar[, colSums(x1Rar != 0) >= 10]

set.seed(1)
x2Rar <- rrarefy(dataMidges$x, min(rowSums(dataMidges$x)))
xMidgeRar <- x2Rar[, colSums(x2Rar != 0) >= 10]



pdf("rowmeans/rowmeansDataRar.pdf",width=8, height=8)


par(mfrow=c(2, 2))
par(mar=c(3, 4, 2, 1))   #c(bottom, left, top, right)
par(mgp=c(2.5, 1, 0))
par(oma=c(0, 0, 0, 0))

boxplot(rowMeans(log(dataMidges$x + 1)) ~ dataMidges$response, ylab="mean log transformed abundance (r)", main="Biting midges example", names=c("Sugarwater", "Antibiotics"), xlab="", cex.axis=1.0, cex.lab=1.1)
text(x = 0.55, y = 1.61, labels = "(A)")
boxplot(rowMeans(log(xMidgeRar + 1)) ~ dataMidges$response, ylab="mean log transformed abundance (r)", main="Biting midges example - Rarefied", names=c("Sugarwater", "Antibiotics"), xlab="", cex.axis=1.0, cex.lab=1.1)
text(x = 0.55, y = 0.93, labels = "(B)")
boxplot(rowMeans(log(dataRice$x + 1)) ~ dataRice$treatment, ylab="mean log transformed abundance (r)", main="Rice example", names=c("Treatment 1","Treatment 2"), xlab="", cex.axis=1.0, cex.lab=1.1)
text(x = 0.55, y = 0.352, labels = "(C)")
boxplot(rowMeans(log(xRiceRar + 1)) ~ dataRice$treatment, ylab="mean log transformed abundance (r)", main="Rice example - Rarefied", names=c("Treatment 1","Treatment 2"), xlab="", cex.axis=1.0, cex.lab=1.1)
text(x = 0.55, y = 0.797, labels = "(D)")



dev.off() 



#################################################################################################



pdf("rowmeans/rowmeansLibsizeData.pdf",width=8, height=8)


par(mfrow=c(2, 2))
par(mar=c(3, 4, 2, 1))   #c(bottom, left, top, right)
par(mgp=c(2.5, 1, 0))
par(oma=c(0, 0, 0, 0))

boxplot(rowSums(dataMidges$x) ~ dataMidges$response, ylab="Library size", main="Biting midges example", names=c("Sugarwater", "Antibiotics"), xlab="", cex.axis=1.0, cex.lab=1.3)
text(x = 0.55, y = 128000, labels = "(A)")
boxplot(rowMeans(log(dataMidges$x + 1)) ~ dataMidges$response, ylab="mean log transformed abundance (r)", main="Biting midges example", names=c("Sugarwater", "Antibiotics"), xlab="", cex.axis=1.0, cex.lab=1.1)
text(x = 0.55, y = 1.61, labels = "(B)")
boxplot(rowSums(dataRice$x) ~ dataRice$treatment, ylab="Library size", main="Rice example", names=c("Treatment 1","Treatment 2"), xlab="", cex.axis=1.0, cex.lab=1.3)
text(x = 0.55, y = 92000, labels = "(C)")
boxplot(rowMeans(log(dataRice$x + 1)) ~ dataRice$treatment, ylab="mean log transformed abundance (r)", main="Rice example", names=c("Treatment 1","Treatment 2"), xlab="", cex.axis=1.0, cex.lab=1.1)
text(x = 0.55, y = 0.352, labels = "(D)")


dev.off() 






