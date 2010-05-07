###################################################
### chunk number 1: setupCols
###################################################
library("RbcBook1")
library("RColorBrewer")
mypal<- brewer.pal(7, "RdBu")
blue <- mypal[7]; red<-mypal[1]


###################################################
### chunk number 2: affydata
###################################################
library("affydata")
data("Dilution")
x <- log2(exprs(Dilution)[, 1:2]) 
x <- x %*% cbind(A=c(1,1), M=c(-1,1))


###################################################
### chunk number 3: scatterplot-subsample
###################################################
par(mfrow=c(1,1))
x <- x[sample(nrow(x), 5e4),  ]


###################################################
### chunk number 4: scatterplot
###################################################
plot(x, pch=".")


###################################################
### chunk number 5: hexbin4
###################################################
library("hexbin")
library("geneplotter")
hb <- hexbin(x, xbins=50)
plot(hb, colramp=colorRampPalette(brewer.pal(9,"YlGnBu")[-1]))


###################################################
### chunk number 6: hexbin5
###################################################
library("prada")
smoothScatter(x, nrpoints=500)


###################################################
### chunk number 7: hexbin6
###################################################
plot(x, col=densCols(x), pch=20)


###################################################
### chunk number 8: ALLh1
###################################################
library("ALL")
data("ALL")
selSamples <- ALL$mol.biol %in% c("ALL1/AF4", "E2A/PBX1")
ALLs <- ALL[, selSamples]
ALLs$mol.biol <- factor(ALLs$mol.biol)
colnames(exprs(ALLs)) <- paste(ALLs$mol.biol, colnames(exprs(ALLs)))


###################################################
### chunk number 9: ALLh2
###################################################
library("genefilter")
meanThr <- log2(100)
g       <- ALLs$mol.biol
s1 <- rowMeans(exprs(ALLs)[, g==levels(g)[1]]) > meanThr
s2 <- rowMeans(exprs(ALLs)[, g==levels(g)[2]]) > meanThr
s3 <- rowttests(ALLs, g)$p.value < 0.0002
selProbes <- (s1 | s2) & s3
ALLhm <- ALLs[selProbes,]


###################################################
### chunk number 10: ALLheatmap
###################################################
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
spcol <- ifelse(ALLhm$mol.biol=="ALL1/AF4", "goldenrod", "skyblue")
heatmap(exprs(ALLhm), col=hmcol, ColSideColors=spcol)


###################################################
### chunk number 11: estrogenProcess
###################################################
esEset <- cache("metaVisualize-esEset", {
  library("estrogen")
  library("limma")
  library("hgu95av2cdf")
  
  datadir <- system.file("extdata", package="estrogen")
  targets <- readTargets("phenoData.txt",path=datadir,sep="")

  covdesc <- list("present or absent","10 or 48 hours")
  names(covdesc) <- names(targets)[-1]
  pdata <- new("phenoData",pData=targets[,-1],varLabels=covdesc)
  rownames(pData(pdata)) <- targets[,1]
  
  esAB <- ReadAffy(filenames=file.path(datadir,targets$filename),phenoData=pdata)
  esEset <- rma(esAB)
})
  
fit <- cache("metaVisualize-fit", {
  pdat <- pData(esEset)
  design <- model.matrix(~factor(estrogen)*factor(time.h),pdat)
  colnames(design) <- c("Intercept","ES","T48","ES:T48")
  lmFit(esEset,design) 
})
stopifnot(all(fit$df.residual==4))  #$

colnames(exprs(esEset)) <- paste(
  c("-", "+")[match(esEset$estrogen, c("absent", "present"))], esEset$time.h)


###################################################
### chunk number 12: predictMArrayLM
###################################################
predict.MArrayLM <- function(f, design=f$design) {
  return(f$coefficients %*% t(design))
}  
esFit <- predict(fit)
res <- exprs(esEset) - esFit


###################################################
### chunk number 13: estroHM1
###################################################
sel <- order(fit$coefficients[, "ES:T48"], decreasing=TRUE)[1:50]
four.groups <- as.integer(factor(colnames(exprs(esEset))))
csc <- brewer.pal(4, "Paired")[four.groups]
heatmap(exprs(esEset)[sel,], col=hmcol, ColSideColors=csc) 


###################################################
### chunk number 14: estroHM2
###################################################
heatmap(res[sel,], col=hmcol, ColSideColors=csc) 


###################################################
### chunk number 15: ALLhc1
###################################################
standardize <- function(z) {
  rowmed <- apply(z, 1, median)
  rowmad <- apply(z, 1, mad)
  rv <- sweep(z, 1, rowmed)
  rv <- sweep(rv, 1, rowmad, "/")
  return(rv)
}

ALLhme <- exprs(ALLhm)
ALLdist1 <- dist(t(standardize(ALLhme)))
ALLhc1 <- hclust(ALLdist1)


###################################################
### chunk number 16: ALLdendro1
###################################################
plot(ALLhc1, xlab="", sub="", main="ALLhc1")


###################################################
### chunk number 17: ALLhc2
###################################################
ALLsub2 <- exprs(ALLs[(s1 | s2), ])
rowMads <- apply(ALLsub2, 1, mad)
ALLsub2 <- ALLsub2[rowMads > 1.4, ]

ALLdist2 <- dist(t(standardize(ALLsub2)))
ALLhc2 <- hclust(ALLdist2)


###################################################
### chunk number 18: ALLdendro2
###################################################
plot(ALLhc2, xlab="", sub="", main="ALLhc2")


###################################################
### chunk number 19: cophenetic1
###################################################
ALLcph1 <- cophenetic(ALLhc1)
cor(ALLdist1, ALLcph1)
plot(ALLdist1, ALLcph1, pch="|", col=blue)


###################################################
### chunk number 20: cophenetic2
###################################################
ALLcph2 <- cophenetic(ALLhc2)
cor(ALLdist2,ALLcph2)
plot(ALLdist2,ALLcph2,pch="|", col=blue)


###################################################
### chunk number 21: cmdScale
###################################################
library(MASS)
cm1 <- cmdscale(ALLdist1, eig=TRUE)
cm1$GOF
samm1 <- sammon(ALLdist1, trace=FALSE)
cm2 <- cmdscale(ALLdist2, eig=TRUE)
cm2$GOF
samm2 <- sammon(ALLdist2, trace=FALSE)


###################################################
### chunk number 22: myPlot eval=FALSE
###################################################
## ALLscol <- c("goldenrod", "skyblue")[as.integer(ALLs$mol.biol)]
## plot(cm1$points, col=ALLscol, ...)


###################################################
### chunk number 23: cmdPlots
###################################################
ALLscol <- c("goldenrod", "skyblue")[as.integer(ALLs$mol.biol)]
myPlot <- function(x, ...)
  plot(x$points, xlab="Component 1", ylab="Component 2", pch=19,  col=ALLscol, ...)
par(mfrow=c(2,2))
myPlot(cm1, main="a) metric / t-test")
myPlot(samm1, main="b) Sammon / t-test")
myPlot(cm2, main="c) metric / MAD")
myPlot(samm2, main="d) Sammon / MAD")
par(mfrow=c(1,1))


###################################################
### chunk number 24: distHM
###################################################
heatmap(as.matrix(ALLdist2), sym=TRUE, col=hmcol,
    distfun=function(x) as.dist(x))


###################################################
### chunk number 25: chrLoc
###################################################
library("geneplotter")
chrLoc <- buildChromLocation("hgu95av2")
## cPlot(chrLoc, main="HG-U95Av2")


###################################################
### chunk number 26: cPlot2a
###################################################
ALLch <- ALLs[s1|s2, ]
m1  <- rowMeans(exprs(ALLch)[, ALLch$mol.biol=="ALL1/AF4"])
m2  <- rowMeans(exprs(ALLch)[, ALLch$mol.biol=="E2A/PBX1"])


###################################################
### chunk number 27: cPlot2b
###################################################
deciles <- quantile(c(m1,m2), probs=seq(0,1,.1))
s1dec <- cut(m1, deciles)
s2dec <- cut(m2, deciles)
gN <- names(s1dec) <- names(s2dec) <- geneNames(ALLch)


###################################################
### chunk number 28: cPlot3
###################################################
colors <- brewer.pal(10, "RdBu")
layout(matrix(1:3,nr=1), widths=c(5,5,2))
cPlot(chrLoc, main="ALL1/AF4")
cColor(gN, colors[s1dec], chrLoc)
cPlot(chrLoc, main="E2A/PBX1")
cColor(gN, colors[s2dec], chrLoc)
image(1,1:10,matrix(1:10,nc=10),col=colors, axes=FALSE,
         xlab="", ylab="")
axis(2, at=(1:10), labels=levels(s1dec), las=1)


###################################################
### chunk number 29: plotChr
###################################################
par(mfrow=c(1,1))
msobj <- Makesense(ALLs, "hgu95av2")
plotChr("22", msobj, 
   col = ifelse(ALLs$mol.biol=="ALL1/AF4", "#EF8A62", "#67A9CF"), log=FALSE)


