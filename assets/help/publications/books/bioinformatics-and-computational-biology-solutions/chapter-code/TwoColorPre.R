###################################################
### chunk number 1: setup
###################################################
library("RbcBook1")


###################################################
### chunk number 2: checkVersions
###################################################
library("arrayQuality")
library("marray")
library("beta7")
stopifnot(package.version("arrayQuality") >= package_version("1.0.9"))
stopifnot(package.version("marray") >= package_version("1.5.29"))
stopifnot(package.version("beta7") >= package_version("0.5.4"))


###################################################
### chunk number 3: GEODo
###################################################
library("AnnBuilder")
samp.6Hs.166 <- cache("samp.6Hs.166", 
   queryGEO(GEO(), "GSM16689"))


###################################################
### chunk number 4: GEOshow eval=FALSE
###################################################
## library("AnnBuilder")
## samp.6Hs.166 <- queryGEO(GEO(), "GSM16689") 


###################################################
### chunk number 5: readbeta7
###################################################
datadir <- system.file("beta7", package="beta7")
TargetInfo <- read.marrayInfo(file.path(datadir, "TargetBeta7.txt"))


###################################################
### chunk number 6: info1
###################################################
TargetInfo@maNotes <- "Files were loaded from beta7 package."


###################################################
### chunk number 7: info2
###################################################
TargetInfo


###################################################
### chunk number 8: Kote13
###################################################
galinfo <- read.Galfile("6Hs.166.gpr", path=datadir)


###################################################
### chunk number 9: oldwd1
###################################################
oldwd <- getwd()
setwd(datadir)


###################################################
### chunk number 10: read.GenePix
###################################################
setwd(datadir)
files <- c("6Hs.166.gpr", "6Hs.187.1.gpr")
mraw <- read.GenePix(files, name.Gb=NULL, name.Rb=NULL)


###################################################
### chunk number 11: oldwd2
###################################################
setwd(oldwd)


###################################################
### chunk number 12: gnurps3
###################################################
library("beta7")
checkTargetInfo(beta7)


###################################################
### chunk number 13: maGeneTable
###################################################
maGeneTable(beta7)[1:4, 1:5]


###################################################
### chunk number 14: whatAnUglyHack
###################################################
beta7nbg <- beta7
beta7nbg@maGb <- beta7nbg@maRb <- 0 * beta7nbg@maRb


###################################################
### chunk number 15: subsett
###################################################
beta7sub <- beta7[1:100,2:3]


###################################################
### chunk number 16: subsetu
###################################################
coord <- maCompCoord(1:maNgr(beta7), 1:maNgc(beta7), maNsr(beta7), 1:3)
ind <- maCoord2Ind(coord, L=maLayout(beta7))


###################################################
### chunk number 17:  eval=FALSE
###################################################
## maQualityPlots(beta7)


###################################################
### chunk number 18:  eval=FALSE
###################################################
## agQuality()


###################################################
### chunk number 19: ZZ1
###################################################
image(beta7[,5], xvar = "maRb",  bar = TRUE)


###################################################
### chunk number 20: ZZ2
###################################################
RGcol <- maPalette(low = "blue", mid = "gray", high = "yellow", k = 50)
image(beta7[, 3], xvar = "maM", col=RGcol)


###################################################
### chunk number 21: ZZ3
###################################################
flags <-  beta7@maW[,1] < -50
image(beta7[,1], xvar="maA", overlay=flags)


###################################################
### chunk number 22: maBoxplotplate
###################################################
par(mar=c(5, 3,3,3), cex.axis=0.7)
boxplot(beta7[, 3], xvar = "maPlate", yvar = "maA", outline=FALSE, las=2)


###################################################
### chunk number 23:  eval=FALSE
###################################################
## boxplot(beta7, main = "beta7 arrays", las=2)


###################################################
### chunk number 24: maBoxplotarrays
###################################################
par(mar=c(5, 3,3,3), cex.axis=0.7) #, cex.main=0.8)
boxplot(beta7, ylim=c(-4,4), main = "beta7 arrays", outline=FALSE, las=2)


###################################################
### chunk number 25: maplot2col
###################################################
plot(beta7nbg[,2], lines.func=NULL, legend.func=NULL)
points(beta7nbg[,2], subset=abs(maM(beta7nbg)[,2]) > 2,
       col="red", pch=18)
points(beta7nbg[,2], subset=maControls(beta7nbg) == "Empty", col="blue", pch=18)


###################################################
### chunk number 26: beta7normDo
###################################################
beta7norm <- cache("beta7norm", maNorm(beta7, norm="p"))


###################################################
### chunk number 27: beta7normShow eval=FALSE
###################################################
## beta7norm <- maNorm(beta7, norm="p")


###################################################
### chunk number 28: boxplotscale
###################################################
beta7norm.scale <- maNormScale(beta7norm)


###################################################
### chunk number 29: twoStepSeparateChanel
###################################################
beta7norm@maW <- matrix(0,0,0)      ## Remove weights
beta7.p <- as(beta7norm, "MAList")  ## convert data to RGList
beta7.pq <- normalizeBetweenArrays(beta7.p, method="quantile")


###################################################
### chunk number 30: plotdensityP
###################################################
plotDensities(beta7.p)


###################################################
### chunk number 31: plotdensityPQ
###################################################
plotDensities(beta7.pq)


###################################################
### chunk number 32: vsn0
###################################################
library("vsn")


###################################################
### chunk number 33: vsn1 eval=FALSE
###################################################
## library("vsn")
## beta7.vsn <- normalizeBetweenArrays(as(beta7, "RGList"), method="vsn")


###################################################
### chunk number 34: vsnDo
###################################################
beta7.vsn <- cache("beta7.vsn", vsn(beta7))


###################################################
### chunk number 35: vsnShow eval=FALSE
###################################################
## beta7.vsn <- vsn(beta7)


###################################################
### chunk number 36: getExprs
###################################################
b7 <- exprs(beta7.vsn) 


###################################################
### chunk number 37: vsnundercover
###################################################
fn <- as.character(maInfo(maTargets(beta7))$FileNames)
colnames(b7) <- paste(rep(fn, each=2), c("green", "red"), sep="\n")
b7 <- b7[sample(nrow(b7), 4000), ]


###################################################
### chunk number 38: plotvsn
###################################################
upPan <- function(...){
  points(..., col="darkblue")
  abline(a=0,b=1,col="red")
}
lowPan <- function(x, y, ...){
  text(mean(par("usr")[1:2]), mean(par("usr")[3:4]),signif(cor(x, y),2),cex=2)
}
pairs(b7[, 1:6], pch=".", lower.panel = lowPan, upper.panel=upPan)


###################################################
### chunk number 39: 
###################################################
library("beta7")


###################################################
### chunk number 40: 
###################################################
library("arrayQuality")


###################################################
### chunk number 41:  eval=FALSE
###################################################
## TargetInfo <- read.marrayInfo("TargetBeta7.txt")


###################################################
### chunk number 42:  eval=FALSE
###################################################
## mraw <- read.GenePix(targets = TargetInfo)


###################################################
### chunk number 43:  eval=FALSE
###################################################
## maQualityPlots(mraw)


###################################################
### chunk number 44:  eval=FALSE
###################################################
## normdata <- maNorm(mraw)


###################################################
### chunk number 45:  eval=FALSE
###################################################
## write.marray(normdata)


###################################################
### chunk number 46:  eval=FALSE
###################################################
## library("convert")
## mdata <- as(normdata, "exprSet")


###################################################
### chunk number 47:  eval=FALSE
###################################################
## LMres <- lmFit(normdata, design = c(1, -1, -1, 1, 1, -1), weights=NULL)


###################################################
### chunk number 48:  eval=FALSE
###################################################
## LMres <- eBayes(LMres)


###################################################
### chunk number 49:  eval=FALSE
###################################################
## restable <- topTable(LMres, number=10,resort.by="M")
## table2html(restable, disp="file")


