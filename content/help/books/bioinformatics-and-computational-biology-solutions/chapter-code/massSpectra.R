###################################################
### chunk number 1: setup
###################################################
library("RbcBook1")
library("Biobase")
stopifnot(package.version("PROcess") >= package_version("1.3.4"))


###################################################
### chunk number 2: raw1
###################################################
library("PROcess")
fdat <- system.file("Test", package="PROcess")
fs <- list.files(fdat, pattern="\\.*csv\\.*", full.names=TRUE)
f1 <- read.files(fs[1])
plot(f1, type="l", xlab="m/z")
title(basename(fs[1]))


###################################################
### chunk number 3: bslnoffSingle
###################################################
bseoff <- bslnoff(f1,method="loess",bw=0.1,xlab="m/z",plot=TRUE)
title(basename(fs[1]))


###################################################
### chunk number 4: isPeak1
###################################################
pkgobj <- isPeak(bseoff,span=81,sm.span=11,
 plot=TRUE, zerothrsh=2,area.w = 0.003, ratio = 0.2, main="a)")


###################################################
### chunk number 5: isPeak2
###################################################
specZoom(pkgobj, xlim=c(5000,10000), main="b)")


###################################################
### chunk number 6: plotCalidef
###################################################
amu.cali <- c(1084,1638,3496,5807,7034)
plotCali <- function(f, main, lab.cali) {
  x <- read.files(f)
  plot(x, main=main, ylim=c(0, max(x[,2])), type="n")
  abline(h=0, col="gray")
  abline(v=amu.cali, col="salmon")
  if(lab.cali)
    axis(3, at=amu.cali,labels=amu.cali,las=3,
         tick=FALSE, col="salmon", cex.axis=0.94)
  lines(x)
  return(invisible(x))
}
dir.cali <- system.file("calibration", package="PROcess")
files    <- dir(dir.cali, full.names=TRUE)
i <- seq(along=files)


###################################################
### chunk number 7: plotCalido
###################################################
opar <- par(mfrow=c(4, 2), mar=c(2,2,3,1))
mapply(plotCali, files, LETTERS[i], i<=2) 
par(opar)


###################################################
### chunk number 8: plotCalishow eval=FALSE
###################################################
## mapply(plotCali, files, LETTERS[i], i<=2) 


###################################################
### chunk number 9: batchoff
###################################################
Mcal <- cache("Mcal", rmBaseline(dir.cali))


###################################################
### chunk number 10: batchoffshow eval=FALSE
###################################################
## Mcal <- rmBaseline(dir.cali)


###################################################
### chunk number 11: normalize
###################################################
M.r <- renorm(Mcal, cutoff = 400)


###################################################
### chunk number 12: cutoffDo
###################################################
cts <- round(10^(seq(2, 3.5, length=15)))
sdsFirst <- cache("sdsFirst", sapply(cts, avesd, Ma=Mcal))


###################################################
### chunk number 13: cutoffShow eval=FALSE
###################################################
## cts <- round(10^(seq(2, 4, length=14)))
## sdsFirst <- sapply(cts, avesd, Ma=Mcal)


###################################################
### chunk number 14: cutoff
###################################################
plot(cts, sdsFirst, xlab="cutpoint", pch=21, bg="red", log="x",
        ylab="average sd")


###################################################
### chunk number 15: getPeaksBatch.do
###################################################
peakfile <- "calipeak.csv"
if(!file.exists(peakfile))
   getPeaks(M.r, peakfile, ratio=0.1)


###################################################
### chunk number 16: getPeaksBatch.show eval=FALSE
###################################################
## peakfile <- "calipeak.csv"
## getPeaks(M.r, peakfile, ratio=0.1)


###################################################
### chunk number 17: QC
###################################################
qualRes <- quality(M.r, peakfile, cutoff=400)


###################################################
### chunk number 18: protobiomarkers
###################################################
bmkfile <- "calibmk.csv"
bmk1 <- pk2bmkr(peakfile,M.r,bmkfile, p.fltr=.5)
mk1 <- round(as.numeric(gsub("M","",names(bmk1))))
mk1


###################################################
### chunk number 19: plotCali2def
###################################################
plotCali2 <- function(...) {
  x <- plotCali(...)
  lines(x[,1]*2, x[,2]+25, col="blue")
}


###################################################
### chunk number 20: plotCali2do
###################################################
opar <- par(mfrow=c(4, 2), mar=c(2,2,3,1))
mapply(plotCali2, files, LETTERS[i], i<=2) 
par(opar)


###################################################
### chunk number 21: plotCali2show eval=FALSE
###################################################
## mapply(plotCali2, files, LETTERS[i], i<=2) 


###################################################
### chunk number 22: brcabseoff
###################################################
library("ProData")
f45c <- system.file("f45c", package="ProData")
fs <- dir(f45c,full.names=TRUE)


###################################################
### chunk number 23: rmBaseline.do
###################################################
M1 <- cache("M1", rmBaseline(f45c))


###################################################
### chunk number 24: rmBaseline.show eval=FALSE
###################################################
## M1 <- rmBaseline(f45c)


###################################################
### chunk number 25: brcacutoff
###################################################
data(f45cbmk)
SpecGrp <- pData(f45cbmk)
fns <- colnames(M1)
gi <-  regexpr("i+[0-9]+", fns)
specName <- substr(fns, gi, gi+attr(gi, "match.length")-1)
mt  <- match(SpecGrp[,2], toupper(specName))
M2 <- M1[, mt]
colnames(M2) <- SpecGrp[,2]


###################################################
### chunk number 26: sdsSecond.do
###################################################
stopifnot(!any(is.na(mt)))
sdsSecond <- cache("sdsSecond", 
  sapply(cts, avesd, Ma=M2[,SpecGrp[,1]=="D"]))


###################################################
### chunk number 27: sdsSecond.show eval=FALSE
###################################################
## sdsSecond <- sapply(cts, avesd, Ma=M2[,SpecGrp[,1]=="D"])


###################################################
### chunk number 28: brcacutoffFig
###################################################
plot(cts, sdsSecond, xlab="cutpoint", pch=21, bg="red", log="x",
        ylab="average sd")


###################################################
### chunk number 29: brcarenorm.cache
###################################################
nM <- cache("nM", renorm(M2, cutoff=1000))


###################################################
### chunk number 30: brcarenorm eval=FALSE
###################################################
## nM <- renorm(M2, cutoff=1000)


###################################################
### chunk number 31: brcapeakcallquality.do
###################################################
peakfile <- "f45cpeak.csv"
if(!file.exists(peakfile))
  getPeaks(nM, peakfile, ratio=.1)
qu <- cache("brcapeakcallquality", 
  quality(nM, peakfile, cutoff=1000))


###################################################
### chunk number 32: brcapeakcallquality.show eval=FALSE
###################################################
## peakfile <- "f45cpeak.csv"
## getPeaks(nM, peakfile, ratio=.1)
## qu <- quality(nM, peakfile, cutoff=1000)


###################################################
### chunk number 33: anybad1
###################################################
bad <- qu[,1] < 0.4 & qu[,2] <0.1 & qu[,3] < 1/2
sum(bad)


###################################################
### chunk number 34: anybad2
###################################################
stopifnot(!any(bad))


###################################################
### chunk number 35: brcanormalize
###################################################
drop <- SpecGrp[,1]=="D"
nM1 <- nM[,!drop]


###################################################
### chunk number 36: gelmap
###################################################
Ma <- binning(nM1, breaks=300)
colnames(Ma) <- SpecGrp[!drop,1]
par(xpd=TRUE)
marks <- c(2666, 5055, 7560, 7934)
sel   <- as.numeric(rownames(Ma))<10000
gelmap(log(Ma+1)[sel,],
        at.mz=marks, at.col=c(25, 90, 135),
        col=gray(10:0/10))
segments(x0=1013*c(1,1), y0=c(55,119),
         x1=10000*c(1,1), y1=c(55,119), col="red")
arrows(marks,c(-5,-5), marks, c(1,1),
        length=.08, angle=20, col="red")


###################################################
### chunk number 37: brcapeakalignDo
###################################################
peakfile1 <- "f45cpeak1.csv"
if( !file.exists(peakfile1) )
   getPeaks(nM1, peakfile1, ratio=.1)
bmkfile <- "f45cbmk.csv"
bmk <- cache("bmk", 
   pk2bmkr(peakfile1,nM1,bmkfile, p.fltr=0.1, eps=.003))


###################################################
### chunk number 38: brcapeakalignShow eval=FALSE
###################################################
## peakfile1 <- "f45cpeak1.csv"
## getPeaks(nM1, peakfile1, ratio=.1)
## bmkfile <- "f45cbmk.csv"
## bmk <- pk2bmkr(peakfile1,nM1,bmkfile, p.fltr=.1, eps=.003)


###################################################
### chunk number 39: brcabmk
###################################################
mks <- round(as.numeric(gsub("M","",names(bmk))))
length(mks)
mks[1:12]


###################################################
### chunk number 40: is.multiple
###################################################
mults <- is.multiple(mks, k=2:5)
## length(mults)
mults[[1]]


