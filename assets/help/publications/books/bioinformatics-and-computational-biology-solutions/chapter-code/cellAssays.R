###################################################
### chunk number 1: setup
###################################################
library("RbcBook1")


###################################################
### chunk number 2: startup
###################################################
library("prada")
library("facsDorit")
library("geneplotter")


###################################################
### chunk number 3: checkVersions
###################################################
stopifnot(package.version("prada") >= package_version("1.1.13"))
stopifnot(package.version("facsDorit") >= package_version("1.1.0"))


###################################################
### chunk number 4: readViability
###################################################
viab <- read.table(system.file("extdata", "BoutrosKiger.tab",
           package="facsDorit"), header=TRUE, as.is=TRUE, sep="\t")
viab[1:2, ]


###################################################
### chunk number 5: checkViab
###################################################
fourcolumns <- c("Kc1", "Kc2", "S2R1", "S2R2")
stopifnot(all(fourcolumns %in% colnames(viab)))
stopifnot("Plate" %in% colnames(viab))


###################################################
### chunk number 6: 
###################################################
viab$Plate384 <- with(viab, (Plate-min(Plate)) %/% 4)


###################################################
### chunk number 7: bpprun
###################################################
rg  <- 1000*c(250, 2000) 
boxplot(log(Kc1) ~ Plate384, data=viab, outline=FALSE, col="#A6CEE3", 
        ylim=log(rg))


###################################################
### chunk number 8: scprun
###################################################
fac    <- factor(viab$Plate384)
colors <- rainbow(nlevels(fac))[as.integer(fac)]
perm   <- sample(nrow(viab))
plot(viab[perm, c("Kc1", "Kc2")], col=colors[perm], pch=".", 
     log="xy", xlim=rg, ylim=rg)


###################################################
### chunk number 9: normViability
###################################################
library("MASS")
expts  <- c("Kc1", "Kc2", "S2R1", "S2R2")
allMedian <- log(median(as.matrix(viab[, expts])))
for(ex in expts) {
  lmRes <- rlm(log(viab[[ex]]) ~ fac)
  viab[[paste("norm", ex, sep="")]] <- residuals(lmRes) + allMedian
}


###################################################
### chunk number 10: select
###################################################
score <- rowMeans(viab[, paste("normKc", 1:2, sep="")])
sel <- order(score)[1:3] 
viab[sel,1:7]


###################################################
### chunk number 11: bpprnr
###################################################
boxplot(normKc1 ~ Plate384, data=viab, outline=FALSE, col="#1F78B4", 
        ylim=log(rg))


###################################################
### chunk number 12: scprnr
###################################################
plot(viab[perm, c("normKc1", "normKc2")], pch=".", col=colors[perm],
     xlim=log(rg), ylim=log(rg))


###################################################
### chunk number 13: altModel eval=FALSE
###################################################
## lmRes <- rlm(log(luc) ~ expt * fac, data=rViab)


###################################################
### chunk number 14: readFCS1
###################################################
sampleDir <- system.file("extdata", "map", package="facsDorit")
B05 <- readFCS(file.path(sampleDir,  "060304MAPK controls.B05"))
B05


###################################################
### chunk number 15: readFCS2
###################################################
description(B05)[c(130,137,139)]


###################################################
### chunk number 16: readFCS3
###################################################
exprs(B05)[1:2,]


###################################################
### chunk number 17: readCytoSet
###################################################
mapk <- readCytoSet(path=sampleDir, phenoData="plateIndex.txt")
pData(mapk)[1:2,]


###################################################
### chunk number 18: readCytoSet3
###################################################
mapk[[1]]
exprs(mapk[[1]])[1:2,]


###################################################
### chunk number 19: FSCvsSSC
###################################################
x <- exprs(B05)[, c("FSC-H","SSC-H")]
plot(x, pch=20, col=densCols(x))


###################################################
### chunk number 20: fitNorm2
###################################################
nfit <- fitNorm2(x, scalefac=2)


###################################################
### chunk number 21: mahalanobis eval=FALSE
###################################################
## mh <- mahalanobis(x, center=nfit$mu, cov=nfit$S)
## plot(log(nfit$p), mh)


###################################################
### chunk number 22: plotNorm2-2
###################################################
plotNorm2(nfit, ellipse=TRUE)


###################################################
### chunk number 23: plotNorm2-3
###################################################
nfit3 <- fitNorm2(x, scalefac=3)
plotNorm2(nfit3, ellipse=TRUE)


###################################################
### chunk number 24: B05
###################################################
B05.sel <- B05[nfit$sel,]


###################################################
### chunk number 25: FL1vsFL4do1
###################################################
myPlot <- function(x) {
  ex <- exprs(x)[,c("FL1-H","FL4-H")]
  plot(ex, pch=20, col=densCols(ex),
       xlim=range(exprs(B05)[,"FL1-H"]),
       ylim=range(exprs(B05)[,"FL4-H"]))
}


###################################################
### chunk number 26: FL1vsFL4do2
###################################################
myPlot(B05)


###################################################
### chunk number 27: FL1vsFL4do3
###################################################
myPlot(B05.sel)


###################################################
### chunk number 28: FL1vsFL4show eval=FALSE
###################################################
## myPlot <- function(x) {
##   ex <- exprs(x)[,c("FL1-H","FL4-H")]
##   plot(ex, pch=20, col=densCols(ex))
## }
## myPlot(B05)
## myPlot(B05.sel)


###################################################
### chunk number 29: platePlot1
###################################################
nrCells <- csApply(mapk, nrow)
plotPlate(nrCells, nrow=8, ncol=12, main="Cell number",
          col=brewer.pal(9, "YlOrBr"), width=6.3)


###################################################
### chunk number 30: Rggobi eval=FALSE
###################################################
## library("Rggobi")
## x <- exprs(B05)
## gg <- ggobi(x)
## gg$setGlyphs(5,1, 1:nrow(x))
## gg$setColors(rep(9, nrow(x)))
## gg$scatterplot("FL1-H", "FL4-H")
## gg$setBrushColor(5)
## gg$setMode("Brush")


###################################################
### chunk number 31: preprocApo
###################################################
preprocess <- function(x) {
  for(i in 1:length(x)) {
    dat <- exprs(x[[i]])
    fn  <- fitNorm2(dat[, c("FSC-H", "SSC-H")], scalefac=1.5)
    x[[i]] <- dat[fn$sel,]
  }
  return(x)
}
apo <- readCytoSet(path=system.file("extdata", "apoptosis",
         package="facsDorit"),  phenoData="plateIndex.txt")
apoP <- preprocess(apo)


###################################################
### chunk number 32: apo-schematic
###################################################
plot(0,0, type="n", xlab="protein expression", ylab="caspase-3 activation",
xaxt="n", yaxt="n")
abline(v=0, h=0)
text(x=c(-0.5, 0.5, -0.5, 0.5), y=c(0.5, 0.5, -0.5, -0.5),
labels=c("untransfected\napoptotic cells\n(UA)", 
         "transfected\napoptotic cells\n(TA)",
         "untransfected\nnon-apoptotic cells\n(UN)",
         "transfected\nnon-apoptotic cells\n(TN)"), adj=0.5, cex=1.2)


###################################################
### chunk number 33: apo-mock
###################################################
calcthr <- function(x) {h <- hubers(x) ; h$mu+2.5*h$s}
mock <- exprs(apoP[[1]])[,c("FL1-H", "FL4-H")]
plot(mock, pch=".",
     xlab="protein expression", xlim=c(0, 1023), 
     ylab="caspase-3 activation", ylim=c(0, 1023))
thrYFP   <- calcthr(mock[,1])
thrCASP3 <- calcthr(mock[,2])
abline(v=thrYFP, h=thrCASP3, col="red", lty=2)


###################################################
### chunk number 34: cide
###################################################
cide <- exprs(apoP[[6]])[,c("FL1-H", "FL4-H")]


###################################################
### chunk number 35: FisherTest
###################################################
ct   <- thresholds(cide, xthr=thrYFP, ythr=thrCASP3)
fisher.test(ct)


###################################################
### chunk number 36: apo-bcl2
###################################################
plot(cide, pch=".",
     xlab="protein expression", ylab="caspase-3 activation",
     xlim=c(0, 1023), ylim=c(0, 1023))
abline(v=thrYFP, h=thrCASP3, col="red", lty=2)


###################################################
### chunk number 37: calcOdds
###################################################
calcOdds <- function(x){
   ct  <- thresholds(x[,c("FL1-H", "FL4-H")], xthr=thrYFP, ythr=thrCASP3)
   f <- fisher.test(ct)
   res <- -log10(f$estimate)
   return(ifelse((f$p.value > 0.01 | is.infinite(res)), 0, res))
}
odds <- csApply(apoP, calcOdds)


###################################################
### chunk number 38: platePlotApoOdds
###################################################
cols <-  brewer.pal(9, "Reds")[c(rep(1,4), 2:9)]
plotPlate(odds, nrow=8, ncol=12, main="log odds ratios", desc=c("act", ""),
          col=cols, width=6.3, na.action="omit", ind=pData(apo)$wellnr)


###################################################
### chunk number 39: dopreprocessMap
###################################################
mapkP <- preprocess(combineFrames(mapk, factor(pData(mapk)$clone)))


###################################################
### chunk number 40: plotMap
###################################################
par(mfrow=c(1,3))
groups <- c("MEK", "DSPP", "ERK")
names(groups) <- c("a)", "b)", "c)")
library("locfit")
for(i in seq(along=groups)) {
  dat <- exprs(mapkP[[groups[i]]])[, c("FL1-H", "FL4-H")]
  lcft <- locfit.robust(x=dat[, 1], y=dat[, 2], deg=1, alpha=1, maxk=512)
  sel <- sample(1:nrow(exprs(mapkP[[groups[i]]])), 2000)
  plot(dat[sel,], pch=20, col=densCols(dat[sel,]), main=paste(names(groups)[i], groups[i]))
  lines(lcft, col="red", lwd=2)
}


###################################################
### chunk number 41: plotMapSimp eval=FALSE
###################################################
## groups <- c("MEK", "DSPP", "ERK")
## for(i in groups) {
##   dat <- exprs(mapkP[[i]])[, c("FL1-H", "FL4-H")]
##   lcft <- locfit.robust(x=dat[, 1], y=dat[, 2], deg=1, alpha=1, maxk=512)
##   plot(dat, pch=20, col=densCols(dat), main=i)
##   lines(lcft, col="red", lwd=2)
## }


###################################################
### chunk number 42: renameGroups
###################################################
names(groups) <- groups


###################################################
### chunk number 43: calcMapSimp
###################################################
sapply(groups, function (i) {
  dat   <- exprs(mapkP[[i]])[, c("FL1-H", "FL4-H")]
  dlcft <- locfit.robust(x=dat[,1], y=dat[,2], deg=1, alpha=1, deriv=1,
                            maxk=512)
  pp    <- preplot(dlcft, newdata=600, band="local")
  c(delta=pp$fit, zscore=pp$fit/pp$se.fit)
})


