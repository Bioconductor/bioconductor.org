###################################################
### chunk number 1: setup
###################################################
library("RbcBook1")


###################################################
### chunk number 2: readaffy eval=FALSE
###################################################
## library("affy")
## Data <- ReadAffy()


###################################################
### chunk number 3: loading
###################################################
library("affydata")
data(Dilution)


###################################################
### chunk number 4: pm
###################################################
pm(Dilution,"1001_at")[1:3,]


###################################################
### chunk number 5: matplot-show eval=FALSE
###################################################
## matplot(pm(Dilution,"1001_at"),type="l", xlab="Probe No.",
##         ylab="PM Probe intensity")
## matplot(t(pm(Dilution,"1001_at")),type="l",xlab="Array No.",
##         ylab="PM Probe intensity")


###################################################
### chunk number 6: matplot1
###################################################
matplot(pm(Dilution,"1001_at"),type="l",ylab="PM Probe intensity",xlab="Probe No.")


###################################################
### chunk number 7: matplot2
###################################################
matplot(t(pm(Dilution,"1001_at")),type="l",ylab="PM Probe intensity",xlab="Array No.",xaxt="n")
axis(1,1:4,1:4)


###################################################
### chunk number 8: pdata
###################################################
pData(Dilution)


###################################################
### chunk number 9: checkCNS
###################################################
stopifnot("sn19" %in% colnames(pData(Dilution)), all(Dilution$sn19==0))


###################################################
### chunk number 10: rma.bg.cache
###################################################
Dilution.bg.rma <- cache("Dilution.bg.rma",bg.correct(Dilution,method="rma"))


###################################################
### chunk number 11: rma.bg.show eval=FALSE
###################################################
## Dilution.bg.rma <- bg.correct(Dilution, method="rma")


###################################################
### chunk number 12: mas5.bg.cache
###################################################
Dilution.bg.mas <- cache("Dilution.bg.mas",bg.correct(Dilution,method="mas"))


###################################################
### chunk number 13: mas5.bg.show eval=FALSE
###################################################
## Dilution.bg.mas <- bg.correct(Dilution,method="mas")


###################################################
### chunk number 14: mas5.norm
###################################################
Dilution.norm.scale <- normalize(Dilution,
                                 method="constant")


###################################################
### chunk number 15: invariantset.norm.cache
###################################################
Dilution.norm.nl <- cache("Dilution.norm.nonlinear",normalize(Dilution, method = "invariantset"))


###################################################
### chunk number 16: invariantset.norm.show eval=FALSE
###################################################
## Dilution.norm.nl <- normalize(Dilution, method="invariantset")


###################################################
### chunk number 17: quantile.norm
###################################################
Dilution.norm.quantile <- normalize(Dilution,
                                    method="quantiles")


###################################################
### chunk number 18: loess.norm.cache
###################################################
Dilution.norm.loess <- cache("Dilution.norm.loess",normalize(Dilution,"loess"))


###################################################
### chunk number 19: loess.norm.show eval=FALSE
###################################################
## Dilution.norm.loess <- normalize(Dilution,
##                                  method="loess")


###################################################
### chunk number 20: contrast.norm.cache
###################################################
Dilution.norm.contrast <- cache("Dilution.norm.contrast",normalize(Dilution,"contrast"))


###################################################
### chunk number 21: contrast.norm.show eval=FALSE
###################################################
## Dilution.norm.contrast <- normalize(Dilution,
##                                     method="contrast")


###################################################
### chunk number 22: vsn.norm
###################################################
library("vsn")
Dil.vsn <- normalize(Dilution, method = "vsn")


###################################################
### chunk number 23: normalizemethods
###################################################
normalize.methods(Dilution)


###################################################
### chunk number 24: expresso1.cache
###################################################
eset <- cache("expresso1", expresso(Dilution,
                 bgcorrect.method="rma",
                 normalize.method="constant",
                 pmcorrect.method="pmonly",
                 summary.method="avgdiff"))


###################################################
### chunk number 25: expresso1.show eval=FALSE
###################################################
## eset <- expresso(Dilution,
##                  bgcorrect.method="rma",
##                  normalize.method="constant",
##                  pmcorrect.method="pmonly",
##                  summary.method="avgdiff")


###################################################
### chunk number 26: expresso2.cache
###################################################
eset <- cache("expresso2", expresso(Dilution,
                 normalize.method="invariantset",
                 bg.correct=FALSE,
                 pmcorrect.method="pmonly",
                 summary.method="liwong"))


###################################################
### chunk number 27: expresso2.show eval=FALSE
###################################################
## eset <- expresso(Dilution,
##                  normalize.method="invariantset",
##                  bg.correct=FALSE,
##                  pmcorrect.method="pmonly",
##                  summary.method="liwong")


###################################################
### chunk number 28: mas5.cache
###################################################
eset <- cache("mas5", mas5(Dilution))


###################################################
### chunk number 29: mas5 eval=FALSE
###################################################
## eset <- mas5(Dilution)


###################################################
### chunk number 30: threestep
###################################################
library("affyPLM")
eset <- threestep(Dilution,
  background.method="IdealMM",
  normalize.method="quantile",
  summary.method="tukey.biweight")


###################################################
### chunk number 31: rma
###################################################
eset <- rma(Dilution)


###################################################
### chunk number 32: gcrma1.cache
###################################################
library("gcrma")
Dil.expr <- cache("gcrma1", gcrma(Dilution))


###################################################
### chunk number 33: gcrma1.show eval=FALSE
###################################################
## library("gcrma")
## Dil.expr <- gcrma(Dilution)


###################################################
### chunk number 34: gcrma2a.cache
###################################################
ai <- cache("gcrma2a", compute.affinities(cdfName(Dilution)))


###################################################
### chunk number 35: gcrma2b.cache
###################################################
Dil.expr<-cache("gcrma2b", gcrma(Dilution,affinity.info=ai))


###################################################
### chunk number 36: gcrma2.show eval=FALSE
###################################################
## ai <- compute.affinities(cdfName(Dilution))
## Dil.expr<-gcrma(Dilution,affinity.info=ai)


###################################################
### chunk number 37: gcrma3.cache
###################################################
Dil.expr2 <- cache("gcrma3",
   gcrma(Dilution,affinity.info=ai,type="affinities"))


###################################################
### chunk number 38: gcrma3.show eval=FALSE
###################################################
## Dil.expr2 <- gcrma(Dilution,affinity.info=ai,type="affinities")


###################################################
### chunk number 39: energyfiles eval=FALSE
###################################################
## library("affypdnn")
## energy.files <- list.files(system.file("exampleData",package = "affypdnn"),
##                            "^pdnn-energy-parameter")
## energyfile <- file.path(system.file("exampleData", package="affypdnn"),
##                         energy.files[1])
## ep <- read.table(energyfile, nrows=80, header=TRUE)
## Wg <- as.vector(ep[33:56, 2]) ## weights (specific)
## Wn <- as.vector(ep[57:80, 2]) ## weights (unspecific)


###################################################
### chunk number 40: expressopdnn eval=FALSE
###################################################
## hgu95av2.pdnn.params <- pdnn.params.chiptype(energyfile, probes.pack = "hgu95av2probe")
## attach(hgu95av2.pdnn.params)
## par.ct <- list(params.chiptype = hgu95av2.pdnn.params)
## eset <- expressopdnn(Dilution[, 1], findparams.param = par.ct)
## detach("hgu95av2.pdnn.params")


###################################################
### chunk number 41: affycompload
###################################################
library("affycomp")
data(dilution.phenodata)
data(spikein.phenodata)
data(hgu133a.spikein.phenodata)


###################################################
### chunk number 42: affycompdataload
###################################################
data(rma.assessment)
data(mas5.assessment)


###################################################
### chunk number 43: affycompplot1Show eval=FALSE
###################################################
## affycompPlot(rma.assessment$MA)


###################################################
### chunk number 44: affycompplot2Show eval=FALSE
###################################################
## affycompPlot(rma.assessment$Signal,mas5.assessment$Signal)


###################################################
### chunk number 45: affycompplot1
###################################################
affycompPlot(rma.assessment$MA)


###################################################
### chunk number 46: affycompplot2
###################################################
affycompPlot(rma.assessment$Signal,mas5.assessment$Signal)


