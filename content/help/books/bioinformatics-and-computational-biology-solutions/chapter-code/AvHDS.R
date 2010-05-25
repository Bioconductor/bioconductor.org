###################################################
### chunk number 1: setup
###################################################
library("RbcBook1")
library("Biobase")
library("multtest")
library("limma")
library("statmod")
library("genefilter")
library("annotate")
library("GO")
library("hgu95av2")
set.seed(4711)
quantilePlot <- function (px, py, quant=0.99, pch = ".", ...)
{
  n <- length(px)
  stopifnot(length(px)==n)
  rpx <- rank(px, na.last = FALSE)
  dm <- 0.05
  mids <- seq(dm, 1 - dm, by = dm)
  within <- function(x, x1, x2) {
    x >= x1 & x <= x2
  }
  quantilewind <- function(mp) {
     quantile(py[within(rpx/n, mp - dm, mp + dm)], quant, na.rm = TRUE)
  }
  plot(rpx, py, pch = pch, col="#999999", ...)
  q <- sapply(mids, quantilewind)
  lines(mids * n, q, col = "black", type = "b", pch = 19)
}


###################################################
### chunk number 2: dataALL
###################################################
library("ALL")
data(ALL)
pdat <- pData(ALL)
subset <- intersect(grep("^B", as.character(pdat$BT)),
          which(pdat$mol %in% c("BCR/ABL","NEG")))
eset <- ALL[,subset]


###################################################
### chunk number 3: prefiltering
###################################################
library("genefilter")
f1 <- pOverA(.25, log2(100))
f2 <- function(x)(IQR(x)>0.5)
ff <- filterfun(f1, f2)
selected <- genefilter(eset, ff)
sum(selected)
esetSub <- eset[selected,]


###################################################
### chunk number 4: multtestDo
###################################################
cl   <- as.numeric(esetSub$mol== "BCR/ABL")
## permutation test, using the Welch statistic
resT <- cache("multtest-resT", mt.maxT(exprs(esetSub), classlabel=cl, B=10000))
ord <- order(resT$index)  ## the original gene order
rawp <- resT$rawp[ord]    ## raw permutation $p$-values #names(resT)
names(rawp) <- geneNames(esetSub)


###################################################
### chunk number 5: multtest eval=FALSE
###################################################
## cl   <- as.numeric(esetSub$mol== "BCR/ABL")
## ## permutation test, using the Welch statistic
## resT <- mt.maxT(exprs(esetSub), classlabel=cl, B=10000)#1e+04)
## ord <- order(resT$index)  ## the original gene order
## rawp <- resT$rawp[ord]    ## raw permutation $p$-values #names(resT)
## names(rawp) <- geneNames(esetSub)


###################################################
### chunk number 6: histpraw
###################################################
hist(rawp, breaks=30, main="", xlab="", col="#B2DF8A")


###################################################
### chunk number 7: WY
###################################################
sum(resT$adjp<0.05)


###################################################
### chunk number 8: BH
###################################################
res <- mt.rawp2adjp(rawp, proc = "BH")


###################################################
### chunk number 9: 
###################################################
res <- mt.rawp2adjp(rawp, proc = "BH")
sum(res$adjp[,"BH"]<0.05)


###################################################
### chunk number 10: 
###################################################
rownames(res$adjp) <- names(rawp)


###################################################
### chunk number 11: calc.presel
###################################################
IQRs <- esApply(eset, 1, IQR)
intensityscore <- esApply(eset, 1, function(x) quantile(x, 0.75))
abs.t <- abs(mt.teststat(exprs(eset), classlabel=cl)) ## Welch statistic


###################################################
### chunk number 12: plotpresel
###################################################
par(mfrow=c(1,2))
quantilePlot(IQRs, abs.t, quant=0.95, xlab="rank(IQR)", ylab="abs.t")
quantilePlot(intensityscore, abs.t, quant=0.95,
        xlab="rank(intensity score)", ylab="abs.t")
par(mfrow=c(1,1))


###################################################
### chunk number 13: GOp1
###################################################
tykin <- unique(lookUp("GO:0004713", "hgu95av2", "GO2ALLPROBES"))
length(tykin)


###################################################
### chunk number 14: GO2
###################################################
gN <- geneNames(esetSub)
whSel <- match(tykin, gN)
whSel <- whSel[!is.na(whSel)]


###################################################
### chunk number 15: GO3
###################################################
resTGO <- cache("AvHDS-resTGO", mt.maxT(exprs(esetSub[whSel,]), classlabel=cl, B=10000))
ordGO <- order(resTGO$index) #$  the original gene order
adjpGO <- resTGO$adjp[ordGO] #$  raw permutation $p$-values #names(resT)
names(adjpGO) <- geneNames(esetSub[whSel,])
otherp <- resT$adjp[ord]     #$
names(otherp) <- geneNames(esetSub)
op2 <- otherp[names(adjpGO)]
ndiff <- sum(adjpGO < 0.05)


###################################################
### chunk number 16: GO4
###################################################
print("GO analysis")
sort(adjpGO)[1:ndiff]
print("All Genes")
sort(op2)[1:ndiff]


###################################################
### chunk number 17: preprocessKidney
###################################################
library("vsn")
library("kidpack")
data(qua)
data(spotanno)
data(hybanno)
data(cloneanno)
data(esetSpot)


###################################################
### chunk number 18: kidney.data
###################################################
pdat <- pData(esetSpot)
table(pdat$type)


###################################################
### chunk number 19: model
###################################################
design <- model.matrix(~ -1 + factor(pdat$type))
colnames(design) <- c("ccRCC", "chRCC", "pRCC")


###################################################
### chunk number 20: dupcorDo
###################################################
dupcor <- cache("AvHDS-dupcor", duplicateCorrelation(exprs(esetSpot),design=design,ndups=2,spacing=4224))
fit <- cache("AvHDS-fit", lmFit(esetSpot, design=design, ndups=2, spacing=4224, correlation=dupcor$cor))


###################################################
### chunk number 21: dupcorShow1 eval=FALSE
###################################################
## dupcor <- duplicateCorrelation(exprs(esetSpot),design=design,ndups=2,spacing=4224)
## fit <- lmFit(esetSpot, design=design, ndups=2, spacing=4224, correlation=dupcor$cor)


###################################################
### chunk number 22: dupcorShow2
###################################################
dupcor$cor


###################################################
### chunk number 23: fit.contrasts
###################################################
contrast.matrix <- makeContrasts(ccRCC-chRCC, ccRCC-pRCC, chRCC-pRCC, levels=design)
contrast.matrix
fit2 <- contrasts.fit(fit, contrast.matrix)


###################################################
### chunk number 24: tmod
###################################################
fit3 <- eBayes(fit2)


###################################################
### chunk number 25: topTable
###################################################
topTable(fit3, coef=3, n=8, adjust.method="fdr")


###################################################
### chunk number 26: decideTests
###################################################
clas <- decideTests(fit3, method="nestedF", adjust.method="fdr", p=0.05)
##get the numbers of significant genes per contrast:
colSums(abs(clas))


###################################################
### chunk number 27: compareDupCor
###################################################
nclones <- 4224
datAvDup <- (exprs(esetSpot)[1:nclones,] + exprs(esetSpot)[nclones+1:nclones,])/2
fitAvDup <- lmFit(datAvDup, design=design)
fit2AvDup <- contrasts.fit(fitAvDup, contrast.matrix)
fit3AvDup <- eBayes(fit2AvDup)


###################################################
### chunk number 28: AvDup
###################################################
plot(log10(fit3AvDup$p.value[,3]), log10(fit3$p.value[,3]), xlab="log10(p_avdup)", ylab="log10(p_dupcor)", xlim=c(-42,0), ylim=c(-42,0), pch=19, cex=0.3)
abline(c(0,1),col="#E31A1C")


###################################################
### chunk number 29: xxx
###################################################
gc()


###################################################
### chunk number 30: loadandnormalizeestrogendata
###################################################
library("estrogen")
library("limma")
library("hgu95av2cdf")

datadir <- system.file("extdata", package="estrogen")
targets <- readTargets("phenoData.txt",path=datadir,sep="")

covdesc <- list("present or absent","10 or 48 hours")
names(covdesc) <- names(targets)[-1]
pdata <- new("phenoData",pData=targets[,-1],varLabels=covdesc)
rownames(pData(pdata)) <- targets[,1]

gc()
esAB <- ReadAffy(filenames=file.path(datadir,targets$filename),phenoData=pdata)
esEset <- rma(esAB)


###################################################
### chunk number 31: 
###################################################
esEset
pData(esEset)



###################################################
### chunk number 32: outlierplots
###################################################
par(mfrow=c(1,2))
par(las=2)
for (i in c("728_at","33379_at")){
        expvals <- exprs(esEset)[i,]
        plot(expvals,axes=F,cex=1.5,
                xlab="Conditions",ylab="log 2 Expression Estimate")
        if(i=="728_at") points(7:8,expvals[7:8],pch=16,cex=1.5)
        if(i=="33379_at") points(1:2,expvals[1:2],pch=16,cex=1.5)
        axis(1,at=1:8,labels=c("et1","et2","Et1","Et2","eT1","eT2","ET1","ET2"))
        axis(2)
        title(i)
}
par(mfrow=c(1,1))


###################################################
### chunk number 33: outliers
###################################################
library("factDesign")
op1 <- outlierPair(exprs(esEset)["728_at",],INDEX=pData(esEset))
op1
madOutPair(exprs(esEset)["728_at",],op1[[3]])

op2 <- outlierPair(exprs(esEset)["33379_at",],INDEX=pData(esEset))
madOutPair(exprs(esEset)["33379_at",],op2[[3]])


###################################################
### chunk number 34: model
###################################################
pdat <- pData(esEset)
design <- model.matrix(~factor(estrogen)*factor(time.h),pdat)
colnames(design) <- c("Intercept","ES","T48","ES:T48")
fit <- lmFit(esEset,design)
fit$coefficients[1:3,]


###################################################
### chunk number 35: 
###################################################
contM <- cbind(es10=c(0,1,0,0),es48=c(0,1,0,1))
fitC <- contrasts.fit(fit,contM)
fitC <- eBayes(fitC)

esClas <- classifyTestsF(fitC, p=0.00001)
print(colSums(abs(esClas)))



###################################################
### chunk number 36: heatmap1
###################################################
both1048 <- which(esClas[,"es10"] != 0 & esClas[,"es48"] != 0)
heatmap(exprs(esEset)[both1048,],Colv=NA,col=cm.colors(256),labCol=targets[,"filename"])


###################################################
### chunk number 37: heatmap2
###################################################
library("stats")
only10 <- which(esClas[,"es10"] != 0 & esClas[,"es48"] == 0)
heatmap(exprs(esEset)[only10,],Colv=NA,col=cm.colors(256),labCol=targets[,"filename"],cexRow=1,cexCol=1)
title("Estrogen target at 10 hours only")


###################################################
### chunk number 38: heatmap3
###################################################
only48 <- which(esClas[,"es10"] == 0 & esClas[,"es48"] != 0)
heatmap(exprs(esEset)[only48,],Colv=NA,col=cm.colors(256),labCol=targets[,"filename"])
title("Estrogen target at 48 hours only")


