###################################################
### chunk: loadData
###################################################
library("limma")
library("CCl4")
dataPath = system.file("extdata", package="CCl4")


###################################################
### chunk: ddp1 eval=FALSE
###################################################
## dir(dataPath)


###################################################
### chunk: readTargets
###################################################
adf = read.AnnotatedDataFrame("samplesInfo.txt", 
    path=dataPath)
adf


###################################################
### chunk: readData1
###################################################
targets = pData(adf)
targets$FileName = row.names(targets)


###################################################
### chunk: readData2
###################################################
RG = read.maimages(targets, path=dataPath, source="genepix")


###################################################
### chunk: RGgenes
###################################################
head(RG$genes)


###################################################
### chunk: imagePlots
###################################################
par(mfrow=c(5,1))
imageplot(log2(RG$Rb[,1]), RG$printer, low="white", 
    high="red")
imageplot(log2(RG$Gb[,1]), RG$printer, low="white", 
    high="green")
imageplot(rank(RG$Rb[,1]), RG$printer, low="white", 
    high="red")
imageplot(rank(RG$Gb[,1]), RG$printer, low="white", 
    high="green")
imageplot(rank(log(RG$R[,1])+log(RG$G[,1])), 
    RG$printer, low="white", high="blue")


###################################################
### chunk: calcma
###################################################
MA = normalizeWithinArrays(RG, method="none", 
    bc.method="none")


###################################################
### chunk: maplot
###################################################
library("geneplotter")
smoothScatter(MA$A[, 1], MA$M[, 1], xlab="A", ylab="M")
abline(h=0, col="red")


###################################################
### chunk: boxplot
###################################################
plotformula = log2(RG$G)~col(RG$G) 
boxplot(plotformula, ylim=c(5,9), outline=FALSE, 
    col="forestgreen", xlab="arrays", 
    ylab=expression(log[2]~G), main="boxplot")


###################################################
### chunk: multidensity
###################################################
multidensity(plotformula, xlim=c(5,9), 
    main="densities", xlab=expression(log[2]~G))


###################################################
### chunk: rin
###################################################
rin = with(MA$targets, ifelse(Cy5=="CCl4", RIN.Cy5, 
    RIN.Cy3))
rin
select = (rin == max(rin))
RGgood = RG[, select]
adfgood = adf[select, ]


###################################################
### chunk: justvsn eval=FALSE
###################################################
## library("vsn")
## ccl4 = justvsn(RGgood, backgroundsubtract=TRUE)


###################################################
### chunk: loadvsn
###################################################
library("vsn")


###################################################
### chunk: justvsn
###################################################
ccl4 = justvsn(RGgood, backgroundsubtract=TRUE)


###################################################
### chunk: meanSdPlot
###################################################
r = assayData(ccl4)$R
g = assayData(ccl4)$G
meanSdPlot(cbind(r, g))


###################################################
### chunk: gprextension
###################################################
rownames(pData(adfgood)) = sub("\\.gpr$", "", 
    rownames(pData(adfgood)))
pData(adfgood)


###################################################
### chunk: channel
###################################################
varMetadata(adfgood)$channel = factor(c("G", "R", "G", "R"),
    levels = c("G", "R", "_ALL_"))


###################################################
### chunk: stickiton
###################################################
phenoData(ccl4) = adfgood
validObject(ccl4)


###################################################
### chunk: calcAM
###################################################
ccl4AM = ccl4
assayData(ccl4AM) = assayDataNew(A=(r+g)/2, M=r-g)


###################################################
### chunk: vmdccl4am
###################################################
varMetadata(phenoData(ccl4AM))$channel[] = "_ALL_"
validObject(ccl4AM)


###################################################
### chunk: maplotsvsnShow eval=FALSE
###################################################
## smoothScatter(assayData(ccl4AM)$A[,2], 
##               assayData(ccl4AM)$M[,2])
## abline(h=0, col="red")


###################################################
### chunk: lm1
###################################################
design = modelMatrix(pData(ccl4AM), ref="DMSO")


###################################################
### chunk: lm2
###################################################
fit = lmFit(assayData(ccl4AM)$M, design)


###################################################
### chunk: eBayes
###################################################
fit = eBayes(fit)


###################################################
### chunk: lm2show
###################################################
class(fit)
names(fit)


###################################################
### chunk: histp
###################################################
hist(fit$p.value, 1000)


###################################################
### chunk: topgenes
###################################################
fit$genes = pData(featureData(ccl4AM))
topTable(fit, number=10, adjust="BH")


###################################################
### chunk: volcano
###################################################
plot(fit$coefficients, -log10(fit$p.value), pch=".")


###################################################
### chunk: write.fit eval=FALSE
###################################################
## write.fit(fit, file="fit.tab")


