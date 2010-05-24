###################################################
### chunk: load0
###################################################
library("BiocCaseStudies")
options(width=56)


###################################################
### chunk: loadAllLibsFirst
###################################################
library("Biobase")
library("RColorBrewer")
library("bioDist")
library("genefilter")
library("class")
library("MLInterfaces")
library("hgu95av2.db")
library("annotate")
library("randomForest")


###################################################
### chunk: phenotypeSubset
###################################################
library("ALL")
data(ALL)
bcell = grep("^B", as.character(ALL$BT))
moltyp = which(as.character(ALL$mol.biol) 
    %in% c("NEG", "BCR/ABL"))
ALL_bcrneg = ALL[, intersect(bcell, moltyp)]
ALL_bcrneg$mol.biol = factor(ALL_bcrneg$mol.biol)


###################################################
### chunk: tables
###################################################
table(ALL_bcrneg$mol.biol)


###################################################
### chunk: nsFilter
###################################################
ALLfilt_bcrneg = nsFilter(ALL_bcrneg, var.cutoff=0.75)$eset


###################################################
### chunk: classBNsub
###################################################
class(ALLfilt_bcrneg)


###################################################
### chunk: rowIQRsFun
###################################################
rowIQRs = function(eSet) {
    numSamp = ncol(eSet)
    lowQ = rowQ(eSet, floor(0.25 * numSamp))
    upQ = rowQ(eSet, ceiling(0.75 * numSamp))
    upQ - lowQ
}


###################################################
### chunk: standardize
###################################################
standardize = function(x) (x - rowMedians(x)) / rowIQRs(x)
exprs(ALLfilt_bcrneg) = standardize(exprs(ALLfilt_bcrneg))


###################################################
### chunk: distEx
###################################################
eucD = dist(t(exprs(ALLfilt_bcrneg)))
eucM = as.matrix(eucD)
dim(eucM)


###################################################
### chunk: distHM
###################################################
library("RColorBrewer")
hmcol = colorRampPalette(brewer.pal(10, "RdBu"))(256)
hmcol = rev(hmcol)
heatmap(eucM, sym=TRUE, col=hmcol, distfun=as.dist)


###################################################
### chunk: distTau
###################################################
spD =  spearman.dist(ALLfilt_bcrneg)
## spD@Size
attr(spD, "Size")
spM = as.matrix(spD)


###################################################
### chunk: spFig
###################################################
heatmap(spM, sym=TRUE, col=hmcol, 
    distfun=function(x) as.dist(x))


###################################################
### chunk: nearestNeighbor
###################################################
closest.top("03002", eucM, 1)


###################################################
### chunk: cM
###################################################
cD = MIdist(ALLfilt_bcrneg)
cM = as.matrix(cD)
closest.top("03002", cM, 1)


###################################################
### chunk: Samples
###################################################
 Negs = which(ALLfilt_bcrneg$mol.biol == "NEG")
 Bcr = which(ALLfilt_bcrneg$mol.biol == "BCR/ABL")

 S1 = sample(Negs, 20, replace=FALSE)
 S2 = sample(Bcr, 20, replace = FALSE)
 TrainInd = c(S1, S2)
 TestInd = setdiff(1:79, TrainInd)


###################################################
### chunk: MLearn
###################################################
kans = MLearn( mol.biol ~ ., data=ALLfilt_bcrneg, 
    knnI(k=1,l=0), TrainInd)
confuMat(kans)
dldans = MLearn( mol.biol ~ ., ALLfilt_bcrneg, dldaI, 
    TrainInd)
confuMat(dldans)
ldaans = MLearn( mol.biol ~ ., ALLfilt_bcrneg, ldaI, 
    TrainInd)
confuMat(ldaans)


###################################################
### chunk: selbyttest
###################################################
Traintt = rowttests(ALLfilt_bcrneg[, TrainInd], "mol.biol")
ordTT = order(abs(Traintt$statistic), decreasing=TRUE) 
fNtt = featureNames(ALLfilt_bcrneg)[ordTT[1:50]]


###################################################
### chunk: selallbytt
###################################################
alltt = rowttests(ALLfilt_bcrneg, "mol.biol")
ordall = order(abs(alltt$statistic), decreasing=TRUE)
fNall = featureNames(ALLfilt_bcrneg)[ordall[1:50]]
intersect(fNall, fNtt)


###################################################
### chunk: knnFS
###################################################
BNf = ALLfilt_bcrneg[fNtt,]
knnf = MLearn( mol.biol ~ ., data=BNf, knnI(k=1,l=0), 
    TrainInd)
 confuMat(knnf)


###################################################
### chunk: repDLDA
###################################################

 dldtt = MLearn( mol.biol ~ ., BNf, dldaI, TrainInd)
 confuMat(dldtt)
 ldatt  = MLearn( mol.biol ~ ., BNf, ldaI, TrainInd)
 confuMat(ldatt)


###################################################
### chunk: dataRedux
###################################################
BNx = ALLfilt_bcrneg[1:1000,]


###################################################
### chunk: xval
###################################################
knnXval1 = MLearn(mol.biol~., data=BNx, knn.cvI(k=1, l=0),
    trainInd=1:ncol(BNx)) 


###################################################
### chunk: xval eval=FALSE
###################################################
## knnXval1 = MLearn(mol.biol~., data=BNx, knnI(k=1, l=0),
##     xvalSpec("LOO")) 


###################################################
### chunk: estErrRates
###################################################
knnCM = confuMat(knnXval1)
knnCM
#overall error rate
(knnCM[1,2] + knnCM[2,1])/sum(knnCM)
#class conditional error rates
knnCM[1,2]/sum(knnCM[1,])
knnCM[2,1]/sum(knnCM[2,])


###################################################
### chunk: xval-table
###################################################
confuMat(knnXval1)


###################################################
### chunk: realXval
###################################################
lk3f1 = MLearn(mol.biol~., data=BNx, knnI(k=1),
    xvalSpec("LOO", fsFun=fs.absT(50)))
confuMat(lk3f1)


###################################################
### chunk: solXval
###################################################
lk3f2 = MLearn(mol.biol~., data=BNx, knnI(k=1),
    xvalSpec("LOO", fsFun=fs.absT(5)))
confuMat(lk3f2)
table(unlist(fsHistory(lk3f2)))


###################################################
### chunk: findK2
###################################################
knnXval2 = MLearn(mol.biol~., data=BNx, knn.cvI(k=2, l=0), 
    trainInd=1:ncol(BNx))
confuMat(knnXval2)


###################################################
### chunk: findK3
###################################################
knnXval3 = MLearn(mol.biol~., data=BNx, knn.cvI(k=3, l=0), 
    trainInd=1:ncol(BNx))
confuMat(knnXval3)


###################################################
### chunk: findK5
###################################################
knnXval5 = MLearn(mol.biol~., data=BNx, knn.cvI(k=5, l=0),
    trainInd=1:ncol(BNx)) 
confuMat(knnXval5)


###################################################
### chunk: loadlibs
###################################################
library("randomForest")
set.seed(123)
rf1 = MLearn( mol.biol~., data=ALLfilt_bcrneg, 
    randomForestI, TrainInd, ntree=1000, mtry=55, 
    importance=TRUE) 


###################################################
### chunk: secondone
###################################################
rf2 = MLearn( mol.biol~., data=ALLfilt_bcrneg, 
    randomForestI, TrainInd, ntree=1000, mtry=10, 
    importance=TRUE) 


###################################################
### chunk: rfpred
###################################################
trainY = ALLfilt_bcrneg$mol.biol[TrainInd]
confuMat(rf1, "train")
confuMat(rf1, "test")


###################################################
### chunk: rfpred2
###################################################
confuMat(rf2, "train")
confuMat(rf2, "test")


###################################################
### chunk: errRates
###################################################
cf1 = confuMat(rf1)
overallErrM1 = (cf1[2,1] + cf1[1,2])/sum(cf1)
overallErrM1
perClass1 = c(cf1[1,2], cf1[2,1])/rowSums(cf1)
perClass1


###################################################
### chunk: errRates2
###################################################
cf2 = confuMat(rf2)
overallErrM2 = (cf2[2,1] + cf2[1,2])/sum(cf2)
overallErrM2
perClass2 = c(cf2[1,2], cf2[2,1])/rowSums(cf2)
perClass2


###################################################
### chunk: knnrecalled
###################################################
cfKNN = confuMat(knnf)
(cfKNN[1,2] + cfKNN[2,1])/sum(cfKNN)
#class conditional error rates
cfKNN[1,2]/sum(cfKNN[1,])
cfKNN[2,1]/sum(cfKNN[2,])



###################################################
### chunk: varImpPlot1
###################################################
opar = par(no.readonly=TRUE, mar=c(7,5,4,2))
par(las=2)
impV1 = getVarImp(rf1)
plot(impV1, n=15, plat="hgu95av2", toktype="SYMBOL")
par(opar)


###################################################
### chunk: varImpPlot2
###################################################
par(las=2, mar=c(7,5,4,2))
impV2 = getVarImp(rf2)
plot(impV2, n=15, plat="hgu95av2", toktype="SYMBOL")
par(opar)


###################################################
### chunk: impV
###################################################
impvars = function(x, which="MeanDecreaseAccuracy", k=10) {
    v1 = order(importance(x)[,which], decreasing=TRUE)
    importance(x)[v1[1:k],]
}
ivm1 = impvars(rf1@RObject, k=20)
ivm2 = impvars(rf2@RObject, k=20)
intersect(row.names(ivm1) , row.names(ivm2))


###################################################
### chunk: refRF
###################################################
rfRev = MLearn( mol.biol~., data=ALLfilt_bcrneg, 
    randomForestI, TestInd, ntree=2000, mtry=10, 
    importance=TRUE) 


###################################################
### chunk: rfRevshow eval=FALSE
###################################################
## rfRev


###################################################
### chunk: confuShow
###################################################
cfR = confuMat(rfRev)
cfR
overallErr = (cfR[2,1] + cfR[1,2])/sum(cfR)
overallErr
perClass = c(cfR[1,2], cfR[2,1])/rowSums(cfR)
perClass


###################################################
### chunk: rfALL
###################################################
rfAll = MLearn( mol.biol~., data=ALLfilt_bcrneg, 
    randomForestI, 1:79, ntree=1000, mtry=10, 
    importance=TRUE)


###################################################
### chunk: rfALLshow eval=FALSE
###################################################
## rfAll@RObject


###################################################
### chunk: threeGclass
###################################################
Bcell = grep("^B", ALL$BT)
ALLs = ALL[,Bcell]
types = c("BCR/ABL", "NEG", "ALL1/AF4")
threeG = ALLs$mol.biol %in% types
ALL3g = ALLs[,threeG]
qrange <- function(x)
    diff(quantile(x, c(0.1, 0.9)))
ALL3gf = nsFilter(ALL3g, var.cutoff=0.75, 
    var.func=qrange)$eset
ALL3gf$mol.biol = factor(ALL3gf$mol.biol)


###################################################
### chunk: mgsplit
###################################################
s1 = table(ALL3gf$mol.biol)
trainN = ceiling(s1/2)
sN = split(1:length(ALL3gf$mol.biol), ALL3gf$mol.biol)
trainInd = NULL
testInd = NULL
set.seed(777)
for(i in 1:3) {
    trI = sample(sN[[i]], trainN[[i]])
    teI = setdiff(sN[[i]], trI)
    trainInd = c(trainInd, trI)
    testInd = c(testInd, teI)
}
trainSet = ALL3gf[, trainInd]
testSet = ALL3gf[, testInd]


###################################################
### chunk: knnMV
###################################################

knn1MV = knn(t(exprs(trainSet)), t(exprs(testSet)), 
    trainSet$mol.biol)
tab1 = table(knn1MV, testSet$mol.biol)
tab1


