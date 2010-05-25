###################################################
### chunk number 1: setup
###################################################
library("RbcBook1")
library("xtable")
library("limma")
library("ALL")
library("class")  # these are for my methods table
library("cluster")
library("e1071")
library("gbm")
library("gpls")
library("ipred")
library("MASS")
library("nnet")
library("pamr")
library("randomForest")
library("rpart")
library("Biobase")
library("MLInterfaces")
library("rpart")
library("MLInterfaces")

fixlat <- function(x) gsub("_", "\\\\\\\\_", x)

setOldClass("nnet.formula")
setOldClass("rpart")
setOldClass("knnP")
setOldClass("svm")
example("planarPlot-methods")


###################################################
### chunk number 2: pickSamps
###################################################
data(ALL)
bio <- which( ALL$mol.biol %in% c("BCR/ABL", "NEG"))
isb <- grep("^B", as.character(ALL$BT))
kp <- intersect(bio,isb)
all2 <- ALL[,kp]
X <- t(exprs(all2))
Y <- as.character(all2$mol.biol)


###################################################
### chunk number 3: coverage
###################################################
.allObs <- function() 
sapply(.packages(), function(x) objects(pos=paste("package",x,sep=":")))

catalogFNames <- function(fnames) {
 ao <- .allObs()
 ao <- ao[!(names(ao) %in% c("MLInterfaces", "datasets"))]
 inpack <- function(n,ao) {
  tmp <- sapply(ao, function(x) length(grep(n,x))>0)
  names(ao)[tmp]
 }
 on <- sapply(fnames, function(x,y)inpack(x,y)[1], ao)
 split(fnames, on)
}

allml <- objects(pos="package:MLInterfaces")
allBi <- grep("B$", allml)
allB <- allml[allBi]
nc <- nchar(allB)
allB <- unlist( substring(allB, 1, nc-1) )
cata <- catalogFNames( allB )
slist <- lapply(cata, function(x) paste(x, collapse=", "))
smat <- data.frame(cbind(names(slist), as.character(unlist(slist))))
names(smat) <- c("Package", "Functions covered")


###################################################
### chunk number 4: setNg
###################################################
Ng <- 500


###################################################
### chunk number 5: grabALL
###################################################
library("MLInterfaces")
library("ALL")
data(ALL)



###################################################
### chunk number 6: filtALL
###################################################
bio <- which( ALL$mol.biol %in% c("BCR/ABL", "NEG"))
isb <- grep("^B", as.character(ALL$BT))
kp <- intersect(bio,isb)
all2 <- ALL[,kp]
tmp <- all2$mol.biol == "BCR/ABL"
tmp <- ifelse(tmp, "BCR/ABL", "NEG") 
pData(all2)$bcrabl <- factor(tmp)


###################################################
### chunk number 7: setNdiff
###################################################
Ndiff <- Ng


###################################################
### chunk number 8: limmaRun
###################################################
library("limma")
#f <- all2$bcrabl
des <- model.matrix(~all2$bcrabl)
fit <- lmFit(all2, des)
fit2 <- eBayes(fit)
Tdiff <- topTable(fit2, coef=2, Ndiff)
all2 <- all2[as.numeric(rownames(Tdiff)),]


###################################################
### chunk number 9: ldaArgs
###################################################
args(ldaB)


###################################################
### chunk number 10: ldaRun
###################################################
l1 <- ldaB( all2, "bcrabl", 1:40 )


###################################################
### chunk number 11: ldaConfuSave
###################################################
cm <- confuMat(l1)


###################################################
### chunk number 12: ldaConfuPrint
###################################################
confuMat(l1)


###################################################
### chunk number 13: knnRun
###################################################
k1 <- knnB( all2, "bcrabl", 1:40 )


###################################################
### chunk number 14: confKnnSave
###################################################
ck <- confuMat(k1)


###################################################
### chunk number 15: confKnnPrint
###################################################
confuMat(k1)


###################################################
### chunk number 16: 
###################################################
set.seed(44242)


###################################################
### chunk number 17: nnrun
###################################################
n1 <- nnetB( all2, "bcrabl", 1:40, size= 5, MaxNWts=10000 )


###################################################
### chunk number 18: nnconf
###################################################
confuMat(n1)


###################################################
### chunk number 19: nn2run
###################################################
n1b <- nnetB( all2, "bcrabl", 1:40, size= 5, MaxNWts=10000 )


###################################################
### chunk number 20: nn2conf
###################################################
confuMat(n1b)


###################################################
### chunk number 21: n2
###################################################
n2 <- nnetB( all2, "bcrabl", 1:40, size= 6, decay=.05, MaxNWts=10000 )


###################################################
### chunk number 22: n2c
###################################################
confuMat(n2)


###################################################
### chunk number 23: rungbm
###################################################
g1 <- gbmB( all2, "bcrabl", 1:40, n.minobsinnode=3, n.trees=1000 )
confuMat(g1)


###################################################
### chunk number 24: runrf
###################################################
rf1 <- randomForestB( all2, "bcrabl", 1:40, importance=TRUE )
confuMat(rf1)


###################################################
### chunk number 25: runsvm
###################################################
s1 <- svmB( all2, "bcrabl", 1:40 )
confuMat(s1)


###################################################
### chunk number 26: xvargs
###################################################
args(xval)


###################################################
### chunk number 27: xvdemo1
###################################################
xvloo <- xval( all2, "bcrabl", knnB, "LOO" )


###################################################
### chunk number 28: xvtables
###################################################
table(given=all2$bcrabl, predicted=xvloo)


###################################################
### chunk number 29: viplot2
###################################################
opar <- par(no.readonly=TRUE)
par(las=1, mar=c(6,6,6,6))
plot(getVarImp(g1), resolveenv=hgu95av2SYMBOL)
#nip <- summary(g1@RObject, cBars=20)
par(opar)


###################################################
### chunk number 30: pairhist
###################################################
gg <- ALL[ c(46,50), kp ]
gn <- geneNames(gg) # gg is from example(planarPlot-methods)
par(mfrow=c(2,2))
plot(density(exprs(gg[,gg$mol.biol=="NEG"])[1,]),xlab="NEG", ylab=gn[1], main=" ")
plot(density(exprs(gg[,gg$mol.biol=="NEG"])[2,]),xlab="NEG", ylab=gn[2], main=" ")
plot(density(exprs(gg[,gg$mol.biol=="BCR/ABL"])[1,]),xlab="BCR/ABL", ylab=gn[1], main=" ")
plot(density(exprs(gg[,gg$mol.biol=="BCR/ABL"])[2,]),xlab="BCR/ABL", ylab=gn[2], main=" ")


###################################################
### chunk number 31: eddApp
###################################################
set.seed(1234)
library("edd")


###################################################
### chunk number 32: eddApp2
###################################################
neg <- edd( gg[, gg$mol.biol=="NEG"] )
as.character(neg)
bcr <- edd( gg[, gg$mol.biol=="BCR/ABL"] )
as.character(bcr)


