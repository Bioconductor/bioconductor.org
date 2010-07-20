###################################################
### chunk number 1: initialize
###################################################
library("RbcBook1")
library("Biobase")
stopifnot(package.version("RbcBook1") >= package_version("0.3.0"))


###################################################
### chunk number 2: affydist1
###################################################
library("affy")
library("SpikeInSubset")
library("RColorBrewer")
data(spikein95)
hist(spikein95,lwd=1,xlab=expression(log[2](intensity)),col=brewer.pal(8,"Dark2"))


###################################################
### chunk number 3: affydist2
###################################################
boxplot(spikein95,names=1:6,col="#B2DF8A",
 ylab=expression(log[2](intensity)),xlab="array")


###################################################
### chunk number 4: affydist3 eval=FALSE
###################################################
## library("affy")
## library("SpikeInSubset")
## data("spikein95")
## hist(spikein95)
## boxplot(spikein95)


###################################################
### chunk number 5: biasCalc
###################################################
data(spikein133)
Index <- which(probeNames(spikein133)%in%colnames(pData(spikein133)))
pms <- pm(spikein133)[Index,]
genenames <- probeNames(spikein133)[Index]
nominal <- t(sapply(probeNames(spikein133)[Index],function(i) pData(spikein133)[,i]))
x <- as.vector(log2(nominal))
y <- as.vector(log2(pms))
avg <- tapply(y,x,mean,na.rm=TRUE)


###################################################
### chunk number 6: bias
###################################################
plot(jitter(x),y,las=1,pch=".",
     ylab="Log (Base 2) PM Intensity",
     xlab="Nominal Log (Base 2) Concentration")
lines(as.numeric(names(avg)),avg,lwd=3,col="red")


###################################################
### chunk number 7: lognormalerrors
###################################################
qqnorm(y[x==0],pch=".")
qqline(y[x==0],col="red")


