###################################################
### chunk number 1: setup
###################################################
library("RbcBook1")
library("hgu133bcdf")


###################################################
### chunk number 2: loaddata
###################################################
library("affy")
library("ALLMLL")
data(MLL.B)
Data <- MLL.B[,c(2,1,3:5,14,6,13)] ##subset for some examples
sampleNames(Data) <- letters[1:8]


###################################################
### chunk number 3: rawImageDo
###################################################
bitmap("AffyQuality-image1a.png",height=3,width=3,pointsize=10,res=300)
par(mar=c(2.0,2.1,1.6,1.1),oma=c(1,1,0,0))
#image(Data[,1],transfo=function(x) x)
palette.gray <- c(rep(gray(0:10/10),times=seq(1,41,by=4)))
image(Data[,1],transfo=function(x) x, col=palette.gray)
dev.off()


###################################################
### chunk number 4: rawlogImageDo
###################################################
bitmap("AffyQuality-image1b.png",height=3,width=3,pointsize=10,res=300)
par(mar=c(2.0,2.1,1.6,1.1),oma=c(1,1,0,0))
#image(Data[,1]) ##the default transform is log
image(Data[,1],col=palette.gray)
dev.off()


###################################################
### chunk number 5: rawimageShow eval=FALSE
###################################################
## palette.gray <- c(rep(gray(0:10/10),times=seq(1,41,by=4)))
## image(Data[,1],transfo=function(x) x, col=palette.gray)
## image(Data[,1],col=palette.gray) ##the default transform is log


###################################################
### chunk number 6: boxplotShow eval=FALSE
###################################################
## library("RColorBrewer")
## cols <- brewer.pal(8,"Set1")  ##[-6] ##-6 takes out yellow
## boxplot(Data,col=cols)


###################################################
### chunk number 7: boxplot
###################################################
library("RColorBrewer")
cols <- brewer.pal(9,"Set1")[-6] ##-6 takes out yellow
boxplot(Data,col=cols)


###################################################
### chunk number 8: histShow eval=FALSE
###################################################
## hist(Data,col=cols,lty=1,xlab="Log (base 2) intensities")
## legend(12,1.0,letters[1:8],lty=1,col=cols)


###################################################
### chunk number 9: hist
###################################################
hist(Data,col=cols,lty=1,xlab="Log (base 2) intensities")
legend(12,1.0,letters[1:8],lty=1,col=cols)


###################################################
### chunk number 10: maplot1 eval=FALSE
###################################################
## #par(mar=c(2.4,2.2,1.6,1.1),oma=c(1,1,0,0))
## par(mfrow=c(2,4))
## MAplot(Data,cex=0.75)
## mtext("M",2,outer=TRUE)
## mtext("A",1,outer=TRUE)


###################################################
### chunk number 11: simpleAffyQC
###################################################
library("simpleaffy")
Data.qc <- qc(Data)


###################################################
### chunk number 12: avgbg
###################################################
avbg(Data.qc)


###################################################
### chunk number 13: scalefactor
###################################################
sfs(Data.qc)


###################################################
### chunk number 14: percentpresent
###################################################
percent.present(Data.qc)


###################################################
### chunk number 15: ratios35
###################################################
ratios(Data.qc)[,1:2]


###################################################
### chunk number 16: RNAdegDo
###################################################
library("AmpAffyExample")
data(AmpData)
## for illustrative purposes
sampleNames(AmpData) <- c("N1","N2","N3","A1","A2","A3") 
RNAdeg<-AffyRNAdeg(AmpData)
plotAffyRNAdeg(RNAdeg,col=c(2,2,2,3,3,3))


###################################################
### chunk number 17: RNAshow eval=FALSE
###################################################
## library("AmpAffyExample")
## data(AmpData)
## sampleNames(AmpData) <- c("N1","N2","N3","A1","A2","A3") #for illustrative purposes
## RNAdeg<-AffyRNAdeg(AmpData)
## plotAffyRNAdeg(RNAdeg,col=c(2,2,2,3,3,3))


###################################################
### chunk number 18: RNAdegsummmary
###################################################
summaryAffyRNAdeg(RNAdeg)


###################################################
### chunk number 19: fitPLM.AMPdata.cache
###################################################
library("affyPLM")
Pset1 <-cache("Pset1", fitPLM(AmpData))


###################################################
### chunk number 20: fitPLM.AMPdata eval=FALSE
###################################################
## library("affyPLM")
## Pset1 <- fitPLM(AmpData)
## show(Pset1)


###################################################
### chunk number 21: plmimageShow eval=FALSE
###################################################
## ##bitmap("AffyQuality-image2.png",height=5,width=6,pointsize=10,res=300)
## ##par(mar=c(2.0,2.1,1.6,1.1))
## par(mfrow=c(2,2))
## image(AmpData[,3])
## #mtext("A",side=3,line=0.2,adj=0,font=2)
## 
## image(Pset1,type="weights",which=3)
## #mtext("B",side=3,line=0.2,adj=0,font=2)
## 
## image(Pset1,type="resids",which=3)
## #mtext("C",side=3,line=0.2,adj=0,font=2)
## 
## image(Pset1,type="sign.resids",which=3)
## #mtext("D",side=3,line=0.2,adj=0,font=2)
## #dev.off()


###################################################
### chunk number 22: plmimageDo
###################################################
bitmap("AffyQuality-image2.png",height=5,width=6,pointsize=10,res=600)
par(mar=c(2.0,2.1,1.6,1.1))
par(mfrow=c(2,2))
image(AmpData[,3])
mtext("A",side=3,line=0.2,adj=0,font=2)

image(Pset1,type="weights",which=3)
mtext("B",side=3,line=0.2,adj=0,font=2)

image(Pset1,type="resids",which=3)
mtext("C",side=3,line=0.2,adj=0,font=2)

image(Pset1,type="sign.resids",which=3)
mtext("D",side=3,line=0.2,adj=0,font=2)
dev.off()


###################################################
### chunk number 23: fitPLM.MLLB.cache
###################################################
Pset2 <- cache("Pset2",fitPLM(MLL.B))


###################################################
### chunk number 24: fitPLM.MLLB eval=FALSE
###################################################
## Pset2 <- fitPLM(MLL.B)


###################################################
### chunk number 25: plmimage2Do
###################################################
bitmap("AffyQuality-image3.png",height=2.5,width=2.5,pointsize=8,res=300)
par(mar=c(2.0,2.1,1.6,1.1))
image(Pset2,which=2,type="resids")
dev.off()


###################################################
### chunk number 26: RLEShow eval=FALSE
###################################################
## Mbox(Pset2,ylim=c(-1,1),col=cols,names=NULL,main="RLE")


###################################################
### chunk number 27: RLE
###################################################
Mbox(Pset2,ylim=c(-1.7,1.1),col=cols,names=NULL,main="RLE",outline=FALSE)


###################################################
### chunk number 28: nuse
###################################################
boxplot(Pset2,ylim=c(0.95,1.5),col=cols,names=NULL, main="NUSE",outline=FALSE)


###################################################
### chunk number 29: compareAffyQCPLMQC-Do
###################################################
MLL.QC <- cache("AQ-MLL.QC", qc(MLL.B))
median.nuse <- cache("AQ-median.nuse", apply(nuse(Pset2,type="values"),2,median))


###################################################
### chunk number 30: compareAffyQCPLMQC
###################################################
par(mar=c(2.0,2.1,4.6,1.1))
par(mfrow=c(1,5)) 
boxplot(median.nuse,ylab="median NUSE",main="NUSE")
points(median.nuse[2],pch=20,cex=2,col="red")
boxplot(avbg(MLL.QC),ylab="Average Background",main="Avg bg")
points(avbg(MLL.QC)[2],pch=20,cex=2,col="red")
boxplot(sfs(MLL.QC),ylab="Scale Factor",main="SF")
points(sfs(MLL.QC)[2],pch=20,cex=2,col="red")
boxplot(percent.present((MLL.QC)),ylab="Percent Present",main="PP")
points(percent.present(MLL.QC)[2],cex=2,pch=20,col="red")
boxplot(2^ratios((MLL.QC))[,2],ylab="GAPDH 3'/5' ratio",main="3'/5'")
points(2^ratios((MLL.QC))[2,2],cex=2,pch=20,col="red")


