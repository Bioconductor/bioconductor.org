###################################################
### chunk number 1: setup
###################################################
library("RbcBook1")


###################################################
### chunk number 2: agrLoad
###################################################
library("cluster")
data(agriculture)


###################################################
### chunk number 3: agrAnal
###################################################
part<-pam(agriculture,k=2)
round(part$clusinfo,2)
hier<-diana(agriculture)


###################################################
### chunk number 4: agrPlot
###################################################
par(mfrow=c(1,2))
plot(part,which.plots=1,labels=3,col.clus=3,lwd=2,main="PAM")
plot(hier,which.plots=2,lwd=2,main="DIANA")


###################################################
### chunk number 5: agrMap
###################################################
heatmap(as.matrix(t(agriculture)),Rowv=NA,labRow=c("GNP","% in Agriculture"),cexRow=1,xlab="Country")


###################################################
### chunk number 6: MSSData
###################################################
mu<-c(1,2,5,6,14,15,18,19)
X<-matrix(rnorm(240*25,0,0.5),nrow=240,ncol=25)
step<-240/length(mu)
for(m in 1:length(mu))
 X[((m-1)*step+1):(m*step),]<-X[((m-1)*step+1):(m*step),]+mu[m]
D<-dist(X,method="euclidean")


###################################################
### chunk number 7: MSSClust
###################################################
library("hopach")
k.sil<-silcheck(X)[1]
k.mss<-msscheck(as.matrix(D))[1]
pam.sil<-pam(X,k.sil)
pam.mss<-pam(X,k.mss)


###################################################
### chunk number 8: MSSvSil
###################################################
image(1:240,1:240,as.matrix(D)[order(pam.sil$clust),order(pam.mss$clust)],col=topo.colors(80),xlab=paste("Silhouette (k=",k.sil,")",sep=""),ylab=paste("MSS (k=",k.mss,")",sep=""),main="PAM Clusters: Comparison of Two Criteria",sub="Ordered Euclidean Distance Matrix")
abline(v=cumsum(pam.sil$clusinfo[,1]),lty=2,lwd=2)
abline(h=cumsum(pam.mss$clusinfo[,1]),lty=3,lwd=2)


###################################################
### chunk number 9: kidpack
###################################################
library("kidpack")
data(eset, package="kidpack")
data(cloneanno, package="kidpack")


###################################################
### chunk number 10: geneSelection
###################################################
library("genefilter")
ff<-pOverA(0.5,log10(3))
subset<-genefilter(abs(exprs(eset)),filterfun(ff))
kidney<-exprs(eset)[subset,]
dim(kidney)
gene.names<-cloneanno[subset,"imageid"]
is.dup <- duplicated(gene.names)
gene.names[is.dup]<-paste(gene.names[is.dup],"B",sep="")
rownames(kidney)<-gene.names
colnames(kidney)<-paste("Sample",1:ncol(kidney),sep="")


###################################################
### chunk number 11: distance
###################################################
gene.dist<-distancematrix(kidney, d="cosangle")
dim(gene.dist)


###################################################
### chunk number 12: hopachGene
###################################################
gene.hobj<-hopach(kidney,dmat=gene.dist)


###################################################
### chunk number 13: hg2
###################################################
gene.hobj$clust$k
table(gene.hobj$clust$sizes)
gene.hobj$clust$labels[1:5]


###################################################
### chunk number 14: finalorder
###################################################
gn.ord <- gene.names[gene.hobj$fin$ord]
Bs<-grep("B", gn.ord)
spaces<-NULL
for(b in Bs){
name<-unlist(strsplit(gene.names[gene.hobj$fin$ord][b], "B"))
spaces<-c(spaces, diff(grep(name, gn.ord)))
}
table(spaces)


###################################################
### chunk number 15: PAM
###################################################
bestk<-silcheck(dissvector(gene.dist),diss=TRUE)[1]
pamobj<-pam(dissvector(gene.dist),k=bestk,diss=TRUE)
round(pamobj$clusinfo,2)


###################################################
### chunk number 16: bootstrap
###################################################
bobj<-boothopach(kidney,gene.hobj,B=100)


###################################################
### chunk number 17: bootplot eval=FALSE
###################################################
## bootplot(bobj, gene.hobj, ord="bootp", main="Renal Cell Cancer", showclusters=FALSE)


###################################################
### chunk number 18: hopachArray
###################################################
array.hobj<-hopach(t(kidney), d="euclid")


###################################################
### chunk number 19: har
###################################################
array.hobj$clust$k


###################################################
### chunk number 20: dplotArray
###################################################
tumortype<-unlist(strsplit(phenoData(eset)$type,"RCC"))
dplot(distancematrix(t(kidney), d="euclid"), array.hobj, labels=tumortype,main="Renal Cell Cancer: Array Clustering")


###################################################
### chunk number 21: output eval=FALSE
###################################################
## gene.acc<-cloneanno[subset,"AccNumber"]
## makeoutput(kidney, gene.hobj, bobj, file="kidney.out", gene.names=gene.acc)


###################################################
### chunk number 22: fuzzy eval=FALSE
###################################################
## gene.desc<-cloneanno[subset,"description"]
## boot2fuzzy(kidney, bobj, gene.hobj, array.hobj, file="kidneyFzy", gene.names=gene.desc)


###################################################
### chunk number 23: hierarchical eval=FALSE
###################################################
## hopach2tree(kidney, file="kidneyTree", hopach.genes=gene.hobj, hopach.arrays=array.hobj, dist.genes=gene.dist, gene.names=gene.desc)


