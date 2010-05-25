###################################################
### chunk number 1: setup
###################################################
library("RbcBook1")


###################################################
### chunk number 2: 
###################################################
library("affy")


###################################################
### chunk number 3: 
###################################################
library("SpikeInSubset")
data(spikein95)


###################################################
### chunk number 4: 
###################################################
pd <- data.frame(population=c(1,1,1,2,2,2),
replicate=c(1,2,3,1,2,3))
rownames(pd) <- sampleNames(spikein95)
vl <- list(population="1 is control, 2 is treatment",
           replicate="arbitrary numbering") 
phenoData(spikein95) <- new("phenoData",pData=pd,varLabels=vl)


###################################################
### chunk number 5: 
###################################################
eset <- rma(spikein95)


###################################################
### chunk number 6: 
###################################################
e <- exprs(eset)
dim(e)


###################################################
### chunk number 7: 
###################################################
pData(eset)


###################################################
### chunk number 8: 
###################################################
Index1 <- which(eset$population==1)
Index2 <- which(eset$population==2)


###################################################
### chunk number 9: rowM
###################################################
d <- rowMeans(e[,Index2])-rowMeans(e[,Index1])


###################################################
### chunk number 10: 
###################################################
a <- rowMeans(e)


###################################################
### chunk number 11: sumAbs
###################################################
sum(abs(d)>1)


###################################################
### chunk number 12: maplot
###################################################
plot(a,d,ylim=c(-1,1),main="A) MA-plot",pch=".")


###################################################
### chunk number 13: rowttests
###################################################
library("genefilter")
tt <- rowttests(e, factor(eset$population))


###################################################
### chunk number 14: volcano
###################################################
lod <- -log10(tt$p.value)
plot(d,lod,cex=.25,main="B) Volcano plot for $t$-test")
abline(h=2)


###################################################
### chunk number 15: volcano2
###################################################
o1 <- order(abs(d),decreasing=TRUE)[1:25]
o2 <- order(abs(tt$statistic),decreasing=TRUE)[1:25]
o <- union(o1,o2)
plot(d[-o],lod[-o],cex=.25,xlim=c(-1,1),ylim=range(lod),main="C) Close up of B)")
points(d[o1],lod[o1],pch=18,col="blue")
points(d[o2],lod[o2],pch=1,col="red")


###################################################
### chunk number 16: 
###################################################
library("limma")
design <- model.matrix(~factor(eset$population))
fit <- lmFit(eset, design)
ebayes <- eBayes(fit)


###################################################
### chunk number 17: volcano3
###################################################
lod <- -log10(ebayes$p.value[,2])
mtstat<- ebayes$t[,2]
o1 <- order(abs(d),decreasing=TRUE)[1:25]
o2 <- order(abs(mtstat),decreasing=TRUE)[1:25]
o <- union(o1,o2)
plot(d[-o],lod[-o],cex=.25,xlim=c(-1,1),ylim=c(0,4),main="D) Volcano plot for moderated $t$-test")
points(d[o1],lod[o1],pch=18,col="blue")
points(d[o2],lod[o2],pch=1,col="red")


###################################################
### chunk number 18: 
###################################################
sum(tt$p.value<=0.01)


###################################################
### chunk number 19: 
###################################################
data(spikein95)
spikedin <- colnames(pData(spikein95))
spikedIndex <- match(spikedin,geneNames(eset))


###################################################
### chunk number 20: 
###################################################
d.rank <- sort(rank(-abs(d))[spikedIndex])
t.rank <- sort(rank(-abs(tt$statistic))[spikedIndex])
mt.rank <- sort(rank(-abs(mtstat))[spikedIndex])
ranks <- cbind(mt.rank,d.rank,t.rank)
rownames(ranks) <- NULL
ranks


###################################################
### chunk number 21: 
###################################################
tab <- topTable(ebayes,coef=2,adjust="fdr",n=10)


###################################################
### chunk number 22: 
###################################################
tab[1:5,]


###################################################
### chunk number 23: 
###################################################
genenames <- as.character(tab$ID)


###################################################
### chunk number 24: annoEset
###################################################
annotation(eset)


###################################################
### chunk number 25: loadlibhgu95av2
###################################################
library("hgu95av2")


###################################################
### chunk number 26: loadlibxml
###################################################
library("XML")
library("annotate")
absts<-pm.getabst(genenames,"hgu95av2")


###################################################
### chunk number 27: absts12
###################################################
absts[[1]][[4]]


###################################################
### chunk number 28: pmtit2
###################################################
## wh 13.01.05: commented out pm.titles since it behaves weird. 
## Revert back to it when it is fixed?
##titl <- pm.titles(absts[2])
titl <- sapply(absts[[2]], articleTitle)
strwrap(titl, simplify=FALSE)


###################################################
### chunk number 29: 
###################################################
pro.res <- sapply(absts, function(x) pm.abstGrep("[Pp]rotein", x))
pro.res[[2]]


###################################################
### chunk number 30: 
###################################################
pmAbst2HTML(absts[[2]],filename="pm.html")


###################################################
### chunk number 31: 
###################################################
ll <- getLL(genenames,"hgu95av2")
sym <- getSYMBOL(genenames,"hgu95av2")


###################################################
### chunk number 32: output
###################################################
tab <- data.frame(sym,tab[,-1])
 htmlpage(ll,filename="report.html",title="HTML report",othernames=tab,table.head=c("Locus ID",colnames(tab)),table.center=TRUE)


###################################################
### chunk number 33: 
###################################################
library("KEGG")
library("GO")
library("annaffy")
atab <- aafTableAnn( genenames, "hgu95av2", aaf.handler() )
saveHTML(atab, file="report2.html")


