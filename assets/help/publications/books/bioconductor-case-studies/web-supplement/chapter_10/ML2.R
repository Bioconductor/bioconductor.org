###################################################
### chunk: loadAllLibsFirst
###################################################
library("Biobase")
library("genefilter")
library("cluster")
library("hgu95av2.db")
library("annotate")
library("MASS")
library("hopach")
library("kohonen")
library("RColorBrewer")


###################################################
### chunk: loadALLdata
###################################################
library("ALL")
data(ALL)
bcell = grep("^B", as.character(ALL$BT))
moltyp = which(as.character(ALL$mol.biol)
    %in% c("NEG", "BCR/ABL"))
ALL_bcrneg = ALL[, intersect(bcell, moltyp)]
ALL_bcrneg$mol.biol = factor(ALL_bcrneg$mol.biol)
ALLfilt_bcrneg = nsFilter(ALL_bcrneg, var.cutoff=0.75)$eset


###################################################
### chunk: restrictByGO
###################################################
GOTFfun = function(GOID) {
    x = hgu95av2GO2ALLPROBES[[GOID]]
    unique(x[ names(x) != "IEA"])
}

GOIDs = c("GO:0003700", "GO:0003702", "GO:0003709", 
    "GO:0016563", "GO:0016564")

TFs = unique(unlist(lapply(GOIDs, GOTFfun)))
inSel = match(TFs, featureNames(ALLfilt_bcrneg), nomatch=0)
es2 = ALLfilt_bcrneg[inSel,]


###################################################
### chunk: manh
###################################################
iqrs = esApply(es2, 1, IQR)
gvals = scale(t(exprs(es2)), rowMedians(es2), 
    iqrs[featureNames(es2)])
manDist = dist(gvals, method="manhattan")
hmcol = colorRampPalette(brewer.pal(10, "RdBu"))(256)
hmcol = rev(hmcol)
heatmap(as.matrix(manDist), sym=TRUE,  col=hmcol,
    distfun=function(x) as.dist(x))


###################################################
### chunk: mds
###################################################
cols = ifelse(es2$mol.biol == "BCR/ABL", "black", 
    "goldenrod")
sam1 = sammon(manDist, trace=FALSE)
plot(sam1$points, col=cols, xlab="Dimension 1", 
    ylab="Dimension 2")


###################################################
### chunk: samans
###################################################
sam2 = sammon(manDist, k=3,  trace=FALSE)


###################################################
### chunk: cmdscale
###################################################
cmd1 = cmdscale(manDist)


###################################################
### chunk: usett
###################################################
rtt = rowttests(ALLfilt_bcrneg, "mol.biol")
ordtt = order(rtt$p.value)
esTT = ALLfilt_bcrneg[ordtt[1:50],]
dTT = dist(t(exprs(esTT)), method="manhattan")
sTT = sammon(dTT, trace=FALSE)


###################################################
### chunk: silcheck
###################################################
mD = as.matrix(manDist)
silEst = silcheck(mD, diss=TRUE)
silEst
mD = as.hdist(mD)
mssEst = msscheck(mD)
mssEst

d2 = as.matrix(dist(t(gvals), method="man"))
silEstG = silcheck(d2, diss=TRUE)
silEstG
mssEstG = msscheck(as.hdist(d2))
mssEstG


###################################################
### chunk: silestSol
###################################################
dsol = as.matrix(dist(gvals), method="maximum")
silcheck(dsol, diss=TRUE)
msscheck(as.hdist(dsol))


###################################################
### chunk: hclust
###################################################
hc1 = hclust(manDist)
hc2 = hclust(manDist, method="single")
hc3 = hclust(manDist, method="ward")
hc4 = diana(manDist)


###################################################
### chunk: dendro
###################################################
par(mfrow=c(4,1))
plot(hc1, ann=FALSE) 
title(main="Complete Linkage", cex.main=2)
plot(hc2, ann=FALSE)
title(main="Single Linkage", cex.main=2)
plot(hc3, ann=FALSE)
title(main="Ward's Method", cex.main=2)
plot(hc4, ann=FALSE, which.plots=2)
title(main="Divisive Clustering", cex.main=2)
par(mfrow=c(1,1))


###################################################
### chunk: exCutree
###################################################
hc13 = cutree(hc1, k=3)
hc23 = cutree(hc2, k=3)
hc33 = cutree(hc3, k=3)
hc43 = cutree(hc4, k=3)


###################################################
### chunk: comp2dists
###################################################
table(hc13, hc33)


###################################################
### chunk: cophenetic1
###################################################
cph1 = cophenetic(hc1)
cor1 = cor(manDist, cph1)
cor1
plot(manDist, cph1, pch="|", col="blue")


###################################################
### chunk: cophensol
###################################################
cph2 = cophenetic(hc2)
cor2 = cor(manDist, cph2)
cor2

cph3 = cophenetic(hc3)
cor3 = cor(manDist, cph3)
cor3

cph4 = cophenetic(hc4)
cor4 = cor(manDist, cph4)
cor4


###################################################
### chunk: kmeans
###################################################
km2 = kmeans(gvals, centers=2, nstart=5)
kmx = kmeans(gvals, centers=2, nstart=25)


###################################################
### chunk: kmvals
###################################################
names(km2)
table(km2$cluster, kmx$cluster)


###################################################
### chunk: ALLfactors
###################################################

sapply(pData(es2), function(x) is.factor(x) || 
    is.logical(x) )



###################################################
### chunk: twowaytesting
###################################################

tt1 = table(es2$mdr, km2$cluster)
test1 = chisq.test(tt1)
test1$p.value


###################################################
### chunk: pamclustering
###################################################

pam2 = pam(manDist, k=2, diss=TRUE)
pam3 = pam(manDist, k=3, diss=TRUE)


###################################################
### chunk: comppamandkm
###################################################
all(names(km2$cluster) == names(pam2$clustering)) 

pam2km = table(km2$cluster, pam2$clustering)
pam2km


###################################################
### chunk: compKMPAM
###################################################
 km3 = kmeans(gvals, centers=3, nstart=25)


###################################################
### chunk: showcompKMPAM
###################################################
 table(km3$cluster, pam3$clustering)


###################################################
### chunk: tryouter
###################################################

 outSamekm3 = outer(km3$cluster, km3$cluster, 
                    FUN = function(x,y) x==y)

 outSamepam3 = outer(pam3$clustering, pam3$clustering, 
                     FUN = function(x,y) x==y)

 inSBoth = outSamekm3 & outSamepam3

 ##then we subtract 79, because an obs is in the same 
 ## cluster as itself this just means that the diagonal 
 ## is TRUE and divide by two
 sameBoth = (sum(inSBoth) - 79)/2

 ##not in the same one, in both are those FALSE entries
 notSBoth = (!outSamekm3) & (!outSamepam3)
 notSameBoth = sum(notSBoth)/2

 ##those that are different, are TRUE in one and FALSE 
 ## in the other or vice versa
 
 diffBoth = ((!outSamekm3) & outSamepam3) | 
    (outSamekm3 & (!outSamepam3))
 discordant = sum(diffBoth)/2



###################################################
### chunk: som
###################################################
set.seed(123)
s1 = som(gvals, grid=somgrid(4,4))
names(s1)
s2 = som(gvals, grid=somgrid(4,4), alpha=c(1,0.1), 
    rlen=1000)
s3 = som(gvals, grid=somgrid(4,4, topo="hexagonal"), 
    alpha=c(1,0.1), rlen=1000)
whGP = table(s3$unit.classif)
whGP


###################################################
### chunk: SOMsol1
###################################################
table(s1$unit.classif)
table(table(s2$unit.classif))
table(table(s3$unit.classif))


###################################################
### chunk: noname
###################################################
intOnes = s1$unit.classif == 13 | s1$unit.classif == 14
gvsub = gvals[intOnes,]
gvclasses = s1$unit.classif[intOnes]
sideC = ifelse(gvclasses==13, "yellow", "blue")
heatmap(t(gvsub), ColSideCol=sideC)



###################################################
### chunk: SOMBDR
###################################################
set.seed(777)
s4 = SOM(gvals, grid=somgrid(4,4, topo="hexagonal"))
SOMgp = knn1(s4$code, gvals, 1:nrow(s4$code))
table(SOMgp)


###################################################
### chunk: useMDS
###################################################
cD = dist(s4$code)
cD


###################################################
### chunk: dropCodes
###################################################
newCodes = s4$code[-c(2,5,6,10, 15, 9, 11, 13, 14),]

SOMgp2 = knn1(newCodes, gvals, 1:nrow(newCodes))
names(SOMgp2) = row.names(gvals)
table(SOMgp2)
cD2 = dist(newCodes)
cmdSOM = cmdscale(cD2)
#plot(cmdSOM)


###################################################
### chunk: km2sol
###################################################

km2sol = kmeans(gvals, centers=4, nstart=25)
table(km2sol$cluster, SOMgp2)



###################################################
### chunk: km3sol
###################################################
dropInds = SOMgp2 %in% c("1", "4", "5")
gvals2 = gvals[!dropInds,]
km3 = kmeans(gvals2, centers=4, nstart=50)
table(km3$cluster, SOMgp2[!dropInds])



###################################################
### chunk: hopachSamples
###################################################
  samp.hobj = hopach(gvals, dmat = manDist)
  samp.hobj$clust$k


###################################################
### chunk: hopachSamsizes
###################################################
samp.hobj$clust$sizes


###################################################
### chunk: hopachGenes
###################################################
gene.dist = distancematrix(t(gvals), d = "cosangle") 
gene.hobj = hopach(t(gvals), dmat = gene.dist) 
gene.hobj$clust$k


###################################################
### chunk: clustsize
###################################################
gene.hobj$clust$sizes


###################################################
### chunk: plotsil
###################################################
silpam2 = silhouette(pam2)


###################################################
### chunk: plotsilp2
###################################################
plot(silpam2, main="")


###################################################
### chunk: plotsil3
###################################################
silpam3 = silhouette(pam3)


###################################################
### chunk: plotsilp3
###################################################
plot(silpam3, main="")


###################################################
### chunk: findNegs
###################################################
silpam3[silpam3[,"sil_width"] < 0,]


###################################################
### chunk: negSol
###################################################
  ans = silpam2[silpam2[, "sil_width"] < 0, ]


###################################################
### chunk: make4grps
###################################################

 dcl4 = cutree(hc4, 4)
 table(dcl4)
 ## we presume the labels are in the order 
 ## given to the \indexTerm{clustering} algorithm

 sild4 =  silhouette(dcl4, manDist)


###################################################
### chunk: usettpc
###################################################
rtt = rowttests(ALLfilt_bcrneg, "mol.biol")
ordtt = order(rtt$p.value)
esTT = ALLfilt_bcrneg[ordtt[1:50],]


###################################################
### chunk: plotraw
###################################################
pairs(t(exprs(esTT)[1:5,]),
    col=ifelse(esTT$mol.biol=="NEG", "green", "blue"))


###################################################
### chunk: dopc
###################################################
pc = prcomp(t(exprs(esTT)))


###################################################
### chunk: plotpc
###################################################
pairs(pc$x[,1:5], col=ifelse(esTT$mol.biol=="NEG", 
    "green", "blue"))


###################################################
### chunk: plotbi
###################################################
biplot(pc)


###################################################
### chunk: dobt eval=FALSE
###################################################
## t.test(exprs(esTT)["33232_at",]~esTT$mol.biol)


###################################################
### chunk:  eval=FALSE
###################################################
## esTT.K = ALLfilt_bcrneg[ordtt[1:K],]


