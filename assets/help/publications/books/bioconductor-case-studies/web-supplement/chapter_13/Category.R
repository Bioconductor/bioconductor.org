###################################################
### chunk: ALL
###################################################
library("ALL")
data("ALL")


###################################################
### chunk: bcrabl
###################################################
bcell = grep("^B", as.character(ALL$BT))
moltyp = which(as.character(ALL$mol.biol) 
    %in% c("NEG", "BCR/ABL"))
ALL_bcrneg = ALL[, intersect(bcell, moltyp)]
ALL_bcrneg$mol.biol = factor(ALL_bcrneg$mol.biol)


###################################################
### chunk: nsFilter
###################################################
library("genefilter")
ALLfilt_bcrneg = nsFilter(ALL_bcrneg, var.cutoff=0.5)$eset


###################################################
### chunk: ans1
###################################################
table(ALLfilt_bcrneg$mol.biol)


###################################################
### chunk: KEGGimat
###################################################
library("GSEABase")
gsc = GeneSetCollection(ALLfilt_bcrneg, 
                  setType=KEGGCollection())
Am = incidence(gsc)
dim(Am)


###################################################
### chunk: EsetfromKEGG
###################################################
nsF = ALLfilt_bcrneg[colnames(Am),]


###################################################
### chunk: compttests
###################################################
rtt = rowttests(nsF, "mol.biol")
rttStat = rtt$statistic


###################################################
### chunk: reducetoInt
###################################################
selectedRows = (rowSums(Am)>10)
Am2 = Am[selectedRows, ]


###################################################
### chunk: pctests
###################################################
tA = as.vector(Am2 %*% rttStat)
tAadj = tA/sqrt(rowSums(Am2))
names(tA) = names(tAadj) = rownames(Am2)


###################################################
### chunk: Aufpassen
###################################################
stopifnot(sum(tAadj<(-8))==1L)


###################################################
### chunk: qqplot
###################################################
qqnorm(tAadj)


###################################################
### chunk: findSmPW
###################################################
library("KEGG.db")
smPW = tAadj[tAadj < (-8)]
pwName = KEGGPATHID2NAME[[names(smPW)]]
pwName


###################################################
### chunk: mnplot
###################################################
KEGGmnplot(names(smPW), nsF, "hgu95av2", nsF$"mol.biol",
  pch=16, col="darkblue")


###################################################
### chunk: heatmap
###################################################
sel = as.integer(nsF$mol.biol)
KEGG2heatmap(names(smPW), nsF, "hgu95av2", 
  col=colorRampPalette(c("white", "darkblue"))(256),
  ColSideColors=c("black", "white")[sel])


###################################################
### chunk: nlevelsnsFmolbiol
###################################################
stopifnot(nlevels(nsF$mol.biol)==2L)


###################################################
### chunk: ttperm
###################################################
library(Category)
set.seed(123)
NPERM = 1000

pvals = gseattperm(nsF, nsF$mol.biol, Am2, NPERM)

pvalCut = 0.025
lowC  = names(which(pvals[, 1]<=pvalCut))
highC = names(which(pvals[, 2]<=pvalCut))


###################################################
### chunk: makeSureWhatWeSayinTheTextIsRight eval=FALSE
###################################################
## stopifnot(identical(1L,length(lowC)))


###################################################
### chunk: fmt
###################################################
fmt = function(x) as.vector(mapply(paste, names(x), x, sep=": "))


###################################################
### chunk: lowC1 eval=FALSE
###################################################
## getPathNames(lowC)


###################################################
### chunk: lowC2
###################################################
fmt(getPathNames(lowC))


###################################################
### chunk: highC1 eval=FALSE
###################################################
## getPathNames(highC)


###################################################
### chunk: highC2
###################################################
fmt(getPathNames(highC))[1:5]


###################################################
### chunk: applypvals
###################################################
apply(pvals, 2, min)
rownames(pvals)[apply(pvals, 2, which.min)]


###################################################
### chunk: pvalcomp
###################################################
permpvs  = pmin(pvals[,1], pvals[,2])
pvsparam = pnorm(tAadj)
pvspara  = pmin(pvsparam, 1-pvsparam)


###################################################
### chunk: pvplot
###################################################
plot(permpvs, pvspara, xlab="Permutation p-values",
     ylab="Parametric p-values")


###################################################
### chunk: filterMissingMAP
###################################################
## depending on which annotation infrastructure we use
## hgu95av2MAP will either be an environment or an
## AnnDbBimap object
fnames = featureNames(ALLfilt_bcrneg)
if(is(hgu95av2MAP, "environment")){
    chrLocs = mget(fnames, hgu95av2MAP)
    mapping = names(chrLocs[sapply(chrLocs, 
        function(x) !all(is.na(x)))])
}else{
    mapping = toTable(hgu95av2MAP[fnames])$probe_id
}
psWithMAP = unique(mapping)
nsF2 = ALLfilt_bcrneg[psWithMAP, ]


###################################################
### chunk: findAmap
###################################################
EGtable = toTable(hgu95av2ENTREZID[featureNames(nsF2)])
entrezUniv = unique(EGtable$gene_id)
chrMat = MAPAmat("hgu95av2", univ=entrezUniv)
rSchr = rowSums(chrMat)


###################################################
### chunk: MAPAmat
###################################################
dim(chrMat)


###################################################
### chunk: reduceCols
###################################################
chrMat = chrMat[rowSums(chrMat) >= 5, ]
dim(chrMat)


###################################################
### chunk: reorderChrMat
###################################################
EGlist = mget(featureNames(nsF2),  hgu95av2ENTREZID)
EGIDs = sapply(EGlist, "[", 1)
idx = match(EGIDs, colnames(chrMat))
chrMat = chrMat[, idx]


###################################################
### chunk: overlap
###################################################
Ams = Am2[union(lowC, highC),]
Amx = Ams %*% t(Ams)
minS = outer(diag(Amx), diag(Amx), pmin)
overlapIndex = Amx/minS


###################################################
### chunk: overlapSetsPlot
###################################################
library("lattice")
myFun = function(x) {
  dd.row = as.dendrogram(hclust(dist(x)))
  row.ord = order.dendrogram(dd.row)
  
  dd.col = as.dendrogram(hclust(dist(t(x))))
  col.ord = order.dendrogram(dd.col)

  colnames(x) = sapply(getPathNames(colnames(x)), "[", 1L)

  levelplot(x[row.ord,col.ord],  scales = list(x = list(rot = 90)),
    xlab="", ylab="", main="overlapIndex",
    col.regions=colorRampPalette(c("white", "darkblue"))(256))
}
print(myFun(overlapIndex))


###################################################
### chunk: ijmin
###################################################
ijord = function(m, i) {
  m[lower.tri(m, diag=TRUE)] = NA
  o=order(m, decreasing=TRUE, na.last=TRUE)[i]
  cbind((o-1L)%/%nrow(m)+1L, (o-1L)%%nrow(m)+1L)
}  

nm = getPathNames(colnames(overlapIndex))
pathnames = as.vector(mapply(paste, names(nm), sapply(nm, "[", 1L), sep=": "))
k=ijord(overlapIndex, 1:2)


###################################################
### chunk: geneoverlap
###################################################
rowSums(Ams)[c("04510", "04512", "04514", "04940")]
Amx["04512", "04510"]
Amx["04940", "04514"]


###################################################
### chunk: lmfits
###################################################
P04512 = Ams["04512",]
P04510 = Ams["04510",]
lm1 = lm(rttStat ~ P04512)
summary(lm1)$coefficients
lm2 = lm(rttStat ~ P04510)
summary(lm2)$coefficients
lm3 = lm(rttStat ~ P04510+P04512)
summary(lm3)$coefficients


###################################################
### chunk: stopifnot
###################################################
stopifnot(
coefficients(summary(lm1))["P04512", "Pr(>|t|)"]<0.05,
coefficients(summary(lm2))["P04510", "Pr(>|t|)"]<0.05,
coefficients(summary(lm3))["P04512", "Pr(>|t|)"]>0.05,
coefficients(summary(lm3))["P04510", "Pr(>|t|)"]<0.05)


###################################################
### chunk: threegroups
###################################################
P04512.Only = ifelse(P04512 != 0 & P04510 == 0, 1, 0)
P04510.Only = ifelse(P04512 == 0 & P04510 != 0, 1, 0)
Both        = ifelse(P04512 != 0 & P04510 != 0, 1, 0)
lm4 = lm(rttStat ~ P04510.Only + P04512.Only + Both) 


###################################################
### chunk: tgshow eval=FALSE
###################################################
## summary(lm4)


###################################################
### chunk: pairpw2
###################################################
P04514 = Ams["04514",]
P04940 = Ams["04940",]
P04514.Only = ifelse(P04514 != 0 & P04940 == 0, 1, 0)
P04940.Only = ifelse(P04514 == 0 & P04940 != 0, 1, 0)
Both        = ifelse(P04514 != 0 & P04940 != 0, 1, 0)
lm5 = lm(rttStat ~ P04514.Only + P04940.Only + Both) 
summary(lm5)


