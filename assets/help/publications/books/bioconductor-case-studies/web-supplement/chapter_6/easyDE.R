###################################################
### chunk: loaddat
###################################################
library("Biobase")
library("genefilter")
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
### chunk: rowSds
###################################################
library("genefilter")
sds = rowSds(exprs(ALL_bcrneg))
sh = shorth(sds)
sh


###################################################
### chunk: histsds
###################################################
hist(sds, breaks=50, col="mistyrose", xlab="standard deviation")
abline(v=sh, col="blue", lwd=3, lty=2)


###################################################
### chunk: unspecific
###################################################
ALLsfilt = ALL_bcrneg[sds>=sh, ]
dim(exprs(ALLsfilt))


###################################################
### chunk: ttest
###################################################
table(ALLsfilt$mol.biol)
tt = rowttests(ALLsfilt, "mol.biol")
names(tt)


###################################################
### chunk: histp
###################################################
hist(tt$p.value, breaks=50, col="mistyrose", xlab="p-value", 
     main="Retained")


###################################################
### chunk: histpunspec
###################################################
ALLsrest = ALL_bcrneg[sds<sh, ]
ttrest = rowttests(ALLsrest, "mol.biol")
hist(ttrest$p.value, breaks=50, col="lightblue", xlab="p-value", 
     main="Removed")


###################################################
### chunk: adjustp
###################################################
library("multtest")
mt = mt.rawp2adjp(tt$p.value, proc="BH")


###################################################
### chunk: list
###################################################
g = featureNames(ALLsfilt)[mt$index[1:10]]


###################################################
### chunk: anno
###################################################
library("hgu95av2.db")
links(hgu95av2SYMBOL[g])


###################################################
### chunk: matplot
###################################################
mb  = ALLsfilt$mol.biol
y   = exprs(ALLsfilt)[g[1],]
ord = order(mb)
plot(y[ord], pch=c(1,16)[mb[ord]], 
    col=c("black", "red")[mb[ord]], 
    main=g[1], ylab=expression(log[2]~intensity), 
    xlab="samples")


