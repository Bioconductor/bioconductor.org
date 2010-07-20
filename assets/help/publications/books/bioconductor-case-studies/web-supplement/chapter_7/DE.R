###################################################
### chunk: ALLdata
###################################################
library("ALL")
data("ALL")


###################################################
### chunk: bcell
###################################################
bcell = grep("^B", as.character(ALL$BT))


###################################################
### chunk: moltyp
###################################################
moltyp = which(as.character(ALL$mol.biol) 
    %in% c("NEG", "BCR/ABL"))


###################################################
### chunk: ALL_bcrneg
###################################################
ALL_bcrneg = ALL[, intersect(bcell, moltyp)]
ALL_bcrneg$mol.biol = factor(ALL_bcrneg$mol.biol)


###################################################
### chunk: meanSdPlot
###################################################
library("vsn")
meanSdPlot(ALL_bcrneg)


###################################################
### chunk: selectbySpread
###################################################
sds = esApply(ALL, 1, sd)
sel = (sds > quantile(sds, 0.8))
ALLset1 = ALL_bcrneg[sel,]


###################################################
### chunk: rowttests1
###################################################
library("genefilter")
tt = rowttests(ALLset1, "mol.biol")


###################################################
### chunk: rowttests2
###################################################
names(tt)


###################################################
### chunk: volcano
###################################################
plot(tt$dm, -log10(tt$p.value), pch=".", 
    xlab = expression(mean~log[2]~fold~change),
    ylab = expression(-log[10](p)))


###################################################
### chunk: sol2
###################################################
sum(tt$p.value<0.05 & abs(tt$dm)>0.5)


###################################################
### chunk: multtest
###################################################
library("multtest")
cl = as.numeric(ALLset1$mol.biol=="BCR/ABL") 
resT = mt.maxT(exprs(ALLset1), classlabel=cl, B=1000)
ord = order(resT$index)  ## the original gene order
rawp = resT$rawp[ord]    ## permutation p-values 


###################################################
### chunk: histpraw
###################################################
hist(rawp, breaks=50, col="#B2DF8A")


###################################################
### chunk: WY
###################################################
sum(resT$adjp<0.05)


###################################################
### chunk: BH
###################################################
res = mt.rawp2adjp(rawp, proc = "BH")
sum(res$adjp[,"BH"]<0.05)


###################################################
### chunk: limma
###################################################
library("limma")
design  = cbind(mean = 1, diff = cl)


###################################################
### chunk: lmFit
###################################################
fit = lmFit(exprs(ALLset1), design)
fit = eBayes(fit)


###################################################
### chunk: eBayes
###################################################
library("hgu95av2.db")
ALLset1Syms = unlist(mget(featureNames(ALLset1), 
    env = hgu95av2SYMBOL)) 
topTable(fit, coef = "diff", adjust.method = "fdr", 
    sort.by = "p", genelist = ALLset1Syms) 


###################################################
### chunk: plotlimma
###################################################
plot(-log10(tt$p.value), -log10(fit$p.value[, "diff"]),
    xlab = "-log10(p) from two-sample t-test",  
    ylab = "-log10(p) from moderated t-test (limma)", 
    pch=".")
abline(c(0, 1), col = "red")


###################################################
### chunk: ebsubs1
###################################################
subs = c(35, 65, 75, 1, 69, 71)
ALLset2 =  ALL_bcrneg[, subs]
table(ALLset2$mol.biol)


###################################################
### chunk: ebsubs2
###################################################
tt2  = rowttests(ALLset2, "mol.biol")
fit2 = eBayes(lmFit(exprs(ALLset2), design=design[subs, ]))


###################################################
### chunk: ebsubsplot
###################################################
plot(-log10(tt2$p.value), -log10(fit2$p.value[, "diff"]),
    xlab = "-log10(p) from two-sample t-test",  
    ylab = "-log10(p) from moderated t-test (limma)", 
    pch=".")
abline(c(0, 1), col = "red")


###################################################
### chunk: selectWeirdgene
###################################################
g = which(tt2$p.value < 1e-4 &  
    fit2$p.value[, "diff"] > 0.02)


###################################################
### chunk: g
###################################################
sel = (ALLset2$mol.bio == "BCR/ABL")+1
col = c("black", "red")[sel]
pch = c(1,16)[sel]
plot(exprs(ALLset2)[g,], pch=pch, col=col,
    ylab="expression")


###################################################
### chunk: ROC1
###################################################
rocs = rowpAUCs(ALLset1, "mol.biol", p=0.2)


###################################################
### chunk: pAUC
###################################################
j = which.max(area(rocs))
plot(rocs[j], main = featureNames(ALLset1)[j])


###################################################
### chunk: exprvals
###################################################
mtyp = ALLset1$mol.biol
sel = rep(1:2, each=rev(table(mtyp)))
plot(exprs(ALLset1)[j, order(mtyp)], pch=c(1,15)[sel], 
    col=c("black", "red")[sel], 
    main=featureNames(ALLset1)[j], 
    ylab=expression(log[2]~expression~level))
legend("bottomleft", col=c("black", "red"), 
    pch=c(1,15), levels(mtyp), bty="n")


###################################################
### chunk: rocvst1
###################################################
nrsel.ttest = function(x, pthresh=0.05) {
    pval = rowttests(x, "mol.biol")$p.value
    return(sum(pval < pthresh))
}


###################################################
### chunk: rocvst2
###################################################
nrsel.pAUC = function(x, pAUCthresh=2.5e-2) {
     pAUC = area(rowpAUCs(x, fac="mol.biol", p=0.1))
     return(sum(pAUC > pAUCthresh))
}


###################################################
### chunk: rocvst3
###################################################
library(BiocCaseStudies)
x = ALLset1[sample(nrow(ALLset1), 1000), ]


###################################################
### chunk: rocvstplot
###################################################
par(mfrow = c(1, 2))
resample(x, "nrsel.ttest")
resample(x, "nrsel.pAUC")


###################################################
### chunk: rocvst4 eval=FALSE
###################################################
## resample(x, "nrsel.ttest")
## resample(x, "nrsel.pAUC")


