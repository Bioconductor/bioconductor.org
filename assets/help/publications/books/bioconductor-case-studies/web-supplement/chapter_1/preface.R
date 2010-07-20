###################################################
### chunk: loaddat
###################################################
library("Biobase")
library("ALL")
library("genefilter")
data("ALL")


###################################################
### chunk: bcells
###################################################
bcell = grep("^B", as.character(ALL$BT))


###################################################
### chunk: moltyp
###################################################
types = c("NEG", "BCR/ABL")
moltyp = which(as.character(ALL$mol.biol) %in% types)


###################################################
### chunk: subset
###################################################
ALL_bcrneg = ALL[, intersect(bcell, moltyp)]


###################################################
### chunk: cleanup
###################################################
ALL_bcrneg$mol.biol = factor(ALL_bcrneg$mol.biol)
ALL_bcrneg$BT = factor(ALL_bcrneg$BT)


###################################################
### chunk: nsfilter
###################################################
varCut = 0.5
filt_bcrneg = nsFilter(ALL_bcrneg, require.entrez=TRUE, 
    require.GOBP=TRUE, remove.dupEntrez=TRUE,
    var.func=IQR, var.cutoff=varCut, 
    feature.exclude="^AFFX")
filt_bcrneg$filter.log
ALLfilt_bcrneg = filt_bcrneg$eset


###################################################
### chunk: af4bcr
###################################################
types = c("ALL1/AF4", "BCR/ABL")
moltyp = which(ALL$mol.biol %in% types)
ALL_af4bcr = ALL[, intersect(bcell, moltyp)]
ALL_af4bcr$mol.biol = factor(ALL_af4bcr$mol.biol)
ALL_af4bcr$BT = factor(ALL_af4bcr$BT)
filt_af4bcr = nsFilter(ALL_af4bcr,require.entrez=TRUE, 
    require.GOBP=TRUE, remove.dupEntrez=TRUE,
    var.func=IQR, var.cutoff=varCut)
ALLfilt_af4bcr = filt_af4bcr$eset


