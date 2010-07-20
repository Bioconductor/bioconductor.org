###################################################
### chunk: Setup
###################################################
library("Biobase")
library("ALL")
library("hgu95av2.db")
library("annotate")
library("genefilter")
library("GOstats")
library("RColorBrewer")
library("xtable")
library("Rgraphviz")


###################################################
### chunk: subset
###################################################
data(ALL)
bcell = grep("^B", as.character(ALL$BT))
types = c("NEG", "BCR/ABL")
moltyp = which(as.character(ALL$mol.biol) %in% types)
ALL_bcrneg = ALL[, intersect(bcell, moltyp)]
ALL_bcrneg$mol.biol = factor(ALL_bcrneg$mol.biol)


###################################################
### chunk: sol1
###################################################
numSamp = length(ALL_bcrneg$mol.biol)
table(ALL_bcrneg$mol.biol)


###################################################
### chunk: sol2
###################################################
annotation(ALL_bcrneg)
length(featureNames(ALL_bcrneg))


###################################################
### chunk: nsfilter
###################################################
varCut = 0.5
filt_bcrneg = nsFilter(ALL_bcrneg, require.entrez=TRUE,
    require.GOBP=TRUE, remove.dupEntrez=TRUE,
    var.func=IQR, var.cutoff=varCut, 
    feature.exclude="^AFFX")
names(filt_bcrneg)
ALLfilt_bcrneg = filt_bcrneg$eset


###################################################
### chunk: nsFiltering-Y
###################################################
chrN = mget(featureNames(ALLfilt_bcrneg), envir=hgu95av2CHR)
onY = sapply(chrN, function(x) any(x == "Y"))
onY[is.na(onY)] = FALSE
ALLfilt_bcrneg = ALLfilt_bcrneg[!onY, ]


###################################################
### chunk: removeMissingSYMBOL
###################################################
hasSymbol = sapply(mget(featureNames(ALLfilt_bcrneg),
    envir=hgu95av2SYMBOL), function(x) 
    !(length(x) == 1 && is.na(x)))
ALLfilt_bcrneg = ALLfilt_bcrneg[hasSymbol, ]


###################################################
### chunk: defineGeneUniverse
###################################################
affyUniverse = featureNames(ALLfilt_bcrneg)
entrezUniverse = unlist(mget(affyUniverse, 
    hgu95av2ENTREZID))


###################################################
### chunk: altUniv
###################################################
## an alternate universe based on the entire chip
chipAffyUniverse = featureNames(ALLfilt_bcrneg)
chipEntrezUniverse = mget(chipAffyUniverse, hgu95av2ENTREZID)
chipEntrezUniverse = unique(unlist(chipEntrezUniverse))


###################################################
### chunk: parametric1
###################################################
ttestCutoff = 0.05
ttests = rowttests(ALLfilt_bcrneg, "mol.biol")

smPV = ttests$p.value < ttestCutoff

pvalFiltered = ALLfilt_bcrneg[smPV, ]
selectedEntrezIds = unlist(mget(featureNames(pvalFiltered),
    hgu95av2ENTREZID))


###################################################
### chunk: sumPV
###################################################
sumpv = sum(smPV)


###################################################
### chunk: withYourData1 eval=FALSE
###################################################
## ## if you are following along with your own data...
## entrezUniverse = unlist(mget(featureNames(yourData),
##     hgu95av2ENTREZID))
## pvalFiltered = yourData[hasSmallPvalue, ]
## selectedEntrezIds = unlist(mget(featureNames(pvalFiltered),
##     hgu95av2ENTREZID))


###################################################
### chunk: standardHyperGeo
###################################################
hgCutoff = 0.001
params = new("GOHyperGParams", geneIds=selectedEntrezIds,
    universeGeneIds=entrezUniverse, annotation="hgu95av2.db", 
    ontology="BP", pvalueCutoff=hgCutoff, conditional=FALSE,
    testDirection="over")


###################################################
### chunk: standardHGTEST
###################################################
hgOver = hyperGTest(params)


###################################################
### chunk: HGTestAns
###################################################
hgOver


###################################################
### chunk: summaryNames
###################################################
df = summary(hgOver)
names(df)


###################################################
### chunk: summary2
###################################################
df = summary(hgOver, pvalue=0.05, categorySize=350)
nrow(df)


###################################################
### chunk: achelp eval=FALSE
###################################################
## ? HyperGResult-accessors


###################################################
### chunk: htmlReportExample
###################################################
htmlReport(hgOver, file="ALL_hgo.html")


###################################################
### chunk: htmlBrowse eval=FALSE
###################################################
## browseURL("ALL_hgo.html")


###################################################
### chunk: termGraphs
###################################################
sigSub = termGraphs(hgOver)


###################################################
### chunk: sigSub
###################################################
numG = length(sigSub)
sizes = sapply(sigSub, numNodes)
sizes


###################################################
### chunk: cellComFig
###################################################
plotGOTermGraph(sigSub[[1]], hgOver, max.nchar=100)


###################################################
### chunk: condHyperGeo
###################################################
paramsCond = params
conditional(paramsCond) = TRUE


###################################################
### chunk: condhgRun
###################################################
hgCond = hyperGTest(paramsCond)


###################################################
### chunk: condhgResult
###################################################
hgCond


###################################################
### chunk: condhgSummary
###################################################
dfcond = summary(hgCond, categorySize=50)
## trim the term names for display purposes
trimTerm = function(x) {
    if (nchar(x) <= 20)
        x
    else
        paste(substr(x, 1, 20), "...", sep="")
}
dfcond$Term = sapply(dfcond$Term, trimTerm)
sizeOrd = order(dfcond$Size, decreasing=TRUE)
dfcond[sizeOrd, c("Count", "Size", "Term")]


###################################################
### chunk: compareResults
###################################################
stdIds = sigCategories(hgOver)
condIds = sigCategories(hgCond)
setdiff(stdIds, condIds)


###################################################
### chunk: cellcomcond
###################################################
terms = nodes(sigSub[[1]])
df = summary(hgCond, pvalue=0.5)[ , c("Term", "Pvalue")]
df$Pvalue = round(df$Pvalue, 3)
df$Term = sapply(df$Term, function(x) {
    if (nchar(x) <= 20) x
    else paste(substr(x, 1, 20), "...", sep="")
})
df[terms, ]


###################################################
### chunk: basicParams
###################################################
params = new("ChrMapHyperGParams",
    conditional=FALSE, testDirection="over",
    universeGeneIds=entrezUniverse,
    geneIds=selectedEntrezIds,
    annotation="hgu95av2", pvalueCutoff=0.05)
paramsCond = params
conditional(paramsCond) = TRUE


###################################################
### chunk: basicTest
###################################################
hgans = hyperGTest(params)
hgansCond = hyperGTest(paramsCond)


###################################################
### chunk: ChromResult
###################################################
summary(hgans, categorySize=10)


###################################################
### chunk: kegg1
###################################################
kparams = new("KEGGHyperGParams",
    geneIds=selectedEntrezIds, 
    universeGeneIds=entrezUniverse,
    annotation="hgu95av2",
    pvalueCutoff=0.05,
    testDirection="over")

kans = hyperGTest(kparams)


###################################################
### chunk: kegg2
###################################################
summary(kans)
kparamsUnder = kparams
testDirection(kparamsUnder) = "under"


###################################################
### chunk: keggUnder
###################################################
kansUnder = hyperGTest(kparamsUnder)


###################################################
### chunk: kegg3
###################################################
summary(kansUnder)


###################################################
### chunk: pfam1
###################################################
pparams = new("PFAMHyperGParams",
    geneIds=selectedEntrezIds,
    universeGeneIds=entrezUniverse,
    annotation="hgu95av2",
    pvalueCutoff=hgCutoff,
    testDirection="over")


pans = hyperGTest(pparams)


###################################################
### chunk: pfam2
###################################################
summary(pans)


