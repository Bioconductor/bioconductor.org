###################################################
### chunk: load0
###################################################
library("BiocCaseStudies")


###################################################
### chunk: ALLdata
###################################################
library("ALL")
data("ALL")
types = c("ALL1/AF4", "BCR/ABL")
bcell = grep("^B", as.character(ALL$BT))
ALL_af4bcr = ALL[, intersect(bcell, 
    which(ALL$mol.biol %in% types))]
ALL_af4bcr$mol.biol = factor(ALL_af4bcr$mol.biol)


###################################################
### chunk: groupsize
###################################################
table(ALL_af4bcr$mol.biol)


###################################################
### chunk: nsfilter
###################################################
qrange <- function(x)
    diff(quantile(x, c(0.1, 0.9)))
library("genefilter")
filt_af4bcr = nsFilter(ALL_af4bcr, require.entrez=TRUE, 
    require.GOBP=TRUE, var.func=qrange, var.cutoff=0.5)
ALLfilt_af4bcr = filt_af4bcr$eset


###################################################
### chunk: loadlibs
###################################################
library("Biobase")
library("annotate")
library("hgu95av2.db")


###################################################
### chunk: rt
###################################################
rt = rowttests(ALLfilt_af4bcr, "mol.biol")
names(rt)


###################################################
### chunk: figrtsol1
###################################################
hist(rt$statistic, breaks=100, col="skyblue", 
     main="", xlab="t-statistic")


###################################################
### chunk: figrtsol2
###################################################
hist(rt$p.value, breaks=100, col="mistyrose",
     main="", xlab="p-value")


###################################################
### chunk: ALLsub
###################################################
sel = order(rt$p.value)[1:400]
ALLsub = ALLfilt_af4bcr[sel,]


###################################################
### chunk: EGdup1
###################################################
EG    = as.character(hgu95av2ENTREZID[featureNames(ALL)])
EGsub = as.character(hgu95av2ENTREZID[featureNames(ALLsub)])


###################################################
### chunk: tabletable
###################################################
table(table(EG))
table(table(EGsub))


###################################################
### chunk: figprofile
###################################################
syms = as.character(hgu95av2SYMBOL[featureNames(ALLsub)])
whFeat = names(which(syms =="CD44"))
ordSamp = order(ALLsub$mol.biol)
CD44 = ALLsub[whFeat, ordSamp]
plot(as.vector(exprs(CD44)), main=whFeat,
    col=c("sienna", "tomato")[CD44$mol.biol], 
    pch=c(15, 16)[CD44$mol.biol], ylab="expression")


###################################################
### chunk: chrtab
###################################################
z = toTable(hgu95av2CHR[featureNames(ALLsub)])
chrtab = table(z$chromosome)
chrtab


###################################################
### chunk: figchrtab
###################################################
chridx = sub("X", "23", names(chrtab))
chridx = sub("Y", "24", chridx)
barplot(chrtab[order(as.integer(chridx))])


###################################################
### chunk: annaffy1
###################################################
library("annaffy")
anncols = aaf.handler(chip="hgu95av2.db")[c(1:3, 8:9, 11:13)]
anntable = aafTableAnn(featureNames(ALLsub), 
    "hgu95av2.db", anncols)
saveHTML(anntable, "ALLsub.html", 
    title="The Features in ALLsub")
localURL = file.path("file:/", getwd(), "ALLsub.html")


###################################################
### chunk: annaffy3
###################################################
browseURL(localURL)


###################################################
### chunk: figduplo0
###################################################
probeSetsPerGene = split(names(EG), EG)
j = probeSetsPerGene$"7013"
j


###################################################
### chunk: figduplo1
###################################################
plot(t(exprs(ALL_af4bcr)[j[c(1,7)], ]), asp=1, pch=16,
    col=ifelse(ALL_af4bcr$mol.biol=="ALL1/AF4", "black", 
    "grey"))


###################################################
### chunk: figduplo2
###################################################
library("lattice")
mat = exprs(ALL_af4bcr)[j,]
mat = mat - rowMedians(mat)
ro = order.dendrogram(as.dendrogram(hclust(dist(mat))))
co = order.dendrogram(as.dendrogram(hclust(dist(t(mat)))))
at = seq(-1, 1, length=21) * max(abs(mat))
lp = levelplot(t(mat[ro, co]),
  aspect = "fill", at = at,
  scales = list(x = list(rot = 90)),
  colorkey = list(space = "left"))
print(lp)


###################################################
### chunk: chr1
###################################################
ps_chr = toTable(hgu95av2CHR)
ps_eg  = toTable(hgu95av2ENTREZID)
chr = merge(ps_chr, ps_eg)
chr = unique(chr[, colnames(chr)!="probe_id"])
head(chr)


###################################################
### chunk: chr2
###################################################
table(table(chr$gene_id))


###################################################
### chunk: chr4
###################################################
chr = chr[!duplicated(chr$gene_id), ]


###################################################
### chunk: EGsub
###################################################
isdiff = chr$gene_id %in% EGsub 
tab = table(isdiff, chr$chromosome)
tab
fisher.test(tab, simulate.p.value=TRUE) 
chisq.test(tab)


###################################################
### chunk: chrloc1
###################################################
chrloc = toTable(hgu95av2CHRLOC[featureNames(ALLsub)])
head(chrloc)


###################################################
### chunk: chrloc2
###################################################
table(table(chrloc$probe_id))


###################################################
### chunk: chrloc3
###################################################
strds = with(chrloc, 
  unique(cbind(probe_id, sign(start_location))))
table(strds[,2])


###################################################
### chunk: getMFchildren
###################################################
library("GO.db")
as.list(GOMFCHILDREN["GO:0008094"])


###################################################
### chunk: getMFoffspring
###################################################
as.list(GOMFOFFSPRING["GO:0008094"])


###################################################
### chunk: loadGOstats
###################################################
library("GOstats")


###################################################
### chunk: GOHyperGparam
###################################################
affyUniverse = featureNames(ALLfilt_af4bcr)
uniId = hgu95av2ENTREZID[affyUniverse]
entrezUniverse = unique(as.character(uniId))

params = new("GOHyperGParams",
    geneIds=EGsub, universeGeneIds=entrezUniverse,
    annotation="hgu95av2", ontology="BP",
    pvalueCutoff=0.001, conditional=FALSE,
    testDirection="over")


###################################################
### chunk: mfhyper
###################################################
mfhyper = hyperGTest(params)


###################################################
### chunk: figmfhyper
###################################################
hist(pvalues(mfhyper), breaks=50, col="mistyrose")


###################################################
### chunk: sigCats
###################################################
sum = summary(mfhyper, p=0.001)
head(sum)


###################################################
### chunk: GOTERM
###################################################
GOTERM[["GO:0032945"]]


###################################################
### chunk: bmCon eval=FALSE
###################################################
## library("biomaRt")
## head(listMarts())


###################################################
### chunk: choseMart eval=FALSE
###################################################
## mart = useMart("ensembl")


###################################################
### chunk: bmDataset eval=FALSE
###################################################
## head(listDatasets(mart))


###################################################
### chunk: bmChoseset eval=FALSE
###################################################
## ensembl = useDataset("hsapiens_gene_ensembl", 
##     mart=mart)


###################################################
### chunk: bmutrDo eval=FALSE
###################################################
## utr = getSequence(id=EGsub, seqType="3utr", 
##     mart=ensembl, type="entrezgene")
## utr[1,]


###################################################
### chunk: bmfilt eval=FALSE
###################################################
## head(listFilters(ensembl))


###################################################
### chunk: bmatt eval=FALSE
###################################################
## head(listAttributes(ensembl))


###################################################
### chunk: bmDom eval=FALSE
###################################################
## domains = getBM(attributes=c("entrezgene", "pfam", 
##     "prosite", "interpro"), filters="entrezgene", 
##     value=EGsub, mart=ensembl) 
## interpro = split(domains$interpro, domains$entrezgene)
## interpro[1]


###################################################
### chunk: loadDBanno
###################################################
library("hgu133a.db")
dbc = hgu133a_dbconn()


###################################################
### chunk: accessEnvs
###################################################
get("201473_at", hgu133aSYMBOL)
mget(c("201473_at","201476_s_at"), hgu133aSYMBOL)
hgu133aSYMBOL$"201473_at"
hgu133aSYMBOL[["201473_at"]]


###################################################
### chunk: GOcharacteristics1
###################################################
goCats = unlist(eapply(GOTERM, Ontology))
gCnums = table(goCats)[c("BP","CC", "MF")]


###################################################
### chunk: xtable eval=FALSE
###################################################
## library("xtable")
## xtable(as.matrix(gCnums), display=c("d", "d"),
##     caption="Number of GO terms per ontology.", 
##     label="ta:GOprops")


###################################################
### chunk: xtabledo
###################################################
library("xtable")
cat(print(xtable(as.matrix(gCnums), display=c("d", "d"),
    caption="Number of \\indexTerm[Gene Ontology (GO)]{GO} terms per ontology.", 
    label="ta:GOprops"), file="table1.tex"))


###################################################
### chunk: GOcharacteristics2
###################################################
query = "select ontology from go_term"
goCats = dbGetQuery(GO_dbconn(), query)
gCnums2 = table(goCats)[c("BP","CC", "MF")]
identical(gCnums, gCnums2)


###################################################
### chunk: GOterms2
###################################################
query = paste("select term from go_term where term",
    "like '%chromosome%'")
chrTerms = dbGetQuery(GO_dbconn(), query)
nrow(chrTerms)
head(chrTerms)


###################################################
### chunk: tfb
###################################################
query = paste("select go_id from go_term where",
    "term = 'transcription factor binding'")
tfb = dbGetQuery(GO_dbconn(), query)
tfbps =  hgu133aGO2ALLPROBES[[tfb$go_id]]
table(names(tfbps))


###################################################
### chunk: tfsol
###################################################
query = paste("select term from go_term where term",
    "like '%transcription factor%'")
tf = dbGetQuery(GO_dbconn(), query)
nrow(tf)
head(tf)


###################################################
### chunk: SYM2EG
###################################################
queryAlias = function(x) {
  it = paste("('", paste(x, collapse="', '"), "'", sep="")
  paste("select _id, alias_symbol from alias",
        "where alias_symbol in", it, ");")
}

queryGeneinfo = function(x) {
  it = paste("('", paste(x, collapse="', '"), "'", sep="")
  paste("select _id, symbol from gene_info where",
        "symbol in", it, ");")
}

queryGenes = function(x) {
  it = paste("('", paste(x, collapse="', '"), "'", sep="")
  paste("select * from genes where _id in", it,  ");")
}

findEGs = function(dbcon, symbols) {
  rs = dbSendQuery(dbcon, queryGeneinfo(symbols))
  a1 = fetch(rs, n=-1)
  stillLeft = setdiff(symbols, a1[,2])

  if( length(stillLeft)>0 ) {
    rs = dbSendQuery(dbcon, queryAlias(stillLeft))
    a2 = fetch(rs, n=-1)
    names(a2) = names(a1)
    a1 = rbind(a1, a2)
  } 
  
  rs = dbSendQuery(dbcon, queryGenes(a1[,1]))
  merge(a1, fetch(rs, n=-1))
}


###################################################
### chunk: findEGs
###################################################
findEGs(dbc, c("ALL1", "AF4", "BCR", "ABL"))


###################################################
### chunk: revMSym
###################################################
s1 = revmap(hgu133aSYMBOL)
s1$BCR


###################################################
### chunk: toTable
###################################################
toTable(hgu133aGO["201473_at"])


