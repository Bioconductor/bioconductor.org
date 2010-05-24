###################################################
### chunk: ReadAffy1 eval=FALSE
###################################################
## library("affy")
## myAB = ReadAffy()


###################################################
### chunk: ReadAffy2 eval=FALSE
###################################################
## myAB = ReadAffy(filenames=c("a1.cel", "a2.cel", "a3.cel"))


###################################################
### chunk: CLL
###################################################
library("CLL")
data("CLLbatch")
CLLbatch 


###################################################
### chunk: sns
###################################################
sampleNames(CLLbatch)


###################################################
### chunk: disease
###################################################
data("disease")
head(disease)


###################################################
### chunk: rownamesdisease
###################################################
rownames(disease) = disease$SampleID


###################################################
### chunk: removeCELsuffix
###################################################
sampleNames(CLLbatch) = sub("\\.CEL$", "", 
    sampleNames(CLLbatch))


###################################################
### chunk: mt
###################################################
mt = match(rownames(disease), sampleNames(CLLbatch))


###################################################
### chunk: ADF
###################################################
vmd = data.frame(labelDescription = c("Sample ID",  
    "Disease status: progressive or stable disease"))


###################################################
### chunk: samplesCLLbatch
###################################################
phenoData(CLLbatch) = new("AnnotatedDataFrame", 
    data = disease[mt, ], varMetadata = vmd)


###################################################
### chunk: dropmissing
###################################################
CLLbatch = CLLbatch[, !is.na(CLLbatch$Disease)]


###################################################
### chunk: qccalc
###################################################
library("affyQCReport")
saqc = qc(CLLbatch)


###################################################
### chunk: qcplot
###################################################
plot(saqc)


###################################################
### chunk: dist2
###################################################
dd = dist2(log2(exprs(CLLbatch)))


###################################################
### chunk: setupdistplot
###################################################
diag(dd) = 0
dd.row <- as.dendrogram(hclust(as.dist(dd)))
row.ord <- order.dendrogram(dd.row)
library("latticeExtra")
legend = list(top=list(fun=dendrogramGrob, 
    args=list(x=dd.row, side="top")))
lp = levelplot(dd[row.ord, row.ord], 
    scales=list(x=list(rot=90)), xlab="", 
    ylab="", legend=legend)



###################################################
### chunk: affyPLM
###################################################
library("affyPLM")
dataPLM = fitPLM(CLLbatch)


###################################################
### chunk: NUSE
###################################################
boxplot(dataPLM, main="NUSE", ylim = c(0.95, 1.22),
    outline = FALSE, col="lightblue", las=3, 
    whisklty=0, staplelty=0)


###################################################
### chunk: RLE
###################################################
Mbox(dataPLM, main="RLE", ylim = c(-0.4, 0.4), 
    outline = FALSE, col="mistyrose", las=3, 
    whisklty=0, staplelty=0)


###################################################
### chunk: dropBadarrays
###################################################
badArray = match("CLL1", sampleNames(CLLbatch))
CLLB = CLLbatch[, -badArray]


###################################################
### chunk: repRLENUSE
###################################################
dataPLMx = fitPLM(CLLB)
boxplot(dataPLM, main="NUSE", ylim = c(0.95, 1.3),
     outline = FALSE, col="lightblue", las=3,
     whisklty=0, staplelty=0)
Mbox(dataPLM, main="RLE", ylim = c(-0.4, 0.4), 
    outline = FALSE, col="mistyrose", las=3, 
    whisklty=0, staplelty=0)


###################################################
### chunk: rma
###################################################
CLLrma = rma(CLLB)


###################################################
### chunk: exprsDemo
###################################################
e = exprs(CLLrma)
dim(e)
dim(CLLrma)


###################################################
### chunk: howmanyprobesets
###################################################
dim(e)[1]
nrow(e)
dim(exprs(CLLrma))[1]
nrow(CLLrma)
length(featureNames(CLLrma))


###################################################
### chunk: pData1
###################################################
pData(CLLrma)[1:3,]


###################################################
### chunk: pData2
###################################################
table(CLLrma$Disease)


###################################################
### chunk: nsFilter
###################################################
CLLf = nsFilter(CLLrma, remove.dupEntrez=FALSE, 
    var.cutof =0.5)$eset


###################################################
### chunk: rowM
###################################################
CLLtt = rowttests(CLLf, "Disease")
names(CLLtt)


###################################################
### chunk: rowMeans
###################################################
a = rowMeans(exprs(CLLf))


###################################################
### chunk: figdvsa
###################################################
par(mfrow=c(1,2))
myPlot = function(...){
    plot(y = CLLtt$dm, pch = ".", ylim = c(-2,2), 
        ylab = "log-ratio", ...)
    abline(h=0, col="blue")
}

myPlot(x = a, xlab="average intensity")
myPlot(x = rank(a), xlab="rank of average intensity")


###################################################
### chunk: eBayesEx
###################################################
library("limma") 
design = model.matrix(~CLLf$Disease) 
CLLlim = lmFit(CLLf, design)
CLLeb = eBayes(CLLlim) 


###################################################
### chunk: compTstats
###################################################
plot(CLLtt$statistic, CLLeb$t[,2], pch=".")


###################################################
### chunk: volcano
###################################################
lod = -log10(CLLtt$p.value)
plot(CLLtt$dm, lod, pch=".", xlab="log-ratio", 
    ylab=expression(-log[10]~p))
abline(h=2)


###################################################
### chunk: volcanoeb
###################################################
plot(CLLtt$dm, -log10(CLLeb$p.value[,2]), pch=".", 
    xlab="log-ratio", ylab=expression(log[10]~p))
abline(h=2)


###################################################
### chunk: volcano2
###################################################
plot(CLLtt$dm, lod, pch=".", xlab="log-ratio", 
    ylab=expression(log[10]~p))
o1 = order(abs(CLLtt$dm), decreasing=TRUE)[1:25]
points(CLLtt$dm[o1], lod[o1], pch=18, col="blue")


###################################################
### chunk: table
###################################################
sum(CLLtt$p.value<=0.01)
sum(CLLeb$p.value[,2]<=0.01)


###################################################
### chunk: toptable
###################################################
tab = topTable(CLLeb, coef=2, adjust.method="BH", n=10)
genenames = as.character(tab$ID)


###################################################
### chunk: annotate1
###################################################
library("annotate")


###################################################
### chunk: annotate2
###################################################
annotation(CLLf)
library("hgu95av2.db")


###################################################
### chunk: entrezGeneAndSymbol
###################################################
ll = getEG(genenames, "hgu95av2")
sym = getSYMBOL(genenames, "hgu95av2")


###################################################
### chunk: output
###################################################
tab = data.frame(sym, signif(tab[,-1], 3))
htmlpage(list(ll), othernames=tab, 
    filename="GeneList1.html",
    title="HTML report", table.center=TRUE,
    table.head=c("Entrez ID",colnames(tab)))


###################################################
### chunk: browse
###################################################
browseURL("GeneList1.html")


###################################################
### chunk: 
###################################################
library("KEGG.db")
library("annaffy") 
atab = aafTableAnn(genenames, "hgu95av2.db", aaf.handler()) 
saveHTML(atab, file="GeneList2.html") 


###################################################
### chunk: annaffy2
###################################################
atab = aafTableAnn(genenames, "hgu95av2.db", 
    aaf.handler()[c(2,5,8,12)])
saveHTML(atab, file="GeneList3.html")


###################################################
### chunk: pmmm
###################################################
pms = pm(CLLB)
mms = mm(CLLB)


###################################################
### chunk: pmmmplot
###################################################
smoothScatter(log2(mms[,1]), log2(pms[,1]), 
    xlab=expression(log[2] * "MM values"),
    ylab=expression(log[2] * "PM values"), asp=1)
abline(a=0, b=1, col="red")


###################################################
### chunk: PMminusMM
###################################################
table(sign(pms-mms))


###################################################
### chunk: mmhist
###################################################
grouping = cut(log2(pms)[,1], breaks=c(-Inf, log2(2000), 
    Inf), labels=c("Low", "High")) 
library(geneplotter)
multidensity(log2(mms)[,1] ~ grouping, main="", xlab="", 
    col=c("red", "blue"), lwd=2)
legend("topright", levels(grouping), lty=1, lwd=2, 
    col=c("red", "blue"))


###################################################
### chunk: rmabgcorrect
###################################################
bgrma = bg.correct.rma(CLLB)
exprs(bgrma) = log2(exprs(bgrma))


###################################################
### chunk: vsn
###################################################
library("vsn")
bgvsn = justvsn(CLLB)


###################################################
### chunk: sel
###################################################
sel = sample(unlist(indexProbes(CLLB, "pm")), 500)
sel = sel[order(exprs(CLLB)[sel, 1])]


###################################################
### chunk: yoyryv
###################################################
yo = exprs(CLLB)[sel,1]
yr = exprs(bgrma)[sel,1]
yv = exprs(bgvsn)[sel,1]


###################################################
### chunk: bgrmavsn
###################################################
par(mfrow=c(1,3))
plot(yo, yr, xlab="Original", ylab="RMA", log="x", 
    type="l", asp=1)
plot(yo, yv, xlab="Original", ylab="VSN", log="x", 
    type="l", asp=1)
plot(yr, yv, xlab="RMA", ylab="VSN", type="l", asp=1)


###################################################
### chunk: vsnrma
###################################################
CLLvsn = vsnrma(CLLB)


###################################################
### chunk: vsndiff
###################################################
CLLvsnf = nsFilter(CLLvsn, remove.dupEntrez=FALSE, 
    var.cutoff=0.5)$eset
CLLvsntt = rowttests(CLLvsnf, "Disease")


###################################################
### chunk: inboth
###################################################
inboth = intersect(featureNames(CLLvsnf), 
                   featureNames(CLLf))


###################################################
### chunk: comparermavsn
###################################################
plot(CLLtt[inboth, "statistic"], 
     CLLvsntt[inboth, "statistic"], 
     pch=".", xlab="RMA", ylab="VSN", asp=1)


###################################################
### chunk: summary1
###################################################
pns = probeNames(CLLB)
indices = split(seq(along=pns), pns)


###################################################
### chunk: summary2
###################################################
length(indices)
indices[["189_s_at"]]


###################################################
### chunk: matplot
###################################################
colors = brewer.pal(8, "Dark2")
Index = indices[["189_s_at"]][seq(along=colors)]
matplot(t(pms[Index, 1:12]), pch="P", log="y", type="b", 
    lty=1, main="189_s_at", xlab="samples", 
    ylab=expression(log[2]~Intensity),
    ylim=c(50,2000), col=colors)
matplot(t(mms[Index, 1:12]), pch="M", log="y", type="b", 
    lty=3, add=TRUE, col=colors)


###################################################
### chunk: summary2
###################################################
newsummary = t(sapply(indices, function(j) 
    rowMedians(t(pms[j,]-mms[j,]))))
dim(newsummary)


###################################################
### chunk: ans
###################################################
colMeans(newsummary<0)*100


