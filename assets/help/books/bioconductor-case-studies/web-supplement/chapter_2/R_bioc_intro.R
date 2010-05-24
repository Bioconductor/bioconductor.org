###################################################
### chunk: aprop
###################################################
apropos("mean")
find("mean")


###################################################
### chunk: helpsearch eval=FALSE
###################################################
## help.search("mean")


###################################################
### chunk: exHelp1print eval=FALSE
###################################################
## apropos("plot")


###################################################
### chunk: exHelp2 eval=FALSE
###################################################
## help.search("mann-whitney")


###################################################
### chunk: exHelp3 eval=FALSE
###################################################
## library("Biobase")
## openVignette("Biobase")


###################################################
### chunk: biocLiteEx eval=FALSE
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite(c("graph", "xtable"))


###################################################
### chunk: updatepackagesEx eval=FALSE
###################################################
## source("http://bioconductor.org/biocLite.R")
## update.packages(repos=biocinstallRepos(), ask=FALSE)


###################################################
### chunk: exSession eval=FALSE
###################################################
## sessionInfo()


###################################################
### chunk: ex0sol1
###################################################
x = c(0.1, 1.1, 2.5, 10)
y = 1:100
z = y < 10
pets = c(Rex="dog", Garfield="cat", Tweety="bird")


###################################################
### chunk: ex0sol2
###################################################
2 * x + c(1,2)


###################################################
### chunk: ex0sol3
###################################################
##logical
y[z]
## integer
y[1:4]
y[-(1:95)]
## character
pets["Garfield"]


###################################################
### chunk: ex0sol4
###################################################
m = matrix(1:12, ncol=4)
m[1,3]


###################################################
### chunk: ex0sol5
###################################################
l = list(name="Paul", sex=factor("male"), age=35)
l$name
l[[3]]


###################################################
### chunk: simpleFun
###################################################
sq1 = function(x) return(x*x)
sq2 = function(x) x*x


###################################################
### chunk: ppc
###################################################
ppc = function(x) paste("^", x, sep="")


###################################################
### chunk: map
###################################################
library("hgu95av2.db")
hgu95av2MAP$"1001_at"


###################################################
### chunk: eApply1
###################################################
myPos = eapply(hgu95av2MAP, function(x) grep("^17p", x, 
    value=TRUE))
myPos = unlist(myPos)
length(myPos)


###################################################
### chunk: eApply2
###################################################
f17p = function(x) grep("^17p", x, value=TRUE)
myPos2 = eapply(hgu95av2MAP, f17p)
myPos2 = unlist(myPos2)
identical(myPos, myPos2)


###################################################
### chunk: ans2
###################################################
myFindMap = function(mapEnv, which) {
    myg = ppc(which)
    a1 = eapply(mapEnv, function(x)
        grep(myg, x, value=TRUE))
    unlist(a1)
}


###################################################
### chunk: envEx
###################################################
e1 = new.env(hash=TRUE)
e1$a = rnorm(10)
e1$b = runif(20)
ls(e1)
xx = as.list(e1)
names(xx)
rm(a, envir=e1)


###################################################
### chunk: ex2Sol1
###################################################
theEnv = new.env(hash=TRUE)
theEnv$locations = myFindMap(hgu95av2MAP, 18)


###################################################
### chunk: ex2Sol2
###################################################
theEnv$strip = function(x) gsub("18", "", x)


###################################################
### chunk: ex2Sol3
###################################################
myExtract = function(env)  env$strip(env$locations)   
myExtract(theEnv)[1:5]


###################################################
### chunk: convert eval=FALSE
###################################################
## library(convert)
## as(object, "ExpressionSet")


###################################################
### chunk: read-table-geneData
###################################################
dataDirectory = system.file("extdata", package="Biobase")
exprsFile = file.path(dataDirectory, "exprsData.txt")
exprs = as.matrix(read.table(exprsFile, header=TRUE, 
    sep="\t", row.names=1, as.is=TRUE))


###################################################
### chunk: exprsFile eval=FALSE
###################################################
## exprsFile = "c:/path/to/exprsData.txt"


###################################################
### chunk: geneData-peak
###################################################
class(exprs)
dim(exprs)
colnames(exprs)
head(exprs)


###################################################
### chunk: pData
###################################################
pDataFile = file.path(dataDirectory, "pData.txt")
pData = read.table(pDataFile,
    row.names=1, header=TRUE, sep="\t")
dim(pData)
rownames(pData)
summary(pData)


###################################################
### chunk: geneCovariate-geneData-name-match
###################################################
all(rownames(pData) == colnames(exprs))


###################################################
### chunk: rtclass
###################################################
class(pData)


###################################################
### chunk: colnames
###################################################
names(pData)


###################################################
### chunk: sapplyClasses
###################################################
sapply(pData, class)


###################################################
### chunk: simpleSubsetting
###################################################
pData[c(15, 20), c("gender", "type")]
pData[pData$score > 0.8,]


###################################################
### chunk: metadata-create
###################################################
metadata = data.frame(labelDescription=c("Patient gender", 
    "Case/control status", "Tumor progress on XYZ scale"),
     row.names=c("gender", "type", "score"))


###################################################
### chunk: AnnotatedDataFrame
###################################################
adf = new("AnnotatedDataFrame", data=pData, 
    varMetadata=metadata)
adf


###################################################
### chunk: AnnotatedDataFrame-subset
###################################################
head(pData(adf))
adf[c("A","Z"), "gender"]
pData(adf[adf$score > 0.8,])


###################################################
### chunk: annotation
###################################################
annotation = "hgu95av2"


###################################################
### chunk: R.MIAME
###################################################
experimentData = new("MIAME", name="Pierre Fermat",
    lab="Francis Galton Lab", 
    contact="pfermat@lab.not.exist",
    title="Smoking-Cancer Experiment",
    abstract="An example ExpressionSet",
    url="www.lab.not.exist",
    other=list(notes="Created from text files"))


###################################################
### chunk: ExpressionSetFinally
###################################################
exampleSet = new("ExpressionSet", exprs=exprs, 
    phenoData=adf, experimentData=experimentData,
    annotation="hgu95av2")


###################################################
### chunk: ExpressionSet-minimal
###################################################
minimalSet = new("ExpressionSet", exprs=exprs)


###################################################
### chunk: helpExpressionSet eval=FALSE
###################################################
## help("ExpressionSet-class")


###################################################
### chunk: showExpressionSet
###################################################
exampleSet


###################################################
### chunk: usingDollar
###################################################
exampleSet$gender[1:5]
exampleSet$gender[1:5] == "Female"


###################################################
### chunk: featureNames
###################################################
featureNames(exampleSet)[1:5]


###################################################
### chunk: sampleNames
###################################################
sampleNames(exampleSet)[1:5]
varLabels(exampleSet)


###################################################
### chunk: exprs
###################################################
mat = exprs(exampleSet)
dim(mat)
adf = phenoData(exampleSet)
adf


###################################################
### chunk: first10
###################################################
vv = exampleSet[1:5, 1:3]
dim(vv)
featureNames(vv)
sampleNames(vv)


###################################################
### chunk: males
###################################################
males = exampleSet[, exampleSet$gender == "Male"]
males


###################################################
### chunk: dotplot
###################################################
x = exprs(exampleSet[, 1])
y = exprs(exampleSet[, 3])
plot(x=x, y=y, log="xy")


###################################################
### chunk: ex5Sol1
###################################################
plot(x=x, y=y, log="xy", 
     xlab="gene expression sample #1", 
     ylab="gene expression sample #3",
     main="scatterplot of expression intensities", 
     pch=20)
abline(a=0, b=1)


###################################################
### chunk: pregcc
###################################################
library("CLL")
library("matchprobes")
library("hgu95av2probe")
library("hgu95av2cdf")
library("RColorBrewer")
data("CLLbatch")
bases = basecontent(hgu95av2probe$sequence)


###################################################
### chunk: matchgcc
###################################################
iab = with(hgu95av2probe, xy2indices(x, y, 
    cdf="hgu95av2cdf"))
probedata = data.frame(
    int=rowMeans(log2(exprs(CLLbatch)[iab, ])),
    gc=bases[, "C"] + bases[, "G"])


###################################################
### chunk: boxplotgcc
###################################################
colorfunction = colorRampPalette(brewer.pal(9, "GnBu"))
mycolors = colorfunction(length(unique(probedata$gc)))
label = expression(log[2]~intensity)
boxplot(int ~ gc, data=probedata, col=mycolors, 
    outline=FALSE, xlab="Number of G and C", 
    ylab=label, main="")


###################################################
### chunk: dens1
###################################################
tab = table(probedata$gc)
gcUse = as.integer(names(sort(tab, decreasing=TRUE)[1:10]))
gcUse


###################################################
### chunk: densplotgcc
###################################################
library("geneplotter")
multidensity(int ~ gc, data=subset(probedata, 
    gc %in% gcUse), xlim=c(6, 11), 
    col=colorfunction(12)[-(1:2)],
    lwd=2, main="", xlab=label)


###################################################
### chunk: ex6Sol1
###################################################
multiecdf(int ~ gc, data=subset(probedata, gc %in% gcUse), 
    xlim=c(6, 11), col=colorfunction(12)[-(1:2)],
    lwd=2, main="", ylab="ECDF")


