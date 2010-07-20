###################################################
### chunk: load0
###################################################
library("BiocCaseStudies")


###################################################
### chunk: loadLibs
###################################################
library("Biobase")
library("graph")
library("Rgraphviz")
library("RColorBrewer")
library("RBGL")
library("yeastExpData")
library(RpsiXML)


###################################################
### chunk: theData
###################################################
data("ccyclered")
data("litG")


###################################################
### chunk: litGNodes
###################################################
nodes(litG)[1:5]


###################################################
### chunk: exA
###################################################
class(ccyclered)
str(ccyclered)
dim(ccyclered)
names(ccyclered)


###################################################
### chunk: connComp
###################################################
cc = connectedComp(litG)
length(cc)
cclens = sapply(cc, length)
table(cclens)


###################################################
### chunk: cc7
###################################################
cc[[7]]


###################################################
### chunk: createSubG
###################################################
ord = order(cclens, decreasing=TRUE)
sg1 = subGraph(cc[[ord[1]]], litG)
sg2 = subGraph(cc[[ord[2]]], litG)


###################################################
### chunk: layoutSubG
###################################################
lsg1 = layoutGraph(sg1, layoutType="neato")
lsg2 = layoutGraph(sg2, layoutType="neato")


###################################################
### chunk: sg1
###################################################
renderGraph(lsg1)


###################################################
### chunk: sg2
###################################################
renderGraph(lsg2)


###################################################
### chunk: extraLayout
###################################################
lay12neato = layoutGraph(sg1, layoutType="dot")
renderGraph(lay12neato, 
     graph.pars=list(nodes=list(fillcolor="steelblue2")))
lay12twopi = layoutGraph(sg2, layoutType="twopi")
renderGraph(lay12twopi, 
     graph.pars=list(nodes=list(fillcolor="steelblue2")))


###################################################
### chunk: sps
###################################################
sps = sp.between(sg1, "YHR129C", "YOL039W")
pth = sps[[1]]$path_detail
pth


###################################################
### chunk: sglayout
###################################################
fill = rep("steelblue2", length(pth))
names(fill) = pth
nodeRenderInfo(lsg1) = list(fill=fill)
edges = paste(pth[-length(pth)], pth[-1], sep="~")
lwd = rep(5, length(edges))
col = rep("steelblue2", length(edges))
names(lwd) = names(col) = edges
edgeRenderInfo(lsg1) = list(col=col, lwd=lwd)


###################################################
### chunk: sg
###################################################
renderGraph(lsg1)


###################################################
### chunk: diam88sg
###################################################
allp = johnson.all.pairs.sp(sg1)


###################################################
### chunk: diamC
###################################################
max(allp)


###################################################
### chunk: howManyDiam
###################################################
sum(allp == max(allp))


###################################################
### chunk: geneToClusterList
###################################################
clusts = with(ccyclered, split(Y.name, Cluster))


###################################################
### chunk: makeClusterGraph
###################################################
cg = new("clusterGraph", clusters = clusts)


###################################################
### chunk: cgConnectedComp
###################################################
ccClust = connectedComp(cg)
length(ccClust)


###################################################
### chunk: intersect
###################################################
commonG = intersection(cg, litG)


###################################################
### chunk: numEdges
###################################################
numEdges(commonG)


###################################################
### chunk: defNodePerm
###################################################
nodePerm = function (g1, g2, B=1000) {
    n1 = nodes(g1)
    sapply(1:B, function(i) {
        nodes(g1) = sample(n1)
        numEdges(intersection(g1, g2))
    })
}


###################################################
### chunk: nodePermDo
###################################################
data("nPdist")
summary(nPdist)


###################################################
### chunk: nPdist
###################################################
barplot(table(nPdist))


###################################################
### chunk: loadlibs
###################################################
library("RpsiXML")
library("ppiStats")
library("apComplex")
library("xtable")


###################################################
### chunk: psi25int
###################################################
fold <- system.file("/extdata/psi25files", package="RpsiXML")
fn <- file.path(fold, "intact_2008_test.xml")
eg <- parsePsimi25Interaction(psimi25file=fn, 
                              psimi25source=INTACT.PSIMI25,
                              verbose=FALSE)


###################################################
### chunk: slotEntries
###################################################
slotNames(eg)


###################################################
### chunk: simpleSlots
###################################################
organismName(eg)
taxId(eg)
releaseDate(eg)


###################################################
### chunk: interactionSlot
###################################################
length(interactions(eg))
class(interactions(eg)[[1]])
interactions(eg)[[1]]
slotNames(interactions(eg)[[1]])


###################################################
### chunk: getBaitPrey
###################################################
egbait = sapply(interactions(eg), bait)
egprey = sapply(interactions(eg), prey)


###################################################
### chunk: baitVec
###################################################
egbait
egprey


###################################################
### chunk: interactors
###################################################
interactors(eg)[1:2]
bt <- egbait[4]
annBt <- interactors(eg)[[bt]]
annBt
xref(annBt)


###################################################
### chunk: psi25complex
###################################################
fn2 = file.path(fold, "intact_complexSample.xml")
comps = parsePsimi25Complex(fn2, INTACT.PSIMI25)
slotNames(comps)


###################################################
### chunk: complex1
###################################################
length(complexes(comps))
class(complexes(comps)[[1]])
slotNames(complexes(comps)[[1]])


###################################################
### chunk: showComplex1
###################################################
complexes(comps)[[1]]


###################################################
### chunk: intactgraph
###################################################
s1 = file.path(fold, "human_stelzl-2005-1_01.xml")
s2 = file.path(fold, "human_stelzl-2005-1_02.xml")
stelzl = separateXMLDataByExpt(xmlFiles=c(s1,s2), 
                               psimi25source=INTACT.PSIMI25, 
                               type="direct", directed=TRUE,
                               abstract=TRUE, verbose=FALSE)
class(stelzl[[1]])
stelzl1 = removeSelfLoops(stelzl[[1]])


###################################################
### chunk: abs
###################################################
abstract(stelzl1)


###################################################
### chunk: complexHG
###################################################
compXML <- file.path(fold, "intact_complexSample.xml")
pcHg <- buildPCHypergraph(compXML, INTACT.PSIMI25, 
                          split.by = "organismName") 
pcHg[[1]]


###################################################
### chunk: egHE
###################################################
complexes(pcHg[[1]])


###################################################
### chunk: activeBait
###################################################
deg = degree(stelzl1)
activeBait = names(which(deg$outDegree > 10 & 
    deg$outDegree<15))
proteins = union(activeBait, unlist(adj(stelzl1, 
    activeBait)))
stelzlSG = subGraph(proteins, stelzl1)


###################################################
### chunk: graphAtt
###################################################
graph.par(list(nodes=list(fill="steelblue", label="")))
baitCol = rep("yellow", length(activeBait))
names(baitCol) = activeBait
nodeRenderInfo(stelzlSG) <- list(fill=baitCol)


###################################################
### chunk: graph-stelzlSG1
###################################################
stelzlSG <- layoutGraph(stelzlSG, layoutType="neato")
renderGraph(stelzlSG)


###################################################
### chunk: graph-stelzlSG2
###################################################
stelzlSG <- layoutGraph(stelzlSG, layoutType="twopi")
renderGraph(stelzlSG)


