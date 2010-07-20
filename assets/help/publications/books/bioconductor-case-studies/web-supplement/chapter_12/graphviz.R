###################################################
### chunk: pipExpData
###################################################
library("yeastExpData")
data("litG")
litG


###################################################
### chunk: ccomp
###################################################
library("RBGL")
cc = connectedComp(litG)
len = sapply(cc, length)
table(len)
ord = order(len, decreasing=TRUE)
g = subGraph(cc[[ord[4]]], litG)


###################################################
### chunk: plotDot
###################################################
library("Rgraphviz")
x = layoutGraph(g)
renderGraph(x)


###################################################
### chunk: noNames
###################################################
graph.par(list(nodes=list(label="", fill="lightgray"),
               edges=list(lwd=3)))


###################################################
### chunk: plotDot2
###################################################
renderGraph(x)


###################################################
### chunk: resetLabels
###################################################
graph.par(list(nodes=list(label=NULL)))


###################################################
### chunk: plotNeato
###################################################
graph.par(list(graph=list("cex.main"=2.5)))
x = layoutGraph(g, layoutType="neato")
renderGraph(x, graph.pars=list(graph=list(main="neato")))


###################################################
### chunk: plotTwopi
###################################################
x = layoutGraph(g, layoutType="twopi")
renderGraph(x, graph.pars=list(graph=list(main="twopi")))


###################################################
### chunk: plotCirco
###################################################
x = layoutGraph(g, layoutType="circo")
renderGraph(x, graph.pars=list(graph=list(main="circo")))


###################################################
### chunk: plotFdp
###################################################
x = layoutGraph(g, layoutType="fdp")
renderGraph(x, graph.pars=list(graph=list(main="fdp")))


###################################################
### chunk: nodePars
###################################################
nodeRenderInfo(x) = list(fill=c(YBR088C="red", 
                                YDR097="red")) 


###################################################
### chunk: edgeePars
###################################################
edgeRenderInfo(x) = list(lty=c("YBR088C~YMR167W"="dotted", 
                               "YDR097C~YMR167W"="dotted")) 


###################################################
### chunk: solPars eval=FALSE
###################################################
## ? renderParameters


###################################################
### chunk: layoutPar
###################################################
x = layoutGraph(x, attrs=list(node=list(shape="ellipse", 
    fixedsize=FALSE)))
renderGraph(x)


###################################################
### chunk: solShapes eval=FALSE
###################################################
## ? layoutParameters


###################################################
### chunk: naiveIMCA
###################################################
data("integrinMediatedCellAdhesion")
IMCAGraph
IMCAGraph = layoutGraph(IMCAGraph)
renderGraph(IMCAGraph)


###################################################
### chunk: longLabel
###################################################
names(labels) = labels = nodes(IMCAGraph)
nc = nchar(labels)
table(nc)
long = labels[order(nc, decreasing=TRUE)][1:4]
long


###################################################
### chunk: linefeed
###################################################
labels[long] = c(paste("Phosphatidyl-\ninositol\n",
    "signaling\nsystem", sep=""), "cell\nproliferation", 
    "cell\nmaintenance", "cell\nmotility") 


###################################################
### chunk: setLabel
###################################################
nodeRenderInfo(IMCAGraph) = list(label=labels)
renderGraph(IMCAGraph)


###################################################
### chunk: redoLayout1
###################################################
attrs = list(node=list(fixedsize=FALSE), 
    graph=list(rankdir="LR"))
width = c(2.5, 1.5, 1.5, 1.5)
height = c(1.5, 1.5, 1.5, 1.5)
names(width) = names(height) = long
nodeAttrs = list(width=width, height=height)
IMCAGraph = layoutGraph(IMCAGraph, attrs=attrs, 
    nodeAttrs=nodeAttrs)
renderGraph(IMCAGraph)


###################################################
### chunk: redoLayout2
###################################################
shape <- rep("rectangle", length(labels))
names(shape) <- labels
shape[long[1]] <- "ellipse" 
shape[c(long[2:4], "F-actin")] <- "plaintext"
nodeRenderInfo(IMCAGraph) <- list(shape=shape)
IMCAGraph <- layoutGraph(IMCAGraph, attrs=attrs, 
                         nodeAttrs=nodeAttrs)
renderGraph(IMCAGraph)


###################################################
### chunk: colorPlot
###################################################
colors = rep("lightgreen", length(nodes(IMCAGraph)))
names(colors) = nodes(IMCAGraph)
transp = c("ITGB", "ITGA", "MYO", "ACTN", "JNK", "p110", 
    "Phosphatidylinositol signaling system", 
    "PI5K", "MYO-P", "cell maintenance", "cell motility", 
    "F-actin", "cell proliferation")
colors[transp] = "transparent"
nodeRenderInfo(IMCAGraph) = list(fill=colors)
renderGraph(IMCAGraph)


###################################################
### chunk: recip
###################################################
IMCAGraph = layoutGraph(IMCAGraph, attrs=attrs, 
    nodeAttrs=nodeAttrs, recipEdges="distinct")
renderGraph(IMCAGraph)


###################################################
### chunk: subgraphs
###################################################
sg1 = subGraph(c("ITGA", "ITGB", "ILK", "CAV"), IMCAGraph)
sg2 = subGraph(c("cell maintenance", "cell motility", 
    "F-actin", "cell proliferation"), IMCAGraph)
sg3 = subGraph(c("ACTN", "VCL", "TLN", "PXN", "TNS", "VASP"), 
    IMCAGraph)
subGList = vector(mode="list", length=3)
subGList[[1]] = list(graph=sg1, attrs=c(rank="source"), 
    cluster=TRUE)
subGList[[2]] = list(graph=sg2, attrs=c(rank="sink"))
subGList[[3]] = list(graph=sg3)
IMCAGraph = layoutGraph(IMCAGraph, attrs=attrs, 
    nodeAttrs=nodeAttrs, subGList=subGList)
renderGraph(IMCAGraph)


###################################################
### chunk: solSubgraph
###################################################
sg4 =  subGraph(c("GRB2", "SOS", "Ha-Ras", "Raf", 
    "MEK", "ERK"), IMCAGraph)
subGList = append(subGList, list(list(graph=sg4)))
IMCAGraph = layoutGraph(IMCAGraph, attrs=attrs, 
    nodeAttrs=nodeAttrs, subGList=subGList)
renderGraph(IMCAGraph)


###################################################
### chunk: litGraph1
###################################################
library("org.Hs.eg.db")
library("biocGraph")
g1   = "1736"
paps = org.Hs.egPMID[[g1]]
genes = mget(paps, org.Hs.egPMID2EG)
names(genes) = paps
len = sapply(genes, length)
table(len)


###################################################
### chunk: litGraph2
###################################################
sel = len < 5 
genes = genes[sel] 
paps = paps[sel]
LLstring = function(i) paste("LL", i, sep=":")
PMstring = function(i) paste("PM", i, sep=":")
nd = c(LLstring(unique(unlist(genes))),
    PMstring(paps))
ed = lapply(nd, function(z) list(edges=integer(0)))
names(ed) = nd


###################################################
### chunk: litGraph3
###################################################
for(i in 1:length(genes)) {
    p = PMstring(names(genes)[i])
    g = LLstring(genes[[i]])
    ed[[p]] = list(edges=match(g, nd))
}
bpg = new("graphNEL", nodes=nd, edgeL=ed, 
    edgemode="directed")


###################################################
### chunk: renderAttrs
###################################################
nt = match(substr(nodes(bpg), 1, 2), c("LL", "PM"))
fill = c("lightblue", "salmon")[nt]
shape = c("ellipse", "rect")[nt]
names(fill) = names(shape) = nodes(bpg)
attrs = list(node=list(fixedsize=TRUE, shape=shape))
nodeRenderInfo(bpg) = list(fill=fill)
graph.par(list(nodes=list(fontsize=10)))


###################################################
### chunk: drawPlot eval=FALSE
###################################################
## imgname = "graph.png"
## png(imgname, width=1024, height=768)
## bpg = layoutGraph(bpg, layoutType="neato", attrs=attrs)
## bpg = renderGraph(bpg)
## dev.off()


###################################################
### chunk: fileConnection eval=FALSE
###################################################
## library(geneplotter)
## fhtml = "index.html"
## con = openHtmlPage(fhtml, paste("PubMed co-citations of",
##     "gene '1736' Please click on the nodes"))


###################################################
### chunk: pmquery
###################################################
pnodes = nodes(bpg)[nt==2]
links = character(length(pnodes))
tooltips = pnodes
links =  paste("http://www.ncbi.nih.gov/entrez/query.fcgi",
    "?tool=bioconductor&cmd=Retrieve&db=PubMed&",
    "list_uids=", gsub("PM:", "", pnodes), sep="")
names(links) = names(tooltips) = pnodes
tags = list(HREF=links, TITLE=tooltips)


###################################################
### chunk: indexhtml
###################################################
imageMap(bpg, con=con, tags=tags, imgname=imgname)
closeHtmlPage(con)


###################################################
### chunk: browseurl eval=FALSE
###################################################
## fhtml
## browseURL(fhtml)


