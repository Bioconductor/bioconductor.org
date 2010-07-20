###################################################
### chunk number 1: loadlibs
###################################################
library("Rgraphviz")
library("RBGL")
library("graph")
library("RColorBrewer")
library("geneplotter")
library("RbcBook1")

histcolor <- brewer.pal(7, "Pastel1")[5]
twocolors <- c(brewer.pal(11,"RdYlGn")[7], brewer.pal(11,"RdYlBu")[7])
sixcolors <- brewer.pal(6, "Dark2")
hundredcolors <- colorRampPalette(
      c(brewer.pal(9, "Pastel1"), brewer.pal(8, "Pastel2")))(100)



###################################################
### chunk number 2: my1stgraph
###################################################
library("graph")
myNodes <- c('s', 'p', 'q', 'r')
myEdges <- list(s=list(edges=c('p', 'q')),
                p=list(edges=c('p', 'q')),
                q=list(edges=c('p', 'r')),
                r=list(edges=c('s')))
g <- new('graphNEL', nodes=myNodes, edgeL=myEdges, edgemode='directed')
g


###################################################
### chunk number 3: plotmy1stgraph
###################################################
plot(g, nodeAttrs=makeNodeAttrs(g, fillcolor = twocolors[1]))


###################################################
### chunk number 4: nodesAndEdges
###################################################
nodes(g)
edges(g)
degree(g)


###################################################
### chunk number 5: prepareAdjAcc
###################################################
edges <- list(a=list(edges=2:3),
              b=list(edges=2:3),
              c=list(edges=c(2,4)),
              d=list(edges=1),
              e=list(edges=6,7),
              f=list(edges=7),
              g=list(edges=7))
g2 <- new('graphNEL', nodes=letters[1:7], edgeL=edges, edgemode='directed')
cc <- connectedComp(g2)
stopifnot(length(cc)==2)
colors <- rep(c("#D9EF8B","#E0F3F8"), listLen(cc))
names(colors) <- unlist(cc)


###################################################
### chunk number 6: plotAdjAcc
###################################################
plot(g2, nodeAttrs=list(fillcolor=colors))


###################################################
### chunk number 7: doAdjAcc
###################################################
adj(g2, "c")
acc(g2, c("b", "f"))


###################################################
### chunk number 8: prepareRandomEGraph
###################################################
set.seed(123)
nodeNames <- sapply(0:99, function(i) sprintf("N%02d", i))
rg <- randomEGraph(nodeNames, edges=50)


###################################################
### chunk number 9: cacheLayoutRandomEGraph
###################################################
rgNattrs <- makeNodeAttrs(rg, label="", fillcolor=hundredcolors)
layoutRandomEGraph <- cache("layoutRandomEGraph",
    agopen(rg, nodeAttrs=rgNattrs, layoutType="neato", name=""))


###################################################
### chunk number 10: plotRandomEGraph
###################################################
plot(layoutRandomEGraph)


###################################################
### chunk number 11: showDegreeRandomGraph eval=FALSE
###################################################
## deg <- degree(rg)
## hist(deg)


###################################################
### chunk number 12: prepareBinomRandomGraph
###################################################
size <- numNodes(rg)-1
prob <- numEdges(rg) / choose(numNodes(rg), 2)


###################################################
### chunk number 13: plotDegreeRandomGraph
###################################################
deg <- degree(rg)
hist(deg, breaks=seq(-0.5, max(deg)+0.5, by=1), col=histcolor, main="")
lines(0:size, numNodes(rg)*dbinom(0:size, size=size, prob=prob),
      pch=16, cex=3, type="b")


###################################################
### chunk number 14: connComp
###################################################
cc <- connComp(rg)
table(listLen(cc))


###################################################
### chunk number 15: connCompHide
###################################################
lcc <- sort(listLen(cc), decreasing=TRUE)
stopifnot(lcc[1]==18, lcc[2]==15)


###################################################
### chunk number 16: connCompExtract
###################################################
wh <- which.max(listLen(cc))
sg <- subGraph(cc[[wh]], rg)


###################################################
### chunk number 17: DFS
###################################################
dfs.res <- dfs(sg, node="N14", checkConn=TRUE)
nodes(sg)[dfs.res$discovered]


###################################################
### chunk number 18: BFS
###################################################
bfs.res <- bfs(sg, "N14")
nodes(sg)[bfs.res]


###################################################
### chunk number 19: plotbdfs1
###################################################
sgNattrs <- makeNodeAttrs(sg, fillcolor=rgNattrs$fillcolor[cc[[wh]]], shape="ellipse")
sgNattrs$label[1:numNodes(sg)] <- paste(" ", sgNattrs$label, " ", sep="")
plot(sg, "neato", nodeAttrs=sgNattrs, main="a) a connected subgraph", cex.main=2)


###################################################
### chunk number 20: plotbdfs2
###################################################
sgNattrs$label[nodes(sg)[dfs.res$discovered]] <- paste(1:numNodes(sg))
plot(sg, "neato", nodeAttrs=sgNattrs, main="b) DFS", cex.main=2)


###################################################
### chunk number 21: plotbdfs3
###################################################
sgNattrs$label[nodes(sg)[bfs.res]] <- paste(1:numNodes(sg))
plot(sg, "neato", nodeAttrs=sgNattrs, main="c) BFS", cex.main=2)


###################################################
### chunk number 22: conComp1
###################################################
sc <- strongComp(g2)
wc <- connComp(g2)


###################################################
### chunk number 23: conComp2
###################################################
stopifnot(length(sc)==4, length(wc)==2)
nattrs1 <- nattrs2 <- makeNodeAttrs(g2, fillcolor="")
for(ic in 1:length(sc))
  nattrs1$fillcolor[sc[[ic]]] <- sixcolors[ic+2]
for(ic in 1:length(wc))
  nattrs2$fillcolor[wc[[ic]]] <- sixcolors[ic]


###################################################
### chunk number 24: plotConComp1
###################################################
plot(g2, "dot", nodeAttrs=nattrs1)


###################################################
### chunk number 25: plotConComp2
###################################################
plot(g2, "dot", nodeAttrs=nattrs2)


###################################################
### chunk number 26: resetnEdges
###################################################
nEdges<-100


###################################################
### chunk number 27: prepareSP
###################################################
set.seed(123)
rg2 <- randomEGraph(nodeNames, edges=nEdges)
fromNode <- "N43"
toNode <- "N81"
sp <- sp.between(rg2, fromNode, toNode)
sp[[1]]$path
sp[[1]]$length


###################################################
### chunk number 28: layoutSP
###################################################
hilite <- brewer.pal(9, "Set1")[1]
layoutShortestPath <- cache("layoutShortestPath", {
  thepath <- sp[[paste(fromNode, ":", toNode, sep="")]]$path

  nattrs <- makeNodeAttrs(rg2, fillcolor="#c0c0c0", label="", shape="circle")
  nattrs$fillcolor[thepath] <- hilite

  lo <- agopen(rg2, nodeAttrs=nattrs, layoutType="neato", name="")

  ## color the edges in the path
  aged <- AgEdge(lo)
  for(i in seq(along=aged)) {
    mt <- match(head(aged[[i]]), thepath) - match(tail(aged[[i]]), thepath)
    if(!is.na(mt)){
      stopifnot(abs(mt)==1)   ## this condition should be fulfilled since 
                              ## 'thepath' is a shortest path
      aged[[i]]@color <- hilite
    }
  }
  lo@AgEdge <- aged
  lo
}) ## end of cache


###################################################
### chunk number 29: plotSP
###################################################
plot(layoutShortestPath)


###################################################
### chunk number 30: dijkstra.sp
###################################################
allsp<-dijkstra.sp(rg2, start=fromNode)
sum(!is.finite(allsp$distances))


###################################################
### chunk number 31: paranoia1
###################################################
## Make sure that our text matches reality.
stopifnot(identical(names(sp[[1]]), c("path", "length", "pweights")))
stopifnot(identical(names(allsp), c("distances", "penult", "start")))


###################################################
### chunk number 32: extractPath
###################################################
i1 <- match(fromNode,  nodes(rg2))
i2 <- match('N15', nodes(rg2))
pft <- RBGL::extractPath(i1, i2, allsp$penult)
nodes(rg2)[pft]


###################################################
### chunk number 33: johnsonAllPairs
###################################################
ap <- johnson.all.pairs.sp(rg2)
table(signif(ap, 3))


###################################################
### chunk number 34: plotHistDist1
###################################################
hist(allsp$distances[is.finite(allsp$distances)],
     breaks=seq(-0.5, max(allsp$distances[is.finite(allsp$distances)])+1, by=1),
     col=histcolor, main="", xlab=paste("distances from", fromNode))


###################################################
### chunk number 35: plotHistDist2
###################################################
stopifnot(identical(ap, t(ap)))
hist(ap, col=histcolor, main="", xlab="all pairwise distances")


###################################################
### chunk number 36: mstree1
###################################################
gr<-new("graphNEL", nodes=paste(1:9),
  edgeL=list("1"=list(edges=c(2,4)),
             "2"=list(edges=c(1,3,5)),
             "3"=list(edges=c(2,6)),
             "4"=list(edges=c(1,7,5)),
             "5"=list(edges=c(2,4,6,8)),
             "6"=list(edges=c(3,5,9)),
             "7"=list(edges=c(4,8)),
             "8"=list(edges=c(7,9,5)),
             "9"=list(edges=c(8,6))))
lgr1 <- lgr2 <- agopen(gr, layoutType="neato", name="",
             nodeAttrs=makeNodeAttrs(gr, shape="circle", color = "#404040"))


###################################################
### chunk number 37: mstree2
###################################################
mst <- mstree.kruskal(gr)
mst$edgeList


###################################################
### chunk number 38: mstree3
###################################################
stopifnot(identical(names(mst), c("edgeList", "weights", "nodes")))
count <- 0
for(i in 1:length(AgEdge(lgr2))) {
  theHead <- lgr2@AgEdge[[i]]@head
  theTail <- lgr2@AgEdge[[i]]@tail
  trHeads <- mst$nodes[mst$edgeList[2,]]
  trTails <- mst$nodes[mst$edgeList[1,]]
  if(any((theHead==trHeads)&(theTail==trTails))) {
    lgr2@AgEdge[[i]]@color <- hilite
    lgr2@AgEdge[[i]]@arrowhead <- "open"
    lgr2@AgEdge[[i]]@arrowsize <- "3"
  } else {
    if(any((theHead==trTails)&(theTail==trHeads))) {
      lgr2@AgEdge[[i]]@color <- hilite
      lgr2@AgEdge[[i]]@arrowtail <- "open"
      lgr2@AgEdge[[i]]@arrowsize <- "3"
    } else {
      count<-count+1
      lgr2@AgEdge[[i]]@color <- "white"
    }
  }
}
stopifnot(count == numEdges(gr)-ncol(mst$edgeList))


###################################################
### chunk number 39: mstreePlot1
###################################################
plot(lgr1)


###################################################
### chunk number 40: mstreePlot2
###################################################
plot(lgr2)


###################################################
### chunk number 41: plotNoAttrs
###################################################
library("Rgraphviz")
plot(g)


###################################################
### chunk number 42: FileDep
###################################################
data(FileDep)
fdnatt <- makeNodeAttrs(FileDep, label="", shape="circle",
  fillcolor=colorRampPalette(
      c(brewer.pal(9, "Pastel1"), brewer.pal(8, "Pastel2")))(numNodes(FileDep)))


###################################################
### chunk number 43: filedepdot
###################################################
plot(FileDep, "dot", nodeAttrs=fdnatt,  main="a) dot layout", name="")


###################################################
### chunk number 44: filedepneato
###################################################
plot(FileDep, "neato", nodeAttrs=fdnatt,  main="b) neato layout", name="")


###################################################
### chunk number 45: filedeptwopi
###################################################
plot(FileDep, "twopi", nodeAttrs=fdnatt,  main="c) twopi layout", name="")


###################################################
### chunk number 46: bipartite
###################################################
nodeN <-  c(paste("a", 1:5, sep=""), paste("b", 1:5, sep=""))
eL <- vector("list", length=10)
names(eL) <- nodeN
bpg <- new("graphNEL", nodes=nodeN, edgeL=eL, edgemode="directed")
bpg <- addEdge(rep(nodeN, c(3,3,1,3,3,0,0,0,0,0)),
        nodeN[c(6,8,9,7,8,10, 10, 6,7,8, 8,9,10)],
           bpg , rep(1, 13))

bpgatt <- makeNodeAttrs(bpg, label="", shape=c(rep("circle", 5), rep("box", 5)),
  fillcolor=c(brewer.pal(9, "YlOrRd")[3:7], brewer.pal(9, "Blues")[3:7]))


###################################################
### chunk number 47: bipartitedot
###################################################
plot(bpg, "dot", nodeAttrs=bpgatt, main="a) dot layout", name="")


###################################################
### chunk number 48: bipartiteneato
###################################################
plot(bpg, "neato", nodeAttrs=bpgatt, main="b) neato layout", name="")


###################################################
### chunk number 49: bipartitetwopi
###################################################
plot(bpg, "twopi", nodeAttrs=bpgatt, main="c) twopi layout", name="")


###################################################
### chunk number 50: getDefaultAttrs
###################################################
defAttrs <- getDefaultAttrs()
names(defAttrs)
names(defAttrs$node)
defAttrs$node$fillcolor


###################################################
### chunk number 51: plotAttrs
###################################################
nodeA <- list(fillcolor="lightblue")
edgeA <- list(color="goldenrod")
attrs <- getDefaultAttrs(list(node=nodeA, edge=edgeA))
plot(g, attrs=attrs)


###################################################
### chunk number 52: nodeAttrs eval=FALSE
###################################################
## cc <- connectedComp(g2)
## colors <- rep(c("#D9EF8B","#E0F3F8"), listLen(cc))
## names(colors) <- unlist(cc)
## plot(g2, nodeAttrs=list(fillcolor=colors))


###################################################
### chunk number 53: edgeAttrs
###################################################
globA   <- list(node=list(width="3", height="1", shape="box"))
nodeA <- list(label=c(p="YNL201C",q="Pph3", r="YBL046W",s="Spt5"),
              shape=c(p="ellipse",q="circle"))
edgeA <- list(label=c("p~q"="strong", "r~s"="weak"),
              color=c("p~q"="red", "r~s"="blue"))
plot(g, attrs=globA, edgeAttrs=edgeA, nodeAttrs=nodeA)


###################################################
### chunk number 54: agopen
###################################################
lg <- agopen(g, attrs=globA, edgeAttrs=edgeA, nodeAttrs=nodeA, name="ex1")


###################################################
### chunk number 55: AgNode
###################################################
ng <- AgNode(lg)
length(ng)
class(ng[[1]])
slotNames(ng[[1]])


###################################################
### chunk number 56: getPoints
###################################################
sapply(ng, function(x) getPoints(x@center))


###################################################
### chunk number 57: AgEdge
###################################################
eg<-AgEdge(lg)
slotNames(eg[[1]])
sapply(eg, function(x) (x@color))


###################################################
### chunk number 58: BezierCurve
###################################################
eg[[1]]@splines[[1]]


###################################################
### chunk number 59: drawNode
###################################################
prop    <- seq(0.2, 0.8, length=numNodes(g))
colors  <- brewer.pal(numNodes(g), "Set2")
names(prop) <- names(colors) <- nodes(g)
drawThermometerNode <- function(node) {
    x  <- getX(getNodeCenter(node))
    y  <- getY(getNodeCenter(node))
    w  <- getNodeLW(node)+getNodeRW(node)
    h  <- getNodeHeight(node)
    nm <- name(node)
    symbols(x, y, thermometers=cbind(w,h,prop[nm]), fg=colors[nm], 
            inches=FALSE, add=TRUE, lwd=3)
  }
lg <- agopen(g, name="")
plot(lg, drawNode=drawThermometerNode)


###################################################
### chunk number 60: getOldPar
###################################################
oldPar<-par(no.readonly=TRUE)


###################################################
### chunk number 61: zeugundsachen
###################################################
par(bg="#606060")

nr <- AgNode(layoutRandomEGraph)
x<-sapply(nr, function(j) getX(getNodeCenter(j)))
y<-sapply(nr, function(j) getY(getNodeCenter(j)))
dg <- degree(rg)+1
colors <- rev(colorRampPalette(brewer.pal(9, "GnBu"))(max(dg)))
plot(x, y, pch=16, col=colors[dg])

for(ed in AgEdge(layoutRandomEGraph))
  for(s in splines(ed)) {
    lines(bezierPoints(s), lty=3, col="#f0f0f0")
  }

legend(max(x), max(y), paste(0:max(degree(rg))), col=colors,
      xjust=1, text.col="white", pch=16)


###################################################
### chunk number 62: setOldPar
###################################################
stopifnot(all(range(degree(rg))==c(0,4)))
par(oldPar)


###################################################
### chunk number 63: imageMap eval=FALSE
###################################################
## imageMap(object, con, tags, imgname, width, height)


