###################################################
### chunk number 1: loadlibs
###################################################
library("RbcBook1")
library("Rgraphviz")
library("RBGL")
library("graph")
library("RColorBrewer")
library("geneplotter")
library("GOstats")
library("GO")

histcolor <- brewer.pal(7, "Pastel1")[5]
twocolors <- c(brewer.pal(11,"RdYlGn")[7], brewer.pal(11,"RdYlBu")[7])
sixcolors <- brewer.pal(6, "Dark2")
hundredcolors <- colorRampPalette(
      c(brewer.pal(9, "Pastel1"), brewer.pal(8, "Pastel2")))(100)

graphVizAttrs <- getDefaultAttrs()
graphVizAttrs$node$fillcolor <- twocolors[1]


###################################################
### chunk number 2: sgraph
###################################################
myNodes <- c('s', 'p', 'q', 'r')
myEdges <- list(s=list(edges=c('p', 'q')),
                p=list(edges=c('p', 'q')),
                q=list(edges=c('p', 'r')),
                r=list(edges=c('s')))
g <- new('graphNEL', nodes=myNodes, edgeL=myEdges, edgemode='directed')
g


###################################################
### chunk number 3: sgraph
###################################################
plot(g, attrs=graphVizAttrs)


###################################################
### chunk number 4: prepareSetOps
###################################################
V <- letters[1:4]
set.seed(4713)
ug1 <- randomGraph(V, 1, .55)
ug2 <- randomGraph(V, 1, .55)

myPlot <- function(g, i, tit) {
  plot(g, nodeAttrs=makeNodeAttrs(g, fillcolor=sixcolors[i]), main=tit)
}


###################################################
### chunk number 5: setOps1
###################################################
par(cex.main=2)
myPlot(ug1, 1, "ug1")


###################################################
### chunk number 6: setOps2
###################################################
par(cex.main=2)
myPlot(complement(ug1), 2, "complement(ug1)")


###################################################
### chunk number 7: setOps3
###################################################
par(cex.main=2)
myPlot(ug2, 3, "ug2")


###################################################
### chunk number 8: setOps4
###################################################
par(cex.main=2)
myPlot(complement(ug2), 4, "complement(ug2)")


###################################################
### chunk number 9: setOps5
###################################################
par(cex.main=2)
myPlot(intersection(ug1,ug2), 5, "intersection(ug1, ug2)")


###################################################
### chunk number 10: setOps6
###################################################
par(cex.main=2)
myPlot(union(ug1,ug2), 6, "union(ug1, ug2)")


