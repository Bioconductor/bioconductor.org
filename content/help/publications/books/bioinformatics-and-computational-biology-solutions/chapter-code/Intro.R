###################################################
### chunk number 1: setup1
###################################################
library("RbcBook1")
library("Rgraphviz")
library("Biobase")
stopifnot(package.version("RbcBook1") >= package_version("0.2.1"))
stopifnot(package.version("Rgraphviz") >= package_version("1.5.6"))
stopifnot(package.version("CoCiteStats") >= package_version("0.5.3"))


###################################################
### chunk number 2: dataIMCA
###################################################
library("graph")
data("integrinMediatedCellAdhesion") 
class(IMCAGraph)


###################################################
### chunk number 3: SOS
###################################################
s <- acc(IMCAGraph, "SOS")


###################################################
### chunk number 4: degree
###################################################
deg <- degree(IMCAGraph)$outDegree
deg[which.max(deg)]


###################################################
### chunk number 5: tfG1
###################################################
library("GO")
library("GOstats")


###################################################
### chunk number 6: tfG2
###################################################
GOTERM$"GO:0003700"


###################################################
### chunk number 7: tfG3
###################################################
tfG   <- GOGraph("GO:0003700", GOMFPARENTS)


###################################################
### chunk number 8: tfG4
###################################################
library("Rgraphviz")
nL    <- nodes(tfG)
ggt   <- unlist(getGOTerm(nL))

labs <- character(length(nL))
names(labs) <- names(nL)
labs[sub("MF.", "", names(ggt))] <- ggt
labs["top"] <- "GO"

nattr <- makeNodeAttrs(tfG, 
   label     = labs[nodes(tfG)],
   shape     = "ellipse",
   fillcolor = "#f2f2f2", 
   fixedsize = FALSE)


###################################################
### chunk number 9: tfGplot
###################################################
plot(tfG, nodeAttrs=nattr)


###################################################
### chunk number 10: tf.children
###################################################
tfch <- GOMFCHILDREN$"GO:0003700"
tfchild <- mget(tfch, GOTERM)


###################################################
### chunk number 11: literatureGraphCalc
###################################################
library("humanLLMappings")

g1   <- "1736"        ##DKC1: 1736
paps <- humanLLMappingsLL2PMID[[g1]]

## throw out the Strausberg one
longpaps <- c("15302935", "12477932")
paps <- paps[!(paps %in% longpaps)]

genes <- lapply(paps, function(p) humanLLMappingsPMID2LL[[p]])
names(genes) <- paps

LLstring <- function(i) paste("LL", i, sep=":")
PMstring <- function(i) paste("PM", i, sep=":")

nd <- c(LLstring(unique(unlist(genes))),
        PMstring(paps))
ed        <- lapply(nd, function(z) list(edges=integer(0)))
names(ed) <- nd

for(i in 1:length(genes)) {
  p <- PMstring(names(genes)[i])
  g <- LLstring(genes[[i]])
  ed[[p]] <- list(edges=match(g, nd))
}
g <- new("graphNEL", nodes=nd, edgeL=ed, edgemode="directed")

nt <- match(substr(nodes(g), 1, 2), c("LL", "PM"))
nattr <-  makeNodeAttrs(g, 
             fillcolor=c("#8da0cb", "#e78ac3")[nt], 
             shape=c("ellipse", "rect")[nt],
             fixedsize=FALSE)


###################################################
### chunk number 12: literatureGraphPlot
###################################################
plot(g, "neato", nodeAttrs=nattr)


