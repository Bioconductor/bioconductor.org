###################################################
### chunk number 1: Setup
###################################################
library("RBGL")
library("Biobase")
require("Rgraphviz", quietly=TRUE)
library("RColorBrewer")
library("RbcBook1")
##
library("yeastExpData")
data(ccyclered)


###################################################
### chunk number 2: makeClusterGraph
###################################################
library("yeastExpData")
data(ccyclered)

clusts <-  split(ccyclered[["Y.name"]],
                 ccyclered[["Cluster"]])

cg1 <- new("clusterGraph", clusters = clusts)

ccClust <- connectedComp(cg1)


###################################################
### chunk number 3: litGraph
###################################################
 data(litG)
 ccLit  <- connectedComp(litG)
 cclens <- listLen(ccLit)
 table(cclens)


###################################################
### chunk number 4: litGraph.undercover
###################################################
ord <- order(cclens, decreasing=TRUE)
nrSingletons <- table(cclens)["1"]


###################################################
### chunk number 5: createSubG
###################################################
ord <- order(cclens, decreasing=TRUE)
sG1 <- subGraph(ccLit[[ord[1]]], litG)
sG2 <- subGraph(ccLit[[ord[2]]], litG)


###################################################
### chunk number 6: layoutSubG
###################################################
## Not sure how to use these attributes in a constructive manner:
attrs <- list(graph=list(ratio="10", nodesep="0"))
lsG1 <- cache("lsG1",
  agopen(sG1, layoutType="neato", nodeAttrs=makeNodeAttrs(sG1), name=""))
lsG2 <- cache("lsG2",
  agopen(sG2, layoutType="neato", nodeAttrs=makeNodeAttrs(sG2, fillcolor="#a6cee3"), name=""))


###################################################
### chunk number 7: sG1
###################################################
plot(lsG1)


###################################################
### chunk number 8: sG2
###################################################
plot(lsG2)


###################################################
### chunk number 9: intersect
###################################################
commonG <- intersection(cg1, litG)


###################################################
### chunk number 10: defNodePerm
###################################################
nodePerm <- function (g1, g2, B=1000) {
 n1 <- nodes(g1)
 sapply(1:B, function(i) {
    nodes(g1) <- sample(n1)
    numEdges(intersection(g1, g2))
  })
}
set.seed(123)


###################################################
### chunk number 11: nodePerm.do
###################################################
## FIXMEPUB: B=1000
## nPdist <- cache("nPdist", nodePerm(litG, cg1, B=100))
##
data(nPdist)


###################################################
### chunk number 12: nodePerm.show eval=FALSE
###################################################
## nPdist <- nodePerm(litG, cg1)


###################################################
### chunk number 13: histnPerm
###################################################
hist(nPdist, col=brewer.pal(7, "Pastel1")[5])


###################################################
### chunk number 14: Setup
###################################################
library("Biobase")
library("annotate")
library("GOstats")
library("xtable")
library("multtest")
library("hgu95av2")
library("genefilter")


###################################################
### chunk number 15: Example
###################################################
## subset of interest: 37+42 samples
data(ALL)
eset <- ALL[, intersect(grep("^B", as.character(ALL$BT)),
          which(as.character(ALL$mol) %in% c("BCR/ABL","NEG")))]

## intensities above 100 in at least 25% of the samples
f1 <- pOverA(.25, log2(100))
f2 <- function(x)(IQR(x)>0.5)
ff <- filterfun(f1, f2)
selected <- genefilter(eset, ff)
sum(selected)
esetSub <- eset[selected,]

cl <- as.numeric(esetSub$mol== "BCR/ABL")

## use multtest - permutation test, using the Welch statistic
##
## FIXMEPUB: B should be bigger - but it is too painful during the
## authoring phase
resT <- mt.maxT(exprs(esetSub), classlabel=cl, B=1000)
ord <- order(resT$index) #the original gene order
rawp <- resT$rawp[ord] # raw permutation $p$-values #names(resT)
names(rawp) <- geneNames(esetSub)

## get the gene names
gnames <- mget(geneNames(esetSub), env = hgu95av2SYMBOL, ifnotfound=NA)

top5 <- resT$index[1:5]
unlist(gnames[top5]) #$ sorted according to adjusted p-values


topsignif <- resT$index[resT$adjp<0.05]
gde <- unlist(gnames[topsignif])
  ## 'gde' for genes that are differentially expressed


###################################################
### chunk number 16: inducedGO
###################################################
gNll <- unlist(mget(names(gde), hgu95av2LOCUSID))
gNll <- as.character(gNll)
gGO  <- makeGOGraph(gNll, "MF")

nAgo <- makeNodeAttrs(gGO, shape="ellipse", label=substr(nodes(gGO), 7, 10))
gGopen <- cache("gGopen",
  agopen(gGO, recipEdges="distinct", layoutType="dot", nodeAttrs=nAgo, name=""))


###################################################
### chunk number 17: GOg1
###################################################
plot(gGopen)


###################################################
### chunk number 18: HyperG
###################################################
gNsLL <- unique(unlist(mget(names(gde),
              env=hgu95av2LOCUSID, ifnotfound=NA)))
gGhyp    <- GOHyperG(gNsLL)
## gGhyp.pv <- gGhyp$pv[nodes(gGO)]


###################################################
### chunk number 19: hypg2
###################################################
gGopenP <- gGopen
agnd    <- AgNode(gGopenP)
pthresh <- 0.1
stopifnot(!any(is.na(gGhyp$pv)))
for(i in seq(along=agnd)) {
  nm <- paste("GO:000", labelText(txtLabel(agnd[[i]])), sep="")
  mt <- match(nm, names(gGhyp$pv))
  if(!is.na(mt)) {
    agnd[[i]]@fillcolor <- ifelse(gGhyp$pv[mt] < pthresh, "#e31a1c", "#edf8fb")
    agnd[[i]]@txtLabel@labelText <- paste(signif(gGhyp$pv[mt], 2))
  } else {
    agnd[[i]]@fillcolor <- "white"
    agnd[[i]]@txtLabel@labelText <- ""
  }
}
gGopenP@AgNode <- agnd


###################################################
### chunk number 20: GOterms
###################################################
gg <- gGhyp$pv[gGhyp$pv < 0.1]
gg <- sort(gg)
gg.terms <- getGOTerm(names(gg))[["MF"]]

ggt <- unlist(gg.terms)
numCh <- nchar(ggt)
ggt2 <- substr(ggt, 1, 25)
ggt3 <- paste(ggt2, ifelse(numCh > 25, "...", ""), sep="")

##get counts
gg.counts <- gGhyp$goCounts[names(gg)]

ggMat <- matrix(c(names(gg.terms), ggt3, round(gg,3), gg.counts),
    byrow=FALSE, nc=4, dimnames=list(1:length(gg), c("GO ID",
   "Term", "p","n")))
print(xtable.matrix(ggMat,
caption="GO terms, $p$-values, and numbers of genes for a selection of GO categories.",
label="ta:GOggterms"), table.placement="bt")


###################################################
### chunk number 21: GOxx
###################################################
gc()
plot(gGopenP)


###################################################
### chunk number 22: PubMedEx1
###################################################
library("CoCiteStats")

gene1 <- "705"
gene2 <- "7216"

##compute gs and ps for use in the text below
ps <- sapply(mget(intersect(get(gene1,humanLLMappingsLL2PMID),
          get(gene2,humanLLMappingsLL2PMID)),
                  humanLLMappingsPMID2LL), length)
gs <- sapply(mget(c(gene1,gene2),humanLLMappingsLL2PMID),length)

twTable <- twowayTable(gene1, gene2, weights=F)
n1. <- twTable["n11"]+twTable["n12"]
n2. <- twTable["n21"]+twTable["n22"]
n.1 <- twTable["n11"]+twTable["n21"]
n.2 <- twTable["n12"]+twTable["n22"]
n <- n1.+ n2.

## hypergeometric distribution
  gene.adj <- function(gene)
    {
      pmids <- get(gene,humanLLMappingsLL2PMID)
      adj.genes <- unique(unlist(mget(pmids, humanLLMappingsPMID2LL)))
      adj.genes <- adj.genes[adj.genes!=gene]
      return(adj.genes)
    }

 hypergeo.pval <- function(gene, geneslist)
 {

       gene.adj <- gene.adj(gene)
        x <- length(intersect(geneslist,gene.adj))
        m <- length(gene.adj)
        n <-
        length(ls(humanLLMappingsLL2PMID)[!(ls(humanLLMappingsLL2PMID)
        %in%  gene.adj)])
        pvals <-
        dhyper(x,m,n,length(geneslist)) + phyper(x, m, n,
                                    length(geneslist), lower.tail=F)
        return(list(scores=x,pvals=pvals))
}

hypergeo.result <- hypergeo.pval(gene1,gene2)


###################################################
### chunk number 23: CoCiteEx1
###################################################

library("xtable")

LL2PMID <- as.list(humanLLMappingsLL2PMID)
LL2PMID2 <- LL2PMID[sapply(LL2PMID,function(x) !all(is.na(x)))]
PMID2LL <- as.list(humanLLMappingsPMID2LL)
PMID2LL2 <- PMID2LL[names(PMID2LL) !="NA"]

### create reduced (no NA) data environ
humanLLMappingsLL2PMID <- l2e(LL2PMID2)
humanLLMappingsPMID2LL <- l2e(PMID2LL2)

pLens <- paperLen(ls(humanLLMappingsLL2PMID))

numPapers <- length(pLens$Counts)
PaperLen <- pLens$Counts


result <- gene.geneslist.sig(gene1, gene2, numPapers, PaperLen, n.resamp=50)
result.mat <- NULL

as.char2 <- function(x, digits=NULL) {
    x <- as.character(round(x, digits))
    y <- vector(length=length(x))
    for ( i in 1:length(x)) {
        y[i] <- paste(c(x[i],
             ifelse(length(strsplit(x[i],"\\.")[[1]])==1, paste(c(".",
             rep(0,digits)),sep="",collapse=""), paste(rep(0,
             digits-nchar(strsplit(x[i],"\\.")[[1]][2])), sep="",
             collapse=""))), sep="",collapse="")
    }
    return(y)
}

for ( i in 1:nrow(result$statistic))
{
  result.mat <- rbind(result.mat, as.char2(result$statistic[i,], 4),
  tapply(as.char2(result$pval[i,],4), 1:length(result$pval[i,]),
  function(y) paste(c("(", y, ")"), sep="", collapse="")))
}

rownames(result.mat) <- c("None","","GS","","PS","","BOTH", "")
colnames(result.mat) <- colnames(result$statistic)
print(xtable(as.data.frame(result.mat),
caption=paste(c("PubMed co-citation: Locuslink ID ", gene1, " and ",
gene2, "."), sep="", collapse=""),  label="ta:CoCiteEx1",align=c("l",
                                                         rep("c",3))),
       table.placement="bt")



###################################################
### chunk number 24: PubMedEx2
###################################################

gene1 <- "10038"
gene2 <- "10039"

ps <- sapply(mget(intersect(get(gene1,humanLLMappingsLL2PMID),get(gene2,humanLLMappingsLL2PMID)),
humanLLMappingsPMID2LL), length)
gs <- sapply(mget(c(gene1,gene2),humanLLMappingsLL2PMID),length)

twTable <- twowayTable(gene1, gene2, weights=F)
n1. <- twTable["n11"]+twTable["n12"]
n2. <- twTable["n21"]+twTable["n22"]
n.1 <- twTable["n11"]+twTable["n21"]
n.2 <- twTable["n12"]+twTable["n22"]
n <- n1.+ n2.


###################################################
### chunk number 25: CoCiteEx2
###################################################

result <- gene.geneslist.sig(gene1, gene2, numPapers, PaperLen, n.resamp=50)
result.mat <- NULL

for ( i in 1:nrow(result$statistic))
{
  result.mat <- rbind(result.mat, as.char2(result$statistic[i,], 4), tapply(as.char2(result$pval[i,],4)
 , 1:length(result$pval[i,]), function(y) paste(c("(", y, ")"), sep="", collapse="")))
}

rownames(result.mat) <- c("None","","GS","","PS","","BOTH", "")
colnames(result.mat) <- colnames(result$statistic)
print(xtable(as.data.frame(result.mat),
caption=paste(c("PubMed co-citation: Locuslink ID ", gene1, " and ",
gene2, "."), sep="", collapse=""),
             label="ta:CoCiteEx2",align=c("l",rep("c",3))), table.placement="bt")


###################################################
### chunk number 26: getLLs
###################################################

intLLs <- unique(unlist(mget(names(gde), hgu95av2LOCUSID)))
intLLc <- as.character(intLLs)

library("humanLLMappings")


###################################################
### chunk number 27: getpapers
###################################################
papersByLL <- mget(intLLc, humanLLMappingsLL2PMID, ifnotfound=NA)
ncit <- sapply(papersByLL, length)
ncit


###################################################
### chunk number 28: getpapersizes
###################################################

 bcrnegpapers <- unique(unlist(papersByLL))
 pcocite <- mget(bcrnegpapers, humanLLMappingsPMID2LL)
 pclens <- sapply(pcocite, length)


###################################################
### chunk number 29: getRels
###################################################
num <- length(papersByLL)
grels <- vector("list", length= num)
names(grels) <- names(papersByLL)

for(i in 1:num) {
  curr <- papersByLL[[i]]
  grels[[i]] <- lapply(papersByLL, function(x) {
    mt <- match(x, curr, 0)
    if(any(mt > 0) )
        curr[mt]
    else
        NULL
})}

for(i in 1:num) grels[[i]] <- grels[[i]][-i]


###################################################
### chunk number 30: getpapers
###################################################
gr2 <- lapply(grels, function(x) {slen<-sapply(x, length); x[slen>0]})
table(unlist(gr2))


###################################################
### chunk number 31: combfun
###################################################
 LL2wts <- function(inList) {
     #first - papers for each LL
     pBLL <- mget(inList, humanLLMappingsLL2PMID, ifnotfound=NA)
     numL <- length(inList)
     ans<-NULL
     for( i in 1:numL) {
#        print(i)
        lls <- mget(as.character(pBLL[[i]]),
                    humanLLMappingsPMID2LL, ifnotfound=NA)
        lens <- sapply(lls, length)
        names(lens)<-NULL
        wts <- rep(1/lens, lens)
        wtsbyg <- split(wts, unlist(lls, use.names=FALSE))
        ans[[i]]<-sapply(wtsbyg, sum)
    }
    ans
 }

vv<- LL2wts(intLLc)



###################################################
### chunk number 32: compEW
###################################################
allLL <- unique(unlist(sapply(vv, names)))
bdrywts <- rep(0, length(allLL))
names(bdrywts) <- allLL
for(wvec in vv) bdrywts[names(wvec)] <- bdrywts[names(wvec)] + wvec

##drop the intLLc LLids
wts <- bdrywts[!(allLL %in% intLLc)]

sum(wts>1)
range(wts[wts>1])



###################################################
### chunk number 33: getchip
###################################################

LL95 <- unlist(as.list(hgu95av2LOCUSID))
bdryLL <- names(wts[wts>1])

onC <- match(bdryLL, LL95, 0)

unlist(mget(names(LL95[onC]), hgu95av2SYMBOL))



###################################################
### chunk number 34: pathstats
###################################################
library("hgu95av2")
library("annotate")
genel <- unlist(eapply(hgu95av2PATH, length))
table(genel)


###################################################
### chunk number 35: uniqLL
###################################################
pathLL <- eapply(hgu95av2PATH2PROBE, function(x) {
    LLs <- getLL(x, "hgu95av2")
    unique(LLs) })

pLens <- sapply(pathLL, length)
range(pLens)

uniqLL <- unique(unlist(pathLL,use.names=FALSE))


###################################################
### chunk number 36: LLpathwBG
###################################################

Amat <- sapply(pathLL, function(x) {
    mtch <- match(x, uniqLL)
    zeros <- rep(0, length(uniqLL))
    zeros[mtch] <- 1
    zeros})




###################################################
### chunk number 37: pathwayGraph
###################################################


pwGmat <- t(Amat) %*% Amat
diag(pwGmat) <- 0
pwG <- as(pwGmat, "graphNEL")



###################################################
### chunk number 38: connComppwG
###################################################
 ccpwG <- connectedComp(pwG)
 sapply(ccpwG, length)


###################################################
### chunk number 39: assertion
###################################################
zzA <- sapply(ccpwG, length)
stopifnot( length(zzA) == 6 )
stopifnot( sum(zzA == 1) == 5)


###################################################
### chunk number 40: pwNames
###################################################
library("KEGG")
for(i in ccpwG) {
   if(length(i) == 1)
    cat(get(i, KEGGPATHID2NAME), "\n")
    }


###################################################
### chunk number 41: pwGDegree
###################################################
hist(degree(pwG), col=brewer.pal(7, "Pastel1")[5])


###################################################
### chunk number 42: setup
###################################################

data(integrinMediatedCellAdhesion)

##nodes that do not correspond to genes have no LocusLink ID
llLens <- sapply(IMCAAttrs$LocusLink, length)
numLL <- sum(llLens > 0)


###################################################
### chunk number 43: hsaProbes
###################################################

hsa04510 <- hgu95av2PATH2PROBE$"04510"
hsaLLs <- getLL(hsa04510, "hgu95av2")



###################################################
### chunk number 44: resolveTechReps
###################################################
LLs <- unlist(sapply(IMCAAttrs$LocusLink, function(x) x[1]))

whProbe <- match(LLs, hsaLLs)

probeNames <- names(hsaLLs)[whProbe]
names(probeNames)<-names(LLs)
pN <- probeNames[!is.na(probeNames)]



###################################################
### chunk number 45: imca.do1
###################################################
nodeA <- makeNodeAttrs(IMCAGraph, width=1, height=.75,
           label=sub("^cell ", "", nodes(IMCAGraph)),
           shape="ellipse")

## the big node in the middle
j <- which( names(nodeA$label) ==  "Phosphatidylinositol signaling system" )
nodeA$label[j]  <- "Phosphatidyl-"
nodeA$shape[j]  <- "ellipse"
nodeA$width[j]  <- 4
nodeA$height[j] <- 2
nodeA$fillcolor[j] <- "white"

## some of the longer names
nc  <- nchar(nodeA$label)
sel <- (nc>4) & (seq(along=nc)!=j)
nodeA$width[sel]  <- nc[sel]/4

## overall attributes
attrs <- IMCAAttrs$defAttrs
attrs$graph$nodesep <- "0.1"
attrs$graph$ranksep <- "0.3"
attrs$graph$margin <- "0"

## all other nodes
nodeA$fillcolor[IMCAAttrs$nodeAttrs$fillcolor=="green"] <- brewer.pal(9, "YlGn")[3]

IMCg <- agopen(IMCAGraph, "", attrs=attrs,
     nodeAttrs=nodeA, subGList=IMCAAttrs$subGList)

IMCg@AgNode[[j]]@txtLabel@labelText <- ""


###################################################
### chunk number 46: imca
###################################################
plot(IMCg)
text( getX(getNodeCenter( AgNode(IMCg)[[j]])),
      getY(getNodeCenter( AgNode(IMCg)[[j]])) - seq(-2,2,length=4)*strheight("P", cex=.4),
      c("Phosphatidyl-", "inositol", "signaling", "system"), adj=c(0.5,1), cex=.4)


###################################################
### chunk number 47: imca.show eval=FALSE
###################################################
## IMCg <- agopen(IMCAGraph, "", attrs=IMCAAttrs$defAttrs,
##      nodeAttrs=IMCAAttrs$nodeAttrs, subGList=IMCAAttrs$subGList)
## plot(IMCg)


###################################################
### chunk number 48: getExprs
###################################################

 exprs <- exprs(eset)[pN,]

 bcrE <- exprs[,eset$mol=="BCR/ABL"]
 countsB <-  apply(bcrE, 1, function(x) {
        table(cut(x, breaks=c(-Inf, 6, 8.5, Inf)))
    })
 dimnames(countsB)[[2]] <- names(pN)

 bcrN <- exprs[,eset$mol=="NEG"]
 countsN <-  apply(bcrN, 1, function(x) {
        table(cut(x, breaks=c(-Inf, 6, 8.5, Inf)))
    })
 dimnames(countsN)[[2]] <- names(pN)

 nNodes <- nodes(IMCAGraph)
 ctsMatB <- matrix(0, ncol=length(nNodes), nrow=3)
 dimnames(ctsMatB) <- list(row.names(countsB), nNodes)
 ctsMatB[, names(pN)] <- countsB

  ctsMatN <- matrix(0, ncol=length(nNodes), nrow=3)
 dimnames(ctsMatN) <- list(row.names(countsN), nNodes)
 ctsMatN[, names(pN)] <- countsN




###################################################
### chunk number 49: colors
###################################################
pal1 <- brewer.pal(7, "RdPu")[c(2,4,6)]
nodeRadius <- median(sapply(agnd[colSums(ctsMatB)>0], getNodeRW))


###################################################
### chunk number 50: pieglyph
###################################################
## Construct a list of drawing functions for the Rgraphviz plot
plotPieChart <- function(curPlot, counts, main) {
  buildDrawing <- function(x) {
    force(x)
    if(all(x==0))
      return(drawAgNode)
    else
      return(function(node) {
             nodeCenter <- getNodeCenter(node)
             pieGlyph(x+0.01, radius=nodeRadius,
                      xpos=getX(nodeCenter),
                      ypos=getY(nodeCenter),
                      col=pal1)
           })
  } ## end of buildDrawing
  drawings <- apply(counts, 2, buildDrawing)
  plot(curPlot, drawNode=drawings, main=main)
  legend(par("usr")[2]*.99, par("usr")[3], 0,
         legend=c("No Data", "0-6", "6-8.5", "8.5+"),
         fill=c("#f0f0f0", pal1), xjust=1, yjust=0)
} ## end of definition 'plotPieChart'
IMCg2 <- IMCg
agnd  <- AgNode(IMCg2)
for(i in seq(along=agnd)) {
  agnd[[i]]@txtLabel@labelText <- ""
  agnd[[i]]@fillcolor <- "#f0f0f0"
}
IMCg2@AgNode <- agnd

## FIXME (wh 4.2.2005): we need proper replacement methods for AgNode and its slots.
## But if done naively this might become  rather inefficient.
## Do we need a different data structure from 'AgNode' lists, e.g.
## an 'AgNodes' object that represents sets of nodes and their properties
## by slots that are character and numeric vectors of the same length?


###################################################
### chunk number 51: plotpie1
###################################################
plotPieChart(IMCg2, ctsMatB,
             main="a) pie chart graph for BCR/ABL")


###################################################
### chunk number 52: plotpie2
###################################################
plotPieChart(IMCg2, ctsMatN,
             main="b) pie chart graph for NEG")


