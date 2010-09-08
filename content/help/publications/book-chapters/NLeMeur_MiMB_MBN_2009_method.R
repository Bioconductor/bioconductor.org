###################################################
### chunk number 1: 
###################################################
options(width=56)


###################################################
### chunk number 2: installBioc eval=FALSE
###################################################
## source("http://bioconductor.org/biocLite.R") 
## biocLite() 


###################################################
### chunk number 3: installgraph eval=FALSE
###################################################
## biocLite(c("graph", "Rgraphviz"))


###################################################
### chunk number 4: loadpackage
###################################################
library("graph")
library("Rgraphviz")


###################################################
### chunk number 5: helpsearch eval=FALSE
###################################################
## class?graphNEL
## apropos("graphNEL")
## help.search("graphNEL")


###################################################
### chunk number 6: vector
###################################################
x <- c(10, 5, 6, 1.5, 9)
typeof(x)


###################################################
### chunk number 7: vector2
###################################################
y <-c(1:10) 
dim(y) <- c(2, 5) 
age <- c(33, 44, 55)
gender <- c("M", "F", "F")
myDF <- data.frame(age, gender)


###################################################
### chunk number 8: vector3
###################################################
is.factor(myDF$gender)
levels(myDF$gender)


###################################################
### chunk number 9: vector4
###################################################
x[1:2]
y[,1:2]
myDF[1,]


###################################################
### chunk number 10: randomgraph
###################################################
set.seed(123) 
g1 <-  randomEGraph(LETTERS[1:15], edges= 20) 
class(g1)
g1 


###################################################
### chunk number 11: graphPropeties
###################################################
library("RBGL")
s1 <-  connectedComp(g1)
s1
sp.between(g1, "A", "O")


###################################################
### chunk number 12: subgraph
###################################################
gsub <- subGraph(s1[[1]], g1)
minCut(gsub)


###################################################
### chunk number 13: rendergraph1
###################################################
g1 <- layoutGraph(g1) 
renderGraph(g1)


###################################################
### chunk number 14: rendergraph2 eval=FALSE
###################################################
## renderGraph(g1, 
##             graph.pars=list(edges= list(col= "lightblue", 
##                               lty= 1))) 


###################################################
### chunk number 15: graphpar
###################################################
default <- graph.par()


###################################################
### chunk number 16: rendergraph2
###################################################
graph.par(list(nodes = list(fill= "lightgray", 
                 textCol= "blue", cex= 1), 
               edges=list(col= "red", lty= 2))) 
renderGraph(g1)


###################################################
### chunk number 17: resetGraphpar
###################################################
graph.par(default) 


###################################################
### chunk number 18: rendergraph3 eval=FALSE
###################################################
## nodeRenderInfo(g1) <- list(fill= c(A= "lightyellow", 
##                              B= "lightblue")) 
## renderGraph(g1) 


###################################################
### chunk number 19: parseXml
###################################################
library("RpsiXML") 
site <- "ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psi25"
burl <- paste(site, "species/bacsu_small.xml", sep="/")
bfile <- paste(readLines(url(burl)), collapse= "")
bacsu <- parsePsimi25Interaction(bfile, 
                                 INTACT.PSIMI25, 
                                 verbose= FALSE) 


###################################################
### chunk number 20: tabulateIntTypes
###################################################
ii = interactions(bacsu)
table(sapply(ii, function(x) x@interactionType))
table(sapply(ii, function(x) x@expPubMed))


###################################################
### chunk number 21: ppiQA
###################################################
filterfun <- function(x) { 
  x@interactionType=="2 hybrid" && is.na(x@bait) && is.na(x@prey)
}
table(sapply(ii, filterfun))


###################################################
### chunk number 22: 
###################################################
outlier= sum(sapply(ii, filterfun)=="TRUE")
y2h=sum(sapply(ii, function(x) x@interactionType=='2 hybrid'))


###################################################
### chunk number 23: perPMID
###################################################
graphs <- separateXMLDataByExpt(bfile, INTACT.PSIMI25, 
                                type="indirect", directed= TRUE,
                                abstract= TRUE, verbose= FALSE) 
names(graphs) 


###################################################
### chunk number 24: QA
###################################################
library("ppiStats")
directedGraph <- sapply(graphs, isDirected)
aSummary <- createSummaryTables(graphs[directedGraph])


###################################################
### chunk number 25: latex summaryTable
###################################################
library("xtable")
table2 <- xtable(aSummary, display = c("s","d","d", "d","f", "f", "d", "f"),
       label = "tab:summaryTable",
       caption="A overview of summary statistics for all AP-MS data sets. 
                Each row displays the result reported in one paper (PMID).
                (VB) viable baits, (VP) viable prey ,  (VBP) viable bait/prey, 
                (TI) the total number of directed interactions.
                The statistics within this table can differentiate between experiment type.
                For example, the last column details the number of interactions per VB; 
                the experiment reported in paper 15045029 
                and 16796675 have remarkably high numbers of average interactions 
                per VB, and so would produce markedly distinct features 
                (i.e. protein complexes) from the others.")

print(table2, file="table/table2.tex", type="latex", size="small",
      tabular.environment="longtable", floating=FALSE, table.placement = "!htpb")


###################################################
### chunk number 26: IshikawaAbstract
###################################################
abstract(graphs[["16796675"]])


###################################################
### chunk number 27: IshikawaSummary
###################################################
Ishikawa <-  sapply(ii, function(x) x@expPubMed == 16796675)
sum(Ishikawa)
IshikawaInts <- ii[Ishikawa]
table(sapply(IshikawaInts, function(x) x@interactionType))


###################################################
### chunk number 28: IshikawaGraph
###################################################
IshikawaGraph <- graphs$`16796675`
IshikawaGraph


###################################################
### chunk number 29: IshikawaGraphProperties
###################################################
degree(IshikawaGraph)


###################################################
### chunk number 30: BaitPrey
###################################################
idViableProteins(IshikawaGraph)$VBP


###################################################
### chunk number 31: Cliques
###################################################
library("RBGL")
xu <- ugraph(IshikawaGraph) 
cls <- maxClique(xu)$maxCliques 
cs <- sapply(cls, length) 
maxC <- cls[cs== max(cs)]
maxC


###################################################
### chunk number 32: Clique
###################################################
ns <- nodes(IshikawaGraph) 
ncols <- rep("lightblue", length(ns))
ncols[ns %in% maxC[[1]]] <- "#FF8033"
nA <- makeNodeAttrs(IshikawaGraph, fillcolor= ncols, 
                    label="", width= 0.4, height= 0.4)
plot(IshikawaGraph, "neato", nodeAttrs= nA)


###################################################
### chunk number 33: CliqueAlone
###################################################
c4sub <- subGraph(maxC[[1]], IshikawaGraph)
ns <- nodes(c4sub)
ncols <- rep("lightblue", length(ns))
ncols[ns %in% maxC[[1]]] <- "#FF8033" 
nA <- makeNodeAttrs(c4sub, fillcolor= ncols)
plot(c4sub, "neato", nodeAttrs= nA)


###################################################
### chunk number 34: biomaRt
###################################################
library("biomaRt")
uniprot <- useMart("uniprot_mart")
uniprot <- useDataset("UNIPROT",mart=uniprot)
filters <- listFilters(uniprot)
attributes <- listAttributes(uniprot)
attributes[c(1:4,23), ]
selectedProt <- maxC[[1]]
getBM(attributes =  c(attributes[c(1:4, 23), 1]),
  filters = "accession", values = selectedProt, mart = uniprot) 


###################################################
### chunk number 35: kogeneName
###################################################
KOgene <- c("arcA", "appY", "oxyR", "soxS", "fnr")


###################################################
### chunk number 36: bESData
###################################################
load("data/bES.rda")
bES
pData(bES)[1:5,1:3]


###################################################
### chunk number 37: loadLimma
###################################################
library("limma")
aerobia <- bES[, bES$growth=="aerobic"]


###################################################
### chunk number 38: designMatrix
###################################################
aGeno <- as.factor(pData(aerobia)$genotype)
aDesign <- model.matrix(~ 0 + aGeno)
colnames(aDesign) <- levels(aGeno)


###################################################
### chunk number 39: contrastMatrix
###################################################
afit <- lmFit(aerobia,  aDesign) 

a.cont.matrix <- makeContrasts("appY"=appY-wt, 
                               "arcA"=arcA-wt, 
                               "fnr"=fnr-wt, 
                               "oxyR"=oxyR-wt, 
                               "soxS"=soxS-wt, levels=aDesign) 
afit2 <- contrasts.fit(afit, a.cont.matrix) 
afit2 <- eBayes(afit2) 
aRes <- c()
for(i in 1:ncol(a.cont.matrix)){
  aDE <- topTable(afit2, coef= i, n=500, adjust="BH", 
                  p.value=0.001)
  intRes <- cbind(colnames(a.cont.matrix)[i], 
                  aDE$genename , round(aDE$logFC, 2))
  aRes <- rbind(aRes, intRes)
}
colnames(aRes) <- c("KO","targets","logFC")


###################################################
### chunk number 40: tabulateDE
###################################################
table(aRes[,1])


###################################################
### chunk number 41: makegraphfromft
###################################################
graphAERO <- ftM2graphNEL(aRes[, c("KO","targets")],  
                        W=as.numeric(aRes[, "logFC"]), 
                          edgemode="directed")


###################################################
### chunk number 42: removeSelfLoop
###################################################
aeroG <- removeSelfLoops(graphAERO)
aeroG <- removeNode("arcA", aeroG)


###################################################
### chunk number 43: graphAero
###################################################
aeroG <- layoutGraph(aeroG, layoutType="dot")

col <- c(rep("orange", numEdges(aeroG)))
names(col) <- edgeNames(aeroG)
low <-  unlist(lapply(edgeWeights(aeroG), function(x) x < 0))
col[low] <- "blue"
edgeRenderInfo(aeroG) <-  list(col=col)  

KOlabel <- nodes(aeroG)%in%KOgene
nodeRenderInfo(aeroG)$labelY[!KOlabel] <- NA
renderGraph(aeroG)


###################################################
### chunk number 44: CategoryAnalysis
###################################################
library("GOstats")
library("Category")
library("org.EcK12.eg.db")

## Define the universe from Affy array
##name troubleshooting for universe
universe <- unique(fData(bES)$genename)
univEG <- unlist(mget(universe, 
                      org.EcK12.egSYMBOL2EG, ifnotfound=NA))
univEG <- univEG[!is.na(univEG)]

## set pvalue cutoff
pCutoff <- 0.001

## Test for each TF
GO <- vector("list", length=length(KOgene))
KEGG <- vector("list", length=length(KOgene))
names(GO) <- names(KEGG) <- KOgene

for(i in 1:length(KOgene)){
adjnode <- unlist(adj(graphAERO, KOgene[i]))
genelist <- c(KOgene[i], adjnode)
genelist <- unique(genelist)
genelist <- unlist(mget(genelist, org.EcK12.egSYMBOL2EG, ifnotfound=NA))
genelist <- genelist[!is.na(genelist)]

a.GOparams <- new("GOHyperGParams", geneIds = genelist, 
                  universeGeneIds = univEG, 
                  annotation = "org.EcK12.eg.db", ontology = "BP", 
                  pvalueCutoff = pCutoff, conditional=TRUE, 
                  testDirection = "over")
GO[[i]] <- hyperGTest(a.GOparams)

a.Kparams <- new("KEGGHyperGParams", geneIds = genelist, 
              universeGeneIds = univEG, annotation = "org.EcK12.eg.db",
               pvalueCutoff = pCutoff, testDirection = "over")
KEGG[[i]] <- hyperGTest(a.Kparams)
}


###################################################
### chunk number 45: GOresults
###################################################
sigGO <- lapply(GO, sigCategories)
sapply(sigGO, length)
lapply(GO["oxyR"], function(x) summary(x)[,c(2,3,7)])


###################################################
### chunk number 46: KEGGresults
###################################################
sigKEGG <-lapply(KEGG, sigCategories)
sapply(sigKEGG, length)


###################################################
### chunk number 47: Ironmetabolism
###################################################
gnodes <- nodes(graphAERO)
geneInKEGG <- geneIdsByCategory(KEGG[[2]])$`01053`
symbol <- unlist(mget(geneInKEGG, 
                      org.EcK12.egSYMBOL, ifnotfound=NA))
siderophore <- gnodes[gnodes%in%symbol]
# level of differentially expressed gene
as.data.frame(aRes[aRes[,2]%in%siderophore,])


