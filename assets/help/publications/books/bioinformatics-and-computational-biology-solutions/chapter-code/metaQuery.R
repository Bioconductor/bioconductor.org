###################################################
### chunk number 1: setup
###################################################
library("RbcBook1")


###################################################
### chunk number 2: Entrezex
###################################################
 ezURL <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
 library("XML")
 z <- xmlTreeParse(paste(ezURL, "einfo.fcgi", sep=""),
 isURL=TRUE, handlers=NULL, asTree=TRUE)

 dbL <- xmlChildren(z[[1]]$children$eInfoResult)$DbList

 dbNames<- xmlSApply(dbL, xmlValue)

 length(dbNames)

 dbNames[1:5]



###################################################
### chunk number 3: data0
###################################################
library("annotate")


###################################################
### chunk number 4: data
###################################################
data(eset, package="Biobase")
affys <- geneNames(eset)[491:500]
affys


###################################################
### chunk number 5: annotation
###################################################
library("hgu95av2")
ids <- getPMID(affys,"hgu95av2")
ids <- unlist(ids,use.names=FALSE)
ids <- unique(ids[!is.na(as.numeric(ids))])
length(ids)
ids[1:8]


###################################################
### chunk number 6: getabsts
###################################################
x <- pubmed(ids)
a <- xmlRoot(x)


###################################################
### chunk number 7: numabsts
###################################################
numAbst <- length(xmlChildren(a))


###################################################
### chunk number 8: searchAbst
###################################################
arts <- xmlApply(a, buildPubMedAbst)
class(arts[[7]])


###################################################
### chunk number 9: showArt
###################################################
arts[[7]]


###################################################
### chunk number 10: abstText
###################################################
## Retrieve the abstract text for this abstract
absts <- sapply(arts, abstText)
class(absts)


###################################################
### chunk number 11: abstWithcDNA
###################################################
found <- grep("cDNA",absts)
goodAbsts <- arts[found]
length(goodAbsts)


###################################################
### chunk number 12: 
###################################################
y <- genbank(ids[1:10], type="uid")
b <- xmlRoot(y)
class(b)


###################################################
### chunk number 13: abst2HTML
###################################################
fname <- tempfile()
pmAbst2HTML(goodAbsts, filename=fname)

fnameBase <- tempfile()
pmAbst2HTML(goodAbsts, filename=fnameBase, frames=TRUE)


###################################################
### chunk number 14: KEGGlibs
###################################################
library("KEGG")
library("KEGGSOAP")


###################################################
### chunk number 15: KEGGSOAPex
###################################################
  KEGGPATHID2NAME$"00020"


###################################################
### chunk number 16: KEGGSOAPexgg
###################################################
 genes <- get.genes.by.pathway("path:eco00020")
 enzymes <- get.enzymes.by.pathway("path:eco00020")


###################################################
### chunk number 17: KEGGSOAPans
###################################################
 enzymes[1:4]


###################################################
### chunk number 18: KEGGSOAPmotif
###################################################

 motifs <- get.motifs.by.gene("eco:b0002", "pfam")

 unlist(motifs[[1]][1:6])

 genes <- get.genes.by.motifs(c("pf:DnaJ", "ps:DNAJ_2"), 1, 10)

 genes[1:3]



###################################################
### chunk number 19: selprobe
###################################################

ps <- ls(env=hgu95av2ACCNUM)

myp <- ps[1001]

myA <- get(myp, hgu95av2ACCNUM)

myseq <- getSEQ(myA)
nchar(myseq)
substr(myseq, 1,10)



###################################################
### chunk number 20: refseq
###################################################


 rsp <- getSEQ("NP_004327")
 substr(rsp, 1, 10)

 rsn <- getSEQ("NM_004336")
 substr(rsn, 1, 10)


###################################################
### chunk number 21: getACC0
###################################################
library("hgu95av2probe")


###################################################
### chunk number 22: getAcc
###################################################

wp <- hgu95av2probe$Probe.Set.Name == myp
myPr <- hgu95av2probe[wp,]

library("Biostrings")
mybs <- NucleotideString(myseq, "DNA")

match1 <- matchDNAPattern(as.character(myPr[1,1]), mybs)
m1m <- as.matrix(match1)
m1m
myPr[1,5] == m1m[1,1]+13-1


