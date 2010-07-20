###################################################
### chunk number 1: loadlib
###################################################
library("RbcBook1")
library("reposTools")
library("annaffy")
library("annotate")
library("hgu95av2")
data(aafExpr)


###################################################
### chunk number 2: getData
###################################################
 library("annaffy")
 data(aafExpr)
 gN <- geneNames(aafExpr)[41:50]
 syms <- unlist(mget(gN, hgu95av2SYMBOL))
 lls <- unlist(mget(gN, hgu95av2LOCUSID))
 syms


###################################################
### chunk number 3: makeTable
###################################################
gl <- list(gN, lls)
repository <- list("affy", "ll")
othernames<-list(syms)
head <- c("Probe ID", "Symbol", "LocusLink")
fName <- "out.html"
htmlpage(gl, fName, title="My Interesting Genes", othernames, head, 
   repository=repository)


###################################################
### chunk number 4: 
###################################################
probeids <- geneNames(aafExpr)
symbols <- aafSymbol(probeids, "hgu95av2")
symbols[55:57]


###################################################
### chunk number 5: 
###################################################
getText(symbols[55:57])


###################################################
### chunk number 6: 
###################################################
gos <- aafGO(probeids, "hgu95av2")
gos[[3]]


###################################################
### chunk number 7: Genbank1 eval=FALSE
###################################################
## gbs <- aafGenBank(probeids, "hgu95av2")
## getURL(gbs[[1]])


###################################################
### chunk number 8: Genbank2
###################################################
gbs <- aafGenBank(probeids, "hgu95av2")
url <- getURL(gbs[[1]])
cat(strbreak(url))


###################################################
### chunk number 9: Locuslink1 eval=FALSE
###################################################
## lls <- aafLocusLink(probeids, "hgu95av2") 
## getURL(lls[[2]])


###################################################
### chunk number 10: Locuslink2 eval=FALSE
###################################################
## lls <- aafLocusLink(probeids, "hgu95av2") 
## url <- getURL(lls[[2]])
## cat(strbreak(url))


###################################################
### chunk number 11: Cytobands1 eval=FALSE
###################################################
## bands <- aafCytoband(probeids, "hgu95av2") 
## getURL(bands[[2]])


###################################################
### chunk number 12: Cytobands2 eval=FALSE
###################################################
## bands <- aafCytoband(probeids, "hgu95av2") 
## url   <- getURL(bands[[2]])
## cat(strbreak(url))


###################################################
### chunk number 13: Pubmed1 eval=FALSE
###################################################
## pmids <- aafPubMed(probeids, "hgu95av2")
## getURL(pmids[[2]])


###################################################
### chunk number 14: Pubmed2 eval=FALSE
###################################################
## pmids <- aafPubMed(probeids, "hgu95av2")
## url   <- getURL(pmids[[2]])
## cat(strbreak(url))


###################################################
### chunk number 15: GO1 eval=FALSE
###################################################
## getURL(getURL(gos[[1]]))


###################################################
### chunk number 16: GO2 eval=FALSE
###################################################
## cat(strbreak(getURL(gos[[1]])))


###################################################
### chunk number 17: GOB1 eval=FALSE
###################################################
## getURL(gos[[1]][[4]])


###################################################
### chunk number 18: GOB2
###################################################
cat(strbreak(getURL(gos[[1]][[4]])))


###################################################
### chunk number 19: KEGG1 eval=FALSE
###################################################
## paths <- aafPathway(probeids, "hgu95av2")
## getURL(paths[[5]])


###################################################
### chunk number 20: KEGG2 eval=FALSE
###################################################
## paths <- aafPathway(probeids, "hgu95av2")
## url <- getURL(paths[[5]])
## cat(strbreak(url))


###################################################
### chunk number 21: loadmulttest
###################################################
library("multtest")


###################################################
### chunk number 22: 
###################################################
class <- as.integer(pData(aafExpr)$covar1) - 1


###################################################
### chunk number 23: 
###################################################
teststat <- mt.teststat(exprs(aafExpr), class)
index <- order(abs(teststat), decreasing = TRUE)
probeids <- geneNames(aafExpr)[index]


###################################################
### chunk number 24: 
###################################################
aaf.handler()


###################################################
### chunk number 25: 
###################################################
anncols <- aaf.handler()[c(1:3,8:9,11:13)]


###################################################
### chunk number 26: 
###################################################
anntable <- aafTableAnn(probeids[1:50], "hgu95av2", anncols)


###################################################
### chunk number 27: 
###################################################
saveHTML(anntable, "ex1.html", title = "Example Table without Data")


###################################################
### chunk number 28: 
###################################################
testtable <- aafTable("t-statistic" = teststat[index[1:50]], signed = TRUE)
table <- merge(anntable, testtable)


###################################################
### chunk number 29: 
###################################################
exprtable <- aafTableInt(aafExpr, probeids = probeids[1:50])
table <- merge(table, exprtable)
saveHTML(table, "example2.html", title = "Example Table with Data")


###################################################
### chunk number 30: 
###################################################
saveText(table, "example2.txt")


###################################################
### chunk number 31: kinasen
###################################################
kinases <- aafSearchText("hgu95av2", "Description", "kinase")
kinases[1:5]
length(kinases)


###################################################
### chunk number 32: 
###################################################
probes <- aafSearchText("hgu95av2", c("Description", "Pathway"),
                        c("ribosome", "polymerase"))
print(length(probes))


###################################################
### chunk number 33: 
###################################################
probes <- aafSearchText("hgu95av2", "Description",
                        c("DNA", "polymerase"), logic = "AND")
print(length(probes))


###################################################
### chunk number 34: 
###################################################
gbs <- c("AF035121", "AL021546", "AJ006123", "AL080082", "AI289489")
aafSearchText("hgu95av2", "GenBank", gbs)


