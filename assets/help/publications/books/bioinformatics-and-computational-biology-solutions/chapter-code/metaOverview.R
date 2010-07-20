###################################################
### chunk number 1: setup
###################################################
library("RbcBook1")
library("GO")
library("KEGG")
library("hgu133a")
library("GOstats")
library("cMAP")
library("annotate")


###################################################
### chunk number 2: lkhgu0
###################################################
library("hgu95av2")
ls("package:hgu95av2")



###################################################
### chunk number 3: getBAD
###################################################

gsyms <- unlist(as.list(hgu95av2SYMBOL))
whBAD <- grep("^BAD$", gsyms)
gsyms[whBAD]


###################################################
### chunk number 4: keggBAD
###################################################

 BADpath <- hgu95av2PATH$"1861_at"
 mget(BADpath, KEGGPATHID2NAME)
 allProbes <- mget(BADpath, hgu95av2PATH2PROBE)
 sapply(allProbes, length)


###################################################
### chunk number 5: lkyea0
###################################################
library("YEAST")
ls("package:YEAST")



###################################################
### chunk number 6: getGIdemo
###################################################
ggi <- getGI("M22490")
ggi
gsq <- getSEQ("M22490")
substring(gsq,1,40)


###################################################
### chunk number 7: getGOdemo
###################################################
getSYMBOL("1000_at", "hgu95av2")
getGO("1000_at", "hgu95av2")[[1]][[1]]


###################################################
### chunk number 8: pmdemo
###################################################
hoxa9 <- "37809_at"
absts <- pm.getabst(hoxa9, "hgu95av2")
substring(abstText(absts[[1]][[1]]),1,60)


###################################################
### chunk number 9: GOcharacteristics
###################################################
 goTerms <- contents(GOTERM)
 goCats <- sapply(goTerms, Ontology)
 gCnums <- table(goCats)[c("BP","CC", "MF")]
 gCmat <- matrix(as.integer(gCnums), dimnames=list(c("BP","CC",
                                      "MF"), "Number of terms"))
 xtable.matrix(gCmat, display=c("d","d"),
  caption="Number of GO terms per ontology.", label="ta:GOprops")


###################################################
### chunk number 10: getMFchildren
###################################################
 get("GO:0008094", GOMFCHILDREN)


###################################################
### chunk number 11: getMFoffspring
###################################################
 get("GO:0008094", GOMFOFFSPRING)


###################################################
### chunk number 12: showEapply
###################################################
 hasChr <- eapply(GOTERM, function(x)
             x[grep("chromosome", Term(x))])
 lens <- sapply(hasChr, length)
 hasChr <- hasChr[lens>0]
 length(hasChr)


###################################################
### chunk number 13: lkGO6
###################################################
GOTerm2Tag <- function(term) {
 GTL <- eapply(GOTERM, function(x) {
     grep(term, x@Term, value=TRUE)})
 Gl <- sapply(GTL, length)
 names(GTL[Gl>0])
}


###################################################
### chunk number 14: lkGO7
###################################################
hasTFA <- GOTerm2Tag("transcription factor binding")
hasTFA


###################################################
### chunk number 15: demoLL
###################################################
gg1 <- get(GOTerm2Tag("^transcription factor binding$"), GOLOCUSID)
table(names(gg1))


###################################################
### chunk number 16: locusid
###################################################

 ll1 <- GOLOCUSID2GO[["7355"]]
 length(ll1)
 sapply(ll1, function(x) x$Ontology)



###################################################
### chunk number 17: getmappings
###################################################

getOntology(ll1, "BP")
getEvidence(ll1)
zz <- dropECode(ll1, code="IEA")
getEvidence(zz)



###################################################
### chunk number 18: loadKEGG
###################################################

library("KEGG")
library("xtable")

tmp1 <- as.list(KEGGPATHID2NAME)

pathnames <- ls(KEGGPATHID2EXTID)
species <- substr(pathnames, 1, 3)
paths <- substr(pathnames, 4, 20)




###################################################
### chunk number 19: spectable
###################################################
  y <- t(as.matrix(table(species)))
  dimnames(y)[[1]] <- "Counts"

  xtable(y, caption="Pathway counts per species",
          display=c("s", rep("d", ncol(y))), label="ta:speccounts")



###################################################
### chunk number 20: pathdemo
###################################################

 KEGGPATHID2NAME$"00140"

 KEGGPATHID2EXTID$hsa00140

 KEGGPATHID2EXTID$sce00140



###################################################
### chunk number 21: KEGGPW
###################################################

KEGGEXTID2PATHID$"5058"

KEGGEXTID2PATHID$"18479"



###################################################
### chunk number 22: assertionKEGG
###################################################
##we say this is length three in text
tmpA1 = KEGGEXTID2PATHID$"5058"
stopifnot(length(tmpA1)==3)



###################################################
### chunk number 23: cMAPInteractions
###################################################

keggproc <- eapply(cMAPKEGGINTERACTION, function(x) x$process)
table(unlist(keggproc))

cartaproc <- eapply(cMAPCARTAINTERACTION, function(x) x$process)
length(table(unlist(cartaproc)))



###################################################
### chunk number 24: assertion
###################################################
zz<-table(unlist(keggproc))
stopifnot( names(zz) == "reaction")



###################################################
### chunk number 25: preloadCM
###################################################
 cMK <- ls(cMAPKEGGPATHWAY)
 cMC <- ls(cMAPCARTAPATHWAY)


###################################################
### chunk number 26: cMAP
###################################################

 cMK <- ls(cMAPKEGGPATHWAY)
 spec <- substr(cMK, 1, 3)
 table(spec)

 cMK[[2]]
 pw2 <- cMAPKEGGPATHWAY[[cMK[2]]]
 names(pw2)
 pw2$name
 pw2$component


###################################################
### chunk number 27: lkint
###################################################
getI1 <- get("63", cMAPKEGGINTERACTION)
unlist(getI1[1:4])
unlist(getI1[[5]][[2]])


###################################################
### chunk number 28: lkint2
###################################################
get("2", cMAPKEGGMOLECULE)[[2]][7:8]


###################################################
### chunk number 29: loadhsa
###################################################
library("hsahomology")
ls("package:hsahomology")


###################################################
### chunk number 30: findESR1
###################################################
esrHG <- hsahomologyLL2HGID$"2099"
hesr <- get(as.character(esrHG), hsahomologyDATA)
sapply(hesr, function(x) x$homoOrg)


###################################################
### chunk number 31: findallH
###################################################

hXp <- eapply(hsahomologyDATA,  function(x) {
 gd <- sapply(x, function(x) if(!is.na(x$homoOrg) &&
             x$homoOrg == "xla") TRUE else FALSE)
  x[gd]})

 lh <- sapply(hXp, length)
 hXp2 <- hXp[lh>0]



