###################################################
### chunk: loadandExt
###################################################
library("hgu95av2")
chrVec = unlist(as.list(hgu95av2CHR))
table(chrVec)
class(chrVec)
#names(chrVec)[1:10]


###################################################
### chunk: splitCHR
###################################################
byChr = split(names(chrVec), chrVec)
sapply(byChr, length)


###################################################
### chunk: chrY
###################################################
byChr[["Y"]]


###################################################
### chunk: map
###################################################
library("hgu95av2")
hgu95av2MAP$"1001_at"


###################################################
### chunk: eApplyex
###################################################
myPos = eapply(hgu95av2MAP, function(x) grep("^17p", x, value=TRUE))
myPos = unlist(myPos)
length(myPos)


###################################################
### chunk: ans2
###################################################
myFindMap = function(mapEnv, which) {
  myg = ppc(which)
  a1 = eapply(mapEnv, function(x) grep(myg, x, value=TRUE))
  unlist(a1)
}


###################################################
### chunk: simpleCbind
###################################################
x = matrix(1:6, nc=2, dimnames=list(letters[1:3], LETTERS[1:2]))
y = matrix(21:26, nc=2, dimnames=list(letters[6:8], LETTERS[3:4]))
cbind(x,y)
rbind(x,y)


###################################################
### chunk: Stackex
###################################################
s1 = list(a=1:3,b= 11:12,c= letters[1:6])
ss = stack(s1)
ss
unsplit(s1, ss[,2])


###################################################
### chunk: chrMapEx
###################################################
mapP = as.list(hgu95av2MAP)
mLens = unlist(eapply(hgu95av2MAP, length))


###################################################
### chunk: summaryChr
###################################################
mlt = table(mLens)
mlt


###################################################
### chunk: 
###################################################
len3 = mLens[mLens==3]
hgu95av2SYMBOL[[names(len3)[1]]]
hgu95av2MAP[[names(len3)[1]]]


###################################################
### chunk: 
###################################################
len2 = names(mLens[mLens==2])
len2EG = unlist(mget(len2, hgu95av2ENTREZID))
len2EG = len2EG[!duplicated(len2EG)]
len2 = len2[!duplicated(len2EG)]
mapP = mget(len2, hgu95av2MAP)
hasX = sapply(mapP, function(x) if( length(grep("^X", x)) == 1) TRUE
else FALSE)

hasY = sapply(mapP, function(x) if( length(grep("^Y", x)) == 1) TRUE
else FALSE)
table(hasX & hasY)


###################################################
### chunk: missingMap
###################################################
missingMap = unlist(eapply(hgu95av2MAP, 
    function(x) any(is.na(x))))
table(missingMap)


###################################################
### chunk: chrMapEx2
###################################################
mapPs = sapply(mapP, function(x) x[1])
mapPs = mapPs[!is.na(mapPs)]

mapByPos = split(names(mapPs), mapPs)
table(sapply(mapByPos, length))


###################################################
### chunk: DBIEx
###################################################
library("RSQLite")
m = dbDriver("SQLite")
con = dbConnect(m, dbname="test")
data(USArrests)
dbWriteTable(con, "USArrests", USArrests, overwrite = TRUE)
dbListTables(con)


###################################################
### chunk: DBIresultSets
###################################################
rs = dbSendQuery(con, "select * from USArrests")
d1 = fetch(rs, n = 5)
d1
dbHasCompleted(rs)
dbListResults(con)
d2 = fetch(rs, n = -1)
dbHasCompleted(rs)
dbClearResult(rs)


###################################################
### chunk: DBIsimple
###################################################
dbListTables(con)
dbListFields(con, "USArrests")


###################################################
### chunk: SQLitegetalltables
###################################################
query = paste("SELECT name FROM sqlite_master WHERE",
    "type='table' ORDER BY name;")
rs = dbSendQuery(con, query)
fetch(rs, n= -1)


###################################################
### chunk: conditional selection
###################################################
rs = dbSendQuery(con, 
    "SELECT * FROM USArrests WHERE Murder > 10")


###################################################
### chunk: DBIcleanup
###################################################
unlink("test")


###################################################
### chunk: setupSQLite
###################################################
 library("RSQLite")
 m = dbDriver("SQLite")
 ##open up our test db    
 testDB = system.file("extdata/hgu95av2-sqlite.db", package="RBioinf")
 con = dbConnect(m, dbname = testDB)

 tabs = dbListTables(con)
 tabs 
 dbListFields(con, tabs[2])


###################################################
### chunk: SFex
###################################################
 query = paste("SELECT go_ont.go_id, go_ont.ont,", 
     "go_ont_name.ont_name FROM go_ont,",
     "go_ont_name WHERE (go_ont.ont = go_ont_name.ont)")
 rs = dbSendQuery(con, query)
 f3 = fetch(rs, n=3)
 f3
 dbClearResult(rs)


###################################################
### chunk: ijsol
###################################################
 query = paste("SELECT acc_num, go_id", 
     "FROM acc, go_probe", 
     "WHERE (acc.affy_id = go_probe.affy_id)")
 rs = dbSendQuery(con, query)
 f3 = fetch(rs, n=3)
 f3
 dbClearResult(rs)


###################################################
### chunk: SFex2
###################################################
query = paste("SELECT g1.*, g2.evi FROM go_probe g1,",
    "go_probe g2 WHERE  (g1.go_id = 'GO:0005737' ",
    "AND g2.go_id = 'GO:0005737') ",
    "AND (g1.affy_id = g2.affy_id) ",
    "AND (g1.evi = 'IDA' AND g2.evi = 'ISS')")
 rs = dbSendQuery(con, query)
 fetch(rs)


###################################################
### chunk: loadDBpkg
###################################################
library("hgu95av2.db")
mycon = hgu95av2_dbconn()


###################################################
### chunk: GOEX
###################################################
colnames(hgu95av2GO)
toTable(hgu95av2GO)[1:10,]
Lkeys(hgu95av2GO)[1:10]
Rkeys(hgu95av2GO)[1:10]


###################################################
### chunk: showLinks
###################################################
links(hgu95av2GO)[1:10,]


###################################################
### chunk: revMapEx
###################################################
is(hgu95av2SYMBOL, "Bimap")
rmMAP = revmap(hgu95av2SYMBOL)
rmMAP$"ABL1"


###################################################
### chunk: revmapList
###################################################
 myl=list(a="w", b="x", c="y")
 revmap(myl)


###################################################
### chunk: SYM2EG
###################################################
queryAlias = function(x) {
  it = paste("('", paste(x, collapse="', '"), "'", sep="")
  paste("select _id, alias_symbol from alias",
        "where alias_symbol in", it, ");")
}

queryGeneinfo = function(x) {
  it = paste("('", paste(x, collapse="', '"), "'", sep="")
  paste("select _id, symbol from gene_info where",
        "symbol in", it, ");")
}

queryGenes = function(x) {
  it = paste("('", paste(x, collapse="', '"), "'", sep="")
  paste("select * from genes where _id in", it,  ");")
}

findEGs = function(dbcon, symbols) {
  rs = dbSendQuery(dbcon, queryGeneinfo(symbols))
  a1 = fetch(rs, n=-1)
  stillLeft = setdiff(symbols, a1[,2])

  if( length(stillLeft)>0 ) {
    rs = dbSendQuery(dbcon, queryAlias(stillLeft))
    a2 = fetch(rs, n=-1)
    names(a2) = names(a1)
    a1 = rbind(a1, a2)
  } 
  
  rs = dbSendQuery(dbcon, queryGenes(a1[,1]))
  ans = merge(a1, fetch(rs, n=-1))
  dbClearResult(rs)
  ans
}


###################################################
### chunk: findEGs
###################################################
findEGs(mycon, c("ALL1", "AF4", "BCR", "ABL"))


###################################################
### chunk: combineGO_hgu95av2
###################################################
GOdbloc = system.file("extdata", "GO.sqlite", package="GO.db")
attachSql = paste("ATTACH '", GOdbloc, "' as go;", sep = "")
dbGetQuery(mycon, attachSql)


###################################################
### chunk: makeQuery
###################################################
sql = paste("SELECT DISTINCT a.go_id AS 'hgu95av2.go_id',",
            "a._id AS 'hgu95av2._id',",
            "g.go_id AS 'GO.go_id', g._id AS 'GO._id',",
            "g.ontology",
            "FROM go_bp_all AS a, go.go_term AS g", 
            "WHERE a.go_id = g.go_id LIMIT 10;")
dataOut = dbGetQuery(mycon, sql)
dataOut


###################################################
### chunk: dbschemaShow eval=FALSE
###################################################
## schema = capture.output(hgu95av2_dbschema())
## head(schema, 18)


###################################################
### chunk: QAlisting
###################################################
qcdata = capture.output(hgu95av2())
head(qcdata, 20)


###################################################
### chunk: mapcounts eval=FALSE
###################################################
## hgu95av2MAPCOUNTS
## hgu95av2_dbInfo()


###################################################
### chunk: getIntermediateDB
###################################################
tryCatch(library("human.db0"), error=function(e) {
      source("http://bioconductor.org/biocLite.R")
      biocLite("human.db0")
      library("human.db0") } )


###################################################
### chunk: gidb eval=FALSE
###################################################
##   source("http://bioconductor.org/biocLite.R")
##   biocLite("human.db0")


###################################################
### chunk: makeSqlite
###################################################
hgu95av2_IDs = system.file("extdata", 
                           "hgu95av2_ID", 
                           package="AnnotationDbi")

#Then specify some of the metadata for my database
myMeta = c("DBSCHEMA" = "HUMANCHIP_DB",
    "ORGANISM" = "Homo sapiens",
    "SPECIES" = "Human",
    "MANUFACTURER" = "Affymetrix",
    "CHIPNAME" = "Affymetrix Human Genome U95 Set Version 2",
    "MANUFACTURERURL" = "http:www.affymetrix.com")



###################################################
### chunk: td
###################################################
tmpout = tempdir()
popHUMANCHIPDB(affy = FALSE, prefix = "hgu95av2Test",
    fileName = hgu95av2_IDs, metaDataSrc = myMeta,
    baseMapType = "gb", outputDir = tmpout,
    printSchema = TRUE)


###################################################
### chunk: makeAnnDbPkg
###################################################
seed <- new("AnnDbPkgSeed",
            Package = "hgu95av2Test.db",
            Version = "1.0.0",
            PkgTemplate = "HUMANCHIP.DB",
            AnnObjPrefix = "hgu95av2Test")

makeAnnDbPkg(seed, 
             file.path(tmpout, "hgu95av2Test.sqlite"),
             dest_dir = tmpout)


###################################################
### chunk: SQLForge
###################################################
makeHUMANCHIP_DB(affy=FALSE,
    prefix="hgu95av2",
    fileName=hgu95av2_IDs,
    baseMapType="gb",
    outputDir = tmpout,
    version="2.1.0",
    manufacturer = "Affymetrix",
    chipName = "Affymetrix Human Genome U95 Set Version 2",
    manufacturerUrl = "http://www.affymetrix.com")


###################################################
### chunk: setupFileName
###################################################
Yeastfn = system.file("extdata", "yeast_small-01.xml", package="RBioinf")


###################################################
### chunk: checkNamespace
###################################################
yeastIntAct = xmlTreeParse(Yeastfn)
nsY = xmlNamespaceDefinitions(xmlRoot(yeastIntAct))
ns = getDefaultNamespace(xmlRoot(yeastIntAct))
namespaces = c(ns = ns)


###################################################
### chunk: numNspc
###################################################
length(nsY)
sapply(nsY, function(x) x[[2]])


###################################################
### chunk: DOM2
###################################################
nullf = function(x, ...) NULL
yeast2 = xmlTreeParse(Yeastfn, 
    handlers = list(sequence = nullf,
    organism = nullf, primaryRef = nullf, 
    secondaryRef = nullf,
    names = nullf), asTree=TRUE)


###################################################
### chunk: DOM1
###################################################
object.size(yeastIntAct)
object.size(yeast2)


###################################################
### chunk: DOM3
###################################################
yeast3 = xmlTreeParse(Yeastfn, useInternalNodes=TRUE)
f1 = getNodeSet(yeast3, "//ns:attributeList", namespaces)
length(f1)


###################################################
### chunk: getNodeSet
###################################################
iaM = getNodeSet(yeast3, 
    "//ns:interactionDetectionMethod//ns:fullName", 
    namespaces)
sapply(iaM, xmlValue)

f4 = getNodeSet(yeast3, "//ns:hostOrganism//ns:fullName", 
    namespaces)
sapply(f4, xmlValue)


###################################################
### chunk: interactors
###################################################
interactors = getNodeSet(yeast3, 
   "//ns:interactorList//ns:interactor",
   namespaces)
length(interactors)
interactions = getNodeSet(yeast3, 
   "//ns:interactionList/ns:interaction",
   namespaces)
length(interactions)


###################################################
### chunk: xpathApply
###################################################
 interactors = xpathApply(yeast3,
     "//ns:interactorList//ns:interactor", 
     xmlValue, namespaces = namespaces)


###################################################
### chunk: simpleEvent
###################################################
entSH = function(name, attrs, ...) {
          cat("Starting", name, "\n")
          level <<- attrs["level"]
          minorVersion <<- attrs["minorVersion"]
    }
e2 = new.env()
e2$level = NULL
e2$minorVersion = NULL
environment(entSH) = e2


###################################################
### chunk: twomoreHandlers
###################################################
hOrg = function(name, attrs, ...) {
         taxid <<- c(attrs["ncbiTaxId"], taxid)
		      }
e3 = new.env()
e3$taxid = NULL
environment(hOrg) = e3

hInt = function(name, attrs, ...)
		     numInt <<- numInt + 1

e3$numInt = 0
environment(hInt) = e3


###################################################
### chunk: xmlEventP
###################################################
s1 = xmlEventParse(Yeastfn,
	handlers = list(entrySet = entSH, hostOrganism = hOrg,
		         interactor = hInt))

environment(s1$entrySet)$level
environment(s1$hostOrganism)$taxid
environment(s1$interactor)$numInt


###################################################
### chunk: HTMLparsing
###################################################
url = paste("http://www.bioconductor.org/checkResults/", 
    "2.1/bioc-LATEST/", sep="")
s1 = htmlTreeParse(url, useInternalNodes=TRUE)
class(s1)


###################################################
### chunk: FindNodes
###################################################
f1 = getNodeSet(s1, "//a[@href]")
length(f1)


###################################################
### chunk: RefinedGetNodes
###################################################
f2 = getNodeSet(s1, "//b/a[@href]")
p2 = sapply(f2, xmlValue)
length(p2)
p2[1:10]


###################################################
### chunk: 
###################################################
pkgs = sapply(f1, xmlGetAttr, "href")
pkg = grep("/packages/2.1/bioc/html/", pkgs, fixed=TRUE)



###################################################
### chunk: Entrezex
###################################################
 ezURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
 t1 = url(ezURL, open="r")
 if( isOpen(t1) ) {
    z = xmlTreeParse(paste(ezURL, "einfo.fcgi", sep=""),
     isURL=TRUE, handlers=NULL, asTree=TRUE) 

    dbL = xmlChildren(z[[1]]$children$eInfoResult)$DbList

    dbNames = xmlSApply(dbL, xmlValue)
    
    length(dbNames)

    dbNames[1:5]
} 


###################################################
### chunk: listMarts
###################################################
library("biomaRt")
head(listMarts())


###################################################
### chunk: selectMart
###################################################
ensM = useMart("ensembl")

ensData = head(listDatasets(ensM))
dim(ensData)

ensMH = useDataset("hsapiens_gene_ensembl", mart=ensM)


###################################################
### chunk: selM eval=FALSE
###################################################
## ensMH = useMart("ensembl", 
##     dataset = "hsapiens_gene_ensembl")


###################################################
### chunk: filters
###################################################
filterSummary(ensMH)
lfilt = listFilters(ensMH, group="GENE:")
nrow(lfilt)
head(lfilt)


###################################################
### chunk: attributes
###################################################
head(attributeSummary(ensMH))
lattr = listAttributes(ensMH, group="PROTEIN:")
lattr


###################################################
### chunk: bioMIntro
###################################################
entrezID = c("983", "3581", "1017") 
rval = getGene(id=entrezID, type="entrezgene", mart = ensMH) 
unique(rval$hgnc_symbol) 


###################################################
### chunk: INTERPRO
###################################################
ensembl = useMart("ensembl", 
    dataset = "hsapiens_gene_ensembl")
    ipro = getBM(attributes=c("entrezgene","interpro",
    "interpro_description"), 
filters = "entrezgene", values = entrezID, 
    mart = ensembl) 
ipro


###################################################
### chunk: GEOqueryEx
###################################################

 library(GEOquery) 

 gds = getGEO("GDS10") 

 eset = GDS2eSet(gds, do.log2 = TRUE) 



###################################################
### chunk: GEOqueryPrintExSetShow eval=FALSE
###################################################
## s1 = experimentData(eset)
## abstract(s1)
## s1@pubMedIds


###################################################
### chunk: KEGGEx
###################################################
library("KEGG")
library("KEGGSOAP")
KEGGPATHID2NAME$"00740"
SoapAns = get.genes.by.pathway("path:sce00740")
SoapAns


###################################################
### chunk: LocalKEGG
###################################################
SA = gsub("^sce:", "", SoapAns)
localAns = KEGGPATHID2EXTID$"sce00740"
setdiff(SA, localAns)


