###################################################
### chunk: storewd
###################################################
getwd()


###################################################
### chunk: simpleread
###################################################
fp = system.file("extdata/test1", package="RBioinf")
f1 = file(fp, open="r")
readLines(f1)
close(f1)


###################################################
### chunk: scanEx
###################################################
scan(fp, what="")


###################################################
### chunk: filepathEx
###################################################
file.path(R.home(), "doc")
system.file(package="RBioinf")


###################################################
### chunk: 
###################################################
cd = setwd(R.home())
list.files(path="doc")
list.files(pattern="Make")
list.files(pattern="Make", full.names=TRUE)


###################################################
### chunk: 
###################################################
setwd("doc")
getwd()
file.info("KEYWORDS")$isdir
file.info("manual")$isdir


###################################################
### chunk: 
###################################################
x = list.files()
x
files = x[!file.info(x)$isdir]
files
setwd(cd)


###################################################
### chunk: head
###################################################
head.file = function(x, n=6, ...) readLines(x, n)


###################################################
### chunk: 
###################################################
tempdir()
tmp1 = tempfile()
tmp2 = tempfile()
tmp1
tmp2


###################################################
### chunk: fileExs
###################################################
tmp1 = tempfile()
file.create(tmp1)
file.exists(tmp1)
file.access(tmp1, 2)
file.remove(tmp1)
file.exists(tmp1)


###################################################
### chunk: expand
###################################################
myhome = path.expand("~")
myhome
toR = file.path(myhome, "bin", "R")
toR


###################################################
### chunk: basename
###################################################
basename(toR)
dirname(toR)


###################################################
### chunk: filecreate
###################################################
z = file.path(tempdir(),"foo")
z
file.create(tmp1)
file.rename(tmp1,z)
file.exists(tmp1, z)
file.copy(z, tmp1)
file.exists(tmp1, z)
file.symlink(z, tmp2)
file.exists(tmp2)
fiz = file.info(z)
fitmp2 = file.info(tmp2)
all.equal(fiz, fitmp2)


###################################################
### chunk: dirandunlink
###################################################
newDir = file.path(tempdir(),"newDir")
newDir
newFile = file.path(newDir,"blah")
newFile
dir.create(newDir)
file.create(newFile)
file.exists(newDir, newFile)
unlink(newDir, recursive=TRUE)


###################################################
### chunk: showConnections
###################################################
showConnections(all=TRUE)


###################################################
### chunk: capabilities
###################################################
capabilities()


###################################################
### chunk: textConnection
###################################################

 zz = textConnection(LETTERS)
     readLines(zz, 2)
     showConnections(all=TRUE)
     scan(zz, "", 4)
     pushBack(c("aa", "bb"), zz)
     scan(zz, "", 4)
     close(zz)



###################################################
### chunk: textConnSink
###################################################
  savedOut = textConnection("foo", "w")
  sink(savedOut, split = TRUE)
  print(1:10)
  cat("That was my first command \n")
  letters[1:4]
  sink()
  close(savedOut)
  cat(foo, sep = "\n")


###################################################
### chunk: test1
###################################################
 p1 = pipe("cal 1 2006")
 p1
 readLines(p1)


###################################################
### chunk: Rcal
###################################################
library("RBioinf")
Rcal


###################################################
### chunk: pipe
###################################################
 ww = system("ls -1", intern=T)
 xx=readLines(pipe("ls -1"))

 all.equal(ww, xx)



###################################################
### chunk: testSeek
###################################################
 fName = file.path(tempdir(), "test1")
 file.create(fName)
 sFile = file(fName, open = "r+w")
 cat(1:10, sep="\n", file = sFile)
 seek(sFile)
 readLines(sFile, 3)
 seek(sFile, 2)
 readLines(sFile)
 close(sFile)
 unlink(fName)


###################################################
### chunk: 
###################################################
a = readLines(con=system.file("CONTENTS", package = "base"),n=2)
a
writeLines(a)


###################################################
### chunk: permIO
###################################################
mydir = tempdir()
for(i in 1:10) {
  fname = paste("perm", i, sep="")
  prm = sample(1:10, replace=FALSE)
  write(prm, file=file.path(mydir, fname))
}


###################################################
### chunk: catcallscal eval=FALSE
###################################################
## cat("10 2005", file="|cal")


###################################################
### chunk: DCF
###################################################
  x = read.dcf(file = system.file("CONTENTS", package = "base"),
                     fields = c("Entry", "Description"))
  head(x, n=3)
  write.dcf(x[1:3,],file="")


###################################################
### chunk: pform eval=FALSE
###################################################
## postForm("http://www.speakeasy.net/main.php", 
##          "some_text" = "Duncan", "choice" = "Ho", 
##          "radbut" = "eep", "box" = "box1, box2" )


###################################################
### chunk: PDBex
###################################################
library("RCurl")
url = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/"
fileNames = getURL(url, 
    .opts = list(customrequest = "NLST *.gz") )
fileNames = strsplit(fileNames, "\n", fixed=TRUE)[[1]]
fileNames = gsub("\r", "", fileNames)
length(fileNames)


###################################################
### chunk: download
###################################################
 fileNames[1]
 download.file(paste("ftp://ftp.wwpdb.org/pub/pdb/data/",
                     "structures/all/pdb/pdb100d.ent.gz", sep=""),
               destfile="pdb100d.ent.gz") 


