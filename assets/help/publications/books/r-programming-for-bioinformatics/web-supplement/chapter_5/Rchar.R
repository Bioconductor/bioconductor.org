###################################################
### chunk: excharv
###################################################
 mychar = c("as", "soon", "as possible")
 mychar
 nchar(mychar)


###################################################
### chunk: zerolenChar
###################################################
 x=character(0) #length 0 character vector
 length(x)
 nchar(x)
 y = ""  #length 0 string
 length(y)
 nchar(y)


###################################################
### chunk: substr
###################################################
substr(x, 2, 4)
substr(x, 2, rep(4,5))
substring(x, 2, rep(4,5))


###################################################
### chunk: DNA2AA
###################################################
 rD = randDNA(102)
 rDtriples = substring(rD, seq(1, 102, by=3), seq(3, 102, 3))
 paste(GENETIC_CODE[rDtriples])


###################################################
### chunk: DNA2AA2
###################################################
 DNA2AA = function(DNAseq) {
    nc = nchar(DNAseq)
    Dtriples = substring(DNAseq, seq(1, nc, by=3), seq(3, nc, 3))
    paste(GENETIC_CODE[Dtriples], collapse="")
 }


###################################################
### chunk: substrgets
###################################################
substring(x, 2, 4) = "abc"
x
x=c("howdy", "dudey friend")
substr(x, 2, 4) = "def"
x
substring(x, 2) <- c("..", "+++")


###################################################
### chunk: paste1
###################################################
paste(1:3, "+", 4:5)
paste(1:3, 1:3, 4:6, sep="+")


###################################################
### chunk: pasteEx
###################################################
paste(1:4, collapse="=")


###################################################
### chunk: strsplitEx
###################################################
strsplit(c("ab", "cde", "XYZ"), c("Y", ""))
strsplit(c("ab", "cde", "XYZ"), c("Y", NULL))


###################################################
### chunk: strtrim
###################################################
 x <- paste(readLines(file.path(R.home(), "COPYING")), collapse = "\n")
strwrap(x, 30, prefix="myidea: ")[1:10]
writeLines(strwrap(x, 30, prefix="myidea: ")[1:5])


###################################################
### chunk: toupper
###################################################
dna2rna = function(inputStr) {
   if(!is.character(inputStr))
      stop("need character input")
   is=toupper(inputStr)
   chartr("T", "U", is)
}

x=c(randDNA(15), randDNA(12))
x
dna2rna(x)


###################################################
### chunk: compSeq
###################################################
compSeq = function(x) 
   chartr("ACTG", "TGAC", x)
compSeq(x)


###################################################
### chunk: isDNA
###################################################
 isDNA = function(x) {
     xU = toupper(x)
     spx = strsplit(xU, NULL)
     sapply( spx, function(z) all( z %in% c("A","C", "G", "T")))
 }


###################################################
### chunk: revString
###################################################
strReverse <- function(x)
             sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")


###################################################
### chunk: revComp
###################################################
reverseComplement = function(x) strReverse(compSeq(x))


###################################################
### chunk: charComps
###################################################
 set.seed(123)
 x = sample(letters[1:10], 5)
 x
 sort(x)
 x < "m"


###################################################
### chunk: matchEx
###################################################
exT = c("Intron", "Exon", "Example", "Chromosome")

match("Exon", exT)

"Example" %in% exT


###################################################
### chunk: pmatchEx
###################################################
pmatch("E", exT)
pmatch("Exo", exT)
pmatch("I", exT)
charmatch("I", exT)


###################################################
### chunk: pmEx2
###################################################
pmatch(c("I", "Int"), exT)
pmatch(c("I", "Int"), exT, duplicates.ok=TRUE)
charmatch(c("I", "Int"), exT)


###################################################
### chunk: pmEx3
###################################################
pmatch(c("ab"), c("ab", "ab"))
charmatch(c("ab"), c("ab", "ab"))


###################################################
### chunk: sw1
###################################################
'I\'m a string'
  "I'm a string"


###################################################
### chunk: ss2
###################################################
 s = "I'm a backslash: \\"
 s


###################################################
### chunk: ss3
###################################################
 cat(s)


###################################################
### chunk: noquote
###################################################
noquote(s)


###################################################
### chunk: ss3
###################################################
 nchar(s)
 strsplit(s, NULL)[[1]]
 nchar("\n")
 charToRaw("\n")


###################################################
### chunk: translate
###################################################
 fn = "c:\\My Documents\\foo.bar"
 fn


###################################################
### chunk: trans2
###################################################
 old = "\\"
 new = "/"
 chartr(old, new, fn)


###################################################
### chunk: trans3
###################################################
 gsub("\\\\", new, fn)
 gsub("\\", new, fn, fixed=TRUE) 


###################################################
### chunk: parse
###################################################
 v1 = parse(text="mean(1:10)")
 v1
 eval(v1)
 deparse(v1)
 deparse(v1[[1]])


###################################################
### chunk: localeShow eval=FALSE
###################################################
##  Sys.getlocale()


###################################################
### chunk: ansX
###################################################
onlyNuc = function(x) {
    ans = rep(TRUE, length(x))
    noNuc = grep("[^ACTGactg]", x)
    ans[noNuc] = FALSE
    ans
}


###################################################
### chunk: Anchor
###################################################
gregexpr("\\<", "my first anchor")
gregexpr("\\>", "my first anchor")


###################################################
### chunk: simpleEx
###################################################
gregexpr("\\b", "once upon a time")
gregexpr("\\>", "once upon a time")
gregexpr("\\<", "once upon a time")
gregexpr("^", "once upon a time")
gregexpr("$", "once upon a time")


###################################################
### chunk: lookahead
###################################################
regexpr("r[^r]", "asffrb", perl=TRUE)
regexpr("r[^r]", "asffr", perl=TRUE)
regexpr("r(?!r)", "asffrb", perl=TRUE)
regexpr("r(?!r)", "asffr", perl=TRUE)


###################################################
### chunk: backrefs
###################################################
gregexpr("([A-Z])\\1", "ABBBZZ")



###################################################
### chunk: greedyRep
###################################################
regexpr("AB{2,4}?", "ABBBBB")
regexpr("AB{2,4}?", "ABBBBB", perl=T)


###################################################
### chunk: Matching1
###################################################
regexpr("foo|foobar", "myfoobar")
regexpr("foo|foobar", "myfoobar", perl=TRUE)


###################################################
### chunk: Matching2
###################################################
testS = "ACTACCACTACCACT"
gregexpr("ACTACCACT", testS)
gregexpr2("ACTACCACT", testS)


###################################################
### chunk: MMDDYYYY
###################################################
regexpr("\\d\\d\\/\\d\\d\\/\\d\\d\\d\\d", 
        "today is 12/01/1977", perl = TRUE)
regexpr("\\d\\d\\/\\d\\d\\/\\d\\d\\d\\d", 
        "today is 21/41/1977", perl = TRUE)


###################################################
### chunk: MMDDYYsol
###################################################
regexpr("(0[1-9]|1[0-2])\\/\\d\\d\\/\\d\\d\\d\\d", 
        "today is 12/01/1977", perl=TRUE)

regexpr("(0[1-9]|1[0-2])\\/\\d\\d\\/\\d\\d\\d\\d", 
        "today is 21/01/1977", perl=TRUE)


###################################################
### chunk: stripWhite
###################################################
strwhite = function(x, lead=TRUE, trail=TRUE) {
  if(lead) 
   x = sub("^[[:blank:]]*", "", x, perl=TRUE)
  if(trail)
   sub("[[:blank:]]*$", "", x, perl=TRUE)
  else 
   x
}


###################################################
### chunk: prositeRE
###################################################
prositeM = "[RK]-x(2,3)-[DE]-x(2,3)-Y."
regexM = gsub("-|\\.", "", prositeM)
regexM = chartr("xX()", "..{}", regexM)


###################################################
### chunk: nowTest
###################################################
testP = "ACRDRACDTUYACRD"
testN = "ACRDRAXXCDTUYACRD"
regexpr(regexM, testP)
regexpr(regexM, testN)


###################################################
### chunk: lcPreSuf
###################################################
library("Biobase")

str1 = c("not now", "not as hard as wow", "not something new")
lcPrefix(str1)
lcSuffix(str1)




###################################################
### chunk: suffT
###################################################
library("Rlibstree")
s1 = "biology"
getLongestSubstring(s1)


###################################################
### chunk: RNAEx
###################################################
st1 = RNAString("UCUCCCAACCCUUGUACCAGUG")
cD = cDNA(st1)


###################################################
### chunk: matchprobes
###################################################
  library("matchprobes")
  seq <- c("CGACTGAGACCAAGACCTACAACAG",
           "CCCGCATCATCTTTCCTGTGCTCTT")

 complementSeq(seq, start=13, stop=13)


###################################################
### chunk: compS
###################################################
cS = function(x, DNA=TRUE) {
    if(DNA) 
        chartr("GCAT", "CGTA", x)
     else
        chartr("GCAU", "CGUA", x)
}


###################################################
### chunk: BSGenomeShow eval=FALSE
###################################################
## library("BSgenome.Hsapiens.UCSC.hg18")
## Hsapiens


###################################################
### chunk: 
###################################################
chr22NoN <- mask(Hsapiens$chr22, "N")
alphabetFrequency(Hsapiens$chr22, freq=TRUE)["N"]


###################################################
### chunk: findTATA
###################################################
TATA="TATAAAA"
mT = matchPattern(TATA, chr22NoN)
countPattern(TATA, chr22NoN)


###################################################
### chunk: fTT2
###################################################
mmT = matchPattern(TATA, chr22NoN, max.mismatch=1)
length(mmT)
mismatch(TATA, mmT[1:3])


###################################################
### chunk: mPDict
###################################################
 library(hgu95av2probe)
 dict <- hgu95av2probe$sequence 
 length(dict)                   
 unique(nchar(dict))            
 dict[1:5]
 pdict <- PDict(dict) 
 
 vindex <- matchPDict(pdict, Hsapiens$chr22) 

 length(vindex)                        
 count_index <- countIndex(vindex)
 sum(count_index)
 table(count_index)


###################################################
### chunk: confirm
###################################################
   dict[count_index == max(count_index)]
   countPattern("CTGTAATCCCAGCACTTTGGGAGGC", Hsapiens$chr22)


###################################################
### chunk: findPalindromes
###################################################
 chr22_pals = findPalindromes(chr22NoN, min.armlength = 40,
                       max.looplength = 20)
 nchar(chr22_pals)
 palindromeArmLength(chr22_pals) 
 palindromeLeftArm(chr22_pals)


 ans = alphabetFrequency(chr22_pals,
       base = TRUE)

head(ans, n = 15)


###################################################
### chunk: ans
###################################################
allThere = function(x) all(x[1:4] > 0)
x1 = apply(ans, 1, allThere)
chr22_pals[x1]


###################################################
### chunk: compPal
###################################################

cpals = findComplementedPalindromes(chr22NoN, min.armlength=40,
                       max.looplength = 20)
cpals


###################################################
### chunk: 
###################################################
  Lpattern <- "CTCCGAG"
  Rpattern <- "GTTCACA"
  LRans = matchLRPatterns(Lpattern, Rpattern, 500, Hsapiens$chr22)
  length(LRans)


###################################################
### chunk: align1
###################################################
   
  aa1 <- AAString("HXBLVYMGCHFDCXVBEHIKQZ")
  aa2 <- AAString("QRNYMYCFQCISGNEYKQN")
  needwunsQS(aa1, aa2, "BLOSUM62", gappen=3)

  ## See how the gap penalty influences the alignment
  needwunsQS(aa1, aa2, "BLOSUM62", gappen=8)


###################################################
### chunk: DNAalign
###################################################
oldD = setwd(system.file("extdata", package="Biostrings"))
Sc = readFASTA("Sc.fa", )[[1]]$seq
Sp = readFASTA("Sp.fa")[[1]]$seq
setwd(oldD)
mat <- matrix(-5L, nrow=4, ncol=4)
for (i in seq_len(4)) mat[i, i] <- 0L
rownames(mat) <- colnames(mat) <- DNA_ALPHABET[1:4]
dnaAlign1 = needwunsQS(Sc, Sp, mat, gappen=1)
nchar(dnaAlign1)


###################################################
### chunk: align2
###################################################
dnaAlign2 = needwunsQS(Sc, Sp, mat, gappen=6)
nchar(dnaAlign2)


###################################################
### chunk: modifyPenalty
###################################################
mat2 = mat
mat2["C", "T"] = mat2["T", "C"] = -1L
dnaAlign3 = needwunsQS(Sc, Sp, mat2, gappen=1)
dnaAlign4 = needwunsQS(Sc, Sp, mat2, gappen=6)
nchar(dnaAlign4)


###################################################
### chunk: longestSharedString
###################################################
library("Rlibstree")
tree = SuffixTree(c(Sc, Sp))
MUM = getLongestCommonSubstring(tree)
nchar(MUM)


###################################################
### chunk: consmat
###################################################
consmat(dnaAlign1)[, 1:20]


