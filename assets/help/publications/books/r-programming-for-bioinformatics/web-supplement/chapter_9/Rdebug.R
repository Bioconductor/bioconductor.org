###################################################
### chunk: setup
###################################################
library("RBioinf")
require("methods")


###################################################
### chunk: setVNames with browser eval=FALSE
###################################################
## setVNames = function(x, nm)
## {
##     browser()
##     names(x) = nm
##     asSimpleVector(x, "numeric")
## }


###################################################
### chunk: call setVNames with browser
###################################################
x = 1:10
x = setVNames(x, letters[1:10])


###################################################
### chunk: 
###################################################
 foo=function(x, y) {
   x = 10
   z = 20
  baz(100)
 }

library("codetools")
 findGlobals(foo)


###################################################
### chunk: codet2
###################################################
findLocals(body(foo))


###################################################
### chunk: codet3
###################################################
 checkUsage(foo, name="foo", all=TRUE)


###################################################
### chunk: convertMode with mode list eval=FALSE
###################################################
##   x = convertMode(1:4, list())


###################################################
### chunk: traceback convertMode call eval=FALSE
###################################################
## traceback()


###################################################
### chunk: warning to error
###################################################
saveopt = options(warn=2)


###################################################
### chunk: warning to error
###################################################
options(saveopt)


###################################################
### chunk: callingHandlers eval=FALSE
###################################################
## withCallingHandlers(expression,
##     warning=function(c) recover())


###################################################
### chunk: set option error to recover
###################################################
 options(error=recover)


###################################################
### chunk: call convertMode with recover eval=FALSE
###################################################
## x = convertMode(1:4, list())


###################################################
### chunk: setVNames without browser eval=FALSE
###################################################
## rm("setVNames")


###################################################
### chunk: call to setVNames for a matrix
###################################################
x = matrix(1:4, nrow=2)
names(setVNames(x, letters[1:4]))


###################################################
### chunk: debug asSimpleVector
###################################################
debug(asSimpleVector)


###################################################
### chunk: call setVNames with debug eval=FALSE
###################################################
## names(setVNames(x, letters[1:4]))


###################################################
### chunk: undebug asSimpleVector
###################################################
undebug(asSimpleVector)


###################################################
### chunk: trace asSimpleVector
###################################################
trace(asSimpleVector)
x = list(1:3, 4:5)
for (i in seq(along=x)) {
    x[[i]] = asSimpleVector(x[[i]], "complex")
}
untrace(asSimpleVector)


###################################################
### chunk: asSimpleVectorShow eval=FALSE
###################################################
##  printWithNumbers(asSimpleVector)


###################################################
### chunk: call to trace for asSimpleVector
###################################################
trace(asSimpleVector, tracer=browser, at=9)


###################################################
### chunk: call asSimpleVector with trace eval=FALSE
###################################################
## names(setVNames(1:4, letters[1:4]))


###################################################
### chunk: call to untrace for asSimpleVector
###################################################
untrace(asSimpleVector)


###################################################
### chunk: subsetAsCharacter generic
###################################################
setGeneric("subsetAsCharacter")


###################################################
### chunk: subsetAsCharacter method
###################################################
setMethod("subsetAsCharacter",
          signature(x="character", i="missing",
                    j="missing"), function(x, i, j) x)


###################################################
### chunk: trace subsetAsCharacter method
###################################################
trace("subsetAsCharacter", tracer = browser,
    signature=c(x = "numeric"))


###################################################
### chunk: call to subsetAsCharacter that is traced eval=FALSE
###################################################
## subsetAsCharacter(1.5, 1:2)


###################################################
### chunk: calls to subsetAsCharacter that are not traced
###################################################
subsetAsCharacter(1+0i, 1:2)
subsetAsCharacter("x")
untrace("subsetAsCharacter")


###################################################
### chunk: profEx
###################################################
  Rprof()
  mad(runif(10000000))
  Rprof(NULL)


###################################################
### chunk: summaryRprof
###################################################
  summaryRprof()


###################################################
### chunk: showgc
###################################################
gc()


###################################################
### chunk: memprof
###################################################
ss = memory.profile()
sort(ss, decreasing=TRUE)
sum(ss)


###################################################
### chunk: RprofEx
###################################################

 library("affy")
 library("affydata")
 data(Dilution)
 Rprof(file="profRMA", memory.profiling = TRUE)

 r1 = rma(Dilution)

 Rprof(NULL)


###################################################
### chunk: RprofExSumm
###################################################

 pS = summaryRprof(file="profRMA", memory="tseries")
 names(pS)



###################################################
### chunk: plotTSprof
###################################################
plot(rownames(pS), pS$dup, type="l", xlab="Time", 
     ylab="Number of calls to duplicate")


###################################################
### chunk: rmaLM1
###################################################
  Rprofmem(file="rma2.out", threshold=100000)
  s2 = rma(Dilution)
  Rprofmem(NULL)  


###################################################
### chunk: rmaLM2Show eval=FALSE
###################################################
## noquote(readLines("rma2.out", n=5))


###################################################
### chunk: rmaLM3
###################################################
length(readLines("rma2.out"))


###################################################
### chunk: traceDil1
###################################################
 tracemem(Dilution)
 s3 <- rma(Dilution)
 tracemem(s3)


