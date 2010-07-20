###################################################
### chunk: loadLibs
###################################################
library("RBioinf")
library("graph")
library("Rgraphviz")
library("methods")


###################################################
### chunk: ffex
###################################################
setClass("Passenger", representation(name="character",
                                     origin="character", 
                                     destination="character"))
setClass("FreqFlyer", representation(ffnumber = "numeric"),
   contains = "Passenger")
getClass("FreqFlyer")
subClassNames("Passenger")
superClassNames("FreqFlyer")


###################################################
### chunk: S4printff
###################################################
setMethod("show", "Passenger",
          function(object) {
              cat("Name: ", object@name, "\n")
              cat("Origin: ", object@origin, "\n")
              cat("Destination:", object@destination, "\n")
          })

p1 = new("Passenger", name="J. Biologist", origin="YXY", destination="TGL")
p1


###################################################
### chunk: showFF
###################################################
setMethod("show", "FreqFlyer",
          function(object) {
              callNextMethod()
              cat("Freq Flyer no: ", object@ffnumber, "\n") 
          })

p2 = new("FreqFlyer", name = "K. Biologist", 
        origin = "YVR", destination="LAX", ffnumber = 1)
p2


###################################################
### chunk: rectEx
###################################################
setClass("Rectangle", 
   representation(h="numeric", w="numeric", area="numeric"))

myr = new("Rectangle", h=10, w=20, area=200)

setGeneric("area", function(shape) standardGeneric("area"))

setMethod("area", signature(shape = "Rectangle"), function(shape) shape@area)

myr@area

area(myr)


###################################################
### chunk: changeClasses
###################################################

setClass("Rectangle", representation(h="numeric", w="numeric"))

setMethod("area", "Rectangle", function(shape) shape@h * shape@w)


###################################################
### chunk: oddS3
###################################################
 x = 1:10
 class(x)
 dim(x) = c(2,5)
 class(x)
 attr(x, "class")
 inherits(x, "integer")


###################################################
### chunk: S3ex1
###################################################
 x=list(name="Josephine Biologist",
    origin = "SEA", destination = "YXY")

 class(x) = "Passenger"

 y = list(name="Josephine Physicist",
    origin = "SEA", destination = "YVR", ffnumber = 10)

 class(y) = c("FreqFlyer", "Passenger")

 inherits(x, "Passenger")
 inherits(x, "FreqFlyer")
 inherits(y, "Passenger")


###################################################
### chunk: isobject
###################################################
x = 1:10
is.object(x)
class(x) = "myint"
is.object(x)


###################################################
### chunk: showClassSetting
###################################################
x = matrix(1:10, nc=2)
class(x) = "matrix"
x
is.object(x)
oldClass(x) = "matrix"
x
is.object(x)


###################################################
### chunk: glmEX
###################################################
   counts = c(18,17,15,20,10,20,25,13,12)
   outcome = gl(3,1,9)
   treatment = gl(3,3)
   d.AD = data.frame(treatment, outcome, counts)
   glm.D93 = glm(counts ~ outcome + treatment, family=poisson())
   class(glm.D93)
   is.list(glm.D93)
   attr(glm.D93, "class")


###################################################
### chunk: glmEX2
###################################################
names(glm.D93)


###################################################
### chunk: VARLS3
###################################################

 ex1VL = c("Sex, M=MALE, F=FEMALE", "Age in years")
 names(ex1VL) = c("Sex", "Age")
 class(ex1VL) = "VARLS3"


###################################################
### chunk: simData
###################################################
 set.seed(123)
 simExprs = matrix(rgamma(10000, 500), nc=10, nr=100)
 simS = sample(c("M", "F"), 10, rep=TRUE)
 simA = sample(30:45, 10, rep=TRUE)
 simPD = data.frame(Sex=simS, Age=simA)


###################################################
### chunk: newS3PD
###################################################

 new.EXPRS3 = function(Class, eData, pData, cDesc) {
     if(!is.matrix(eData) )
         stop("invalid expression data")
     if(!is.data.frame(pData) )
         stop("invalid phenotypic data")
     if(!inherits(cDesc, "VARLS3"))
         stop("invalid cov description")
     ncE = ncol(eData)
     nrP = nrow(pData)
     if( ncE != nrP )
         stop("incorrect dimensions")
     pD = list(pData=pData, varLabels=cDesc)
     class(pD) = "PHENODS3"
     ans = list(exprs=eData, phenoData = pD)
     class(ans) = class(Class)
     ans
}


###################################################
### chunk: makeNewS3
###################################################
 myES3 = new.EXPRS3("EXPRS3", simExprs, simPD, ex1VL)


###################################################
### chunk: simpleGeneric
###################################################
fun = function(x, ...) UseMethod("fun")
fun.default = function(x, ...) print("In the default method")
fun(2)


###################################################
### chunk: showNextMethod
###################################################
fun.Foo = function(x) {
    print("start of fun.Foo")
    NextMethod()
    print("end of fun.Foo")
}
fun.Bar = function(x) {
    print("start of fun.Bar")
    NextMethod()
    print("end of fun.Bar")
}


###################################################
### chunk: shwNM
###################################################
x = 1
class(x) = c("Foo", "Bar")
fun(x)


###################################################
### chunk: ES3printmethod
###################################################
  print.PHENODS3 = function(object) {
      dm = dim(object$pData)
      cat("instance of PHENODS3 with", dm[2], "variables")
      cat("and", dm[1], "cases\n")
      vL = object$varLabels
      cat("\t varLabels\n")
      nm = names(vL)
      for(i in seq(along=vL) )
         cat("\t\t", nm[[i]], ": ", vL[[i]], "\n", sep="")
   }

  print.EXPRS3 = function(object) {
      dm = dim(object$exprs)
      cat("instance of EXPRS3\n")
      cat("number of genes:", dm[1], "\n")
      cat("number of samples:", dm[2], "\n")
      print(object$phenoData)
   }


###################################################
### chunk: S3methods
###################################################
methods("mean")


###################################################
### chunk: S3glmmethods
###################################################
 methods(class="glm")


###################################################
### chunk: S3genericEx
###################################################
 fun.Foo = function(x, ...) print(ls(all=TRUE))

 y=1
 class(y) = c("Foo", "Zip", "Zoom")
 fun(y)


###################################################
### chunk: S3generic2
###################################################
 fun.Foo = function(x, ...) {
  print(paste(".Generic =", .Generic))
  print(paste(".Class =", paste(.Class, collapse=", ")))
  print(paste(".Method =", .Method))
}
 fun(y)


###################################################
### chunk: assmet
###################################################
methods("$<-")


###################################################
### chunk: matrixRepMethod
###################################################
 "$<-.matrix" = function(x, name, value) {
      if( ! name %in% row.names(x)) stop("bad name")
      x[name,] = value
      x
  }


###################################################
### chunk: reptest
###################################################
x = matrix(1:10, nr=5)
rownames(x) = letters[1:5]
oldClass(x) = "matrix"

x$c = c(100, 100)
x


###################################################
### chunk: setClass
###################################################
  setClass("A", representation(s1="numeric"), 
               prototype=prototype(s1=0))
  myA = new("A")
  myA

  m2 = new("A", s1=10)
  m2


###################################################
### chunk: superClass
###################################################
   setClass("B", contains="A", representation(s2="character"),
          prototype=list(s2="hi"))
   myB = new("B")
   myB


###################################################
### chunk: removeClass
###################################################
setClass("Ohno", representation(y="numeric"))
getClass("Ohno")
removeClass("Ohno")
tryCatch(getClass("Ohno"), error=function(x) "Ohno is gone")


###################################################
### chunk: aSlots
###################################################
getSlots("A")
slotNames("A")


###################################################
### chunk: showExtends
###################################################
extends("B")
extends("B", "A")
extends("A", "B")
superClassNames("B")
subClassNames("A")


###################################################
### chunk: simpleUseShow1 eval=FALSE
###################################################
## getClass("matrix")


###################################################
### chunk: simpleUseShow2
###################################################
extends("matrix")


###################################################
### chunk: coercion
###################################################
myb = new("B")
as(myb, "A")


###################################################
### chunk: replacementcoerce
###################################################
mya = new("A", s1 = 20)
as(myb, "A") <- mya
myb


###################################################
### chunk: am2mat eval=FALSE
###################################################
## setAs(from="graphAM", to="matrix",
##       function(from) {
##           if ("weight" %in% names(edgeDataDefaults(from))) {
##               tm <- t(from@adjMat)
##               tm[tm != 0] <- unlist(edgeData(from, attr="weight"))
##               m <- t(tm)
##           } else {
##               m <- from@adjMat         
##           }
##           rownames(m) <- colnames(m)
##           m
##       })


###################################################
### chunk: initialize
###################################################
  setClass("Ex1", representation(s1="numeric"),
           prototype=prototype(s1=rnorm(10)))
  b = new("Ex1")
  b


###################################################
### chunk: ex2
###################################################
 b2 = new("Ex1")


###################################################
### chunk: constructorFun
###################################################
makeex = function() {
    obj = new("Ex1")
    obj@s1 = rnorm(10)
    }
 b2 = makeex()


###################################################
### chunk: protoInherits
###################################################
bb = getClass("B")
bb@prototype


###################################################
### chunk: initializeEx
###################################################
setClass("W", representation(c1 = "character"))
setClass("WA", contains=(c("A", "W")))
a1 = new("A", s1=20)
w1 = new("W", c1 = "hi")
new("WA", a1, w1)


###################################################
### chunk: init2
###################################################
setClass("XX", representation(a1 = "numeric", 
    b1 = "character"), 
    prototype(a1 = 8, b1 = "hi there"))

new("XX")

setMethod("initialize", "XX", function(.Object, ..., b1) {
   callNextMethod(.Object, ..., b1 = b1, a1 = nchar(b1))
})

 new("XX", b1="yowser")


###################################################
### chunk: MMEX
###################################################
 setClass("Capital",
          representation=representation(
            string="character"))
 setClass("CountedCapital",
          contains="Capital",
          representation=representation(
            length="numeric"))
 setMethod("initialize",
           "Capital", 
           function(.Object, ..., string=character(0)) {
             string <- toupper(string)
             callNextMethod(.Object, ..., string=string)
           })
 setMethod("initialize",
           "CountedCapital",
           function(.Object, ...) {
             .Object <- callNextMethod()
             .Object@length <- nchar(.Object@string)
             .Object
           })
 new("Capital", string="MiXeD")
 new("CountedCapital", string="MiXeD")
 new("CountedCapital", string=c("MiXeD", "MeSsaGe"))


###################################################
### chunk: noslots
###################################################
setClass("seq", contains="numeric", 
         prototype=prototype(numeric(3)))
s1 = new("seq")
s1
slotNames(s1)


###################################################
### chunk: initNoSlot
###################################################
 setMethod("initialize", "seq", function(.Object) {
   .Object[1]=10; .Object})

 new("seq")


###################################################
### chunk: showExt
###################################################

tryCatch(setMethod("[", signature("integer"), 
                   function(x, i, j, drop) print("howdy")), 
         error = function(e)
         print("we failed"))

setClass("Myint", representation("integer"))
setMethod("[", signature("Myint"), 
                   function(x, i, j, drop) print("howdy"))
x = new("Myint", 4:5)
x[3]


###################################################
### chunk: funExt
###################################################
setClass("DBFunc", "function")
setMethod("$", signature = c("DBFunc", "character"),
  function(x, name) x(name))


###################################################
### chunk: funExt2
###################################################
mytestFun = function(arg) print(arg)

mtF = new("DBFunc", mytestFun)
mtF$y
is(mtF, "function")


###################################################
### chunk: S4attr
###################################################
mya = new("A", s1 = 20)
class(mya)
attributes(mya)


###################################################
### chunk: S4attr2
###################################################
attr(mya, "s1") <- "L" 
mya


###################################################
### chunk: cUex
###################################################
setClassUnion("lorN", c("list", "NULL"))
subClassNames("lorN")
superClassNames("lorN")

isVirtualClass("lorN")
isClassUnion("lorN")


###################################################
### chunk: access
###################################################
setClass("Foo", representation(a="ANY"))
setGeneric("a", function(object) standardGeneric("a"))
setMethod("a", "Foo", function(object) object@a)
b = new("Foo", a=10)
a(b)


###################################################
### chunk: setOldClass
###################################################
setOldClass("mymatrix")
getClass("mymatrix")


###################################################
### chunk: useClass
###################################################
setClass("myS4mat", representation(m = "mymatrix"))
x=matrix(1:10, nc=2)
class(x) = "mymatrix"
m4 = new("myS4mat", m=x)


###################################################
### chunk: findalloldclasses
###################################################
head(subClassNames(getClass("oldClass")))


###################################################
### chunk: simpleGeneric
###################################################
 setGeneric("foo", 
   function(object, x) standardGeneric("foo") )

 setMethod("foo", signature("numeric", "character"),
   function(object, x) print("Hi, I'm method one"))


###################################################
### chunk: genSig
###################################################
  setGeneric("genSig", signature=c("x"), 
             function(x, y=1) standardGeneric("genSig"))
   
  setMethod("genSig", signature("numeric"), 
            function(x, y=20) print(y))

  genSig(10)


###################################################
### chunk: genericReturns
###################################################
 setGeneric("foo", function(x,y,...) {
  y = standardGeneric("foo")
   print("I'm back")
  y
 })

 setMethod("foo", "numeric", function(x,y,...) {print("I'm gone")}
 )

 foo(1)


###################################################
### chunk: getGenericsEx
###################################################
library("Biobase")
allG = getGenerics()
allGs = split(allG@.Data, allG@package)
allGBB = allGs[["Biobase"]]
length(allGBB)


###################################################
### chunk: bypkg
###################################################
allGbb = getGenerics("package:Biobase")
length(allGbb)


###################################################
### chunk: Rdots
###################################################
 setGeneric("bar", function(x, y, ...) standardGeneric("bar"))

 setMethod("bar", signature("numeric", "numeric"),
    function(x, y, d) print("Method1"))

 ##removes the method above
 setMethod("bar", signature("numeric", "numeric"),
    function(x, y, z) print("Method2"))

 bar(1,1,z=20)
 bar(2,2,30)
 tryCatch(bar(2,4,d=20), error=function(e) 
          print("no method1"))


###################################################
### chunk: replacement
###################################################

 setGeneric("a<-", function(x, value)
            standardGeneric("a<-"))

 setReplaceMethod("a", "Foo",
  function(x, value) {
    x@a = value
    x
  })

  a(b) = 32

  b


###################################################
### chunk: BiobEx eval=FALSE
###################################################
## setMethod("$", "eSet", function(x, name) {
##     eval(substitute(phenoData(x)$NAME_ARG, 
##                     list(NAME_ARG = name)))
## })


###################################################
### chunk: ldots
###################################################
 cnew = function(x, ...) {
   if(nargs() < 3)
    c2(x, ...)
   else
    c2(x, cnew(...))
 }


###################################################
### chunk: ldots2
###################################################
setGeneric("c2", function(x, y) standardGeneric("c2"))


###################################################
### chunk: plus
###################################################
setMethod("c2", signature("numeric", "numeric"), function(x, y) x + y) 
cnew(1,2,3,4)


###################################################
### chunk: isvsinherits
###################################################
x = 1
class(x) = c("C1", "C2")
is(x, "C2")
inherits(x, "C2")


###################################################
### chunk: isworks
###################################################
setClass("A")
setClass("B", representation(s="numeric"), contains = "A")
y = new("B")
is(y, "A")
inherits(y, "A")


###################################################
### chunk: setOldC
###################################################
setOldClass(c("C1", "C2"))
is(x, "C2")


###################################################
### chunk: showasS4
###################################################
x = 1
setClass("A", representation(s1="numeric"))
setMethod("+", c("A", "A"), function (e1, e2) print("howdy"))
class(x) = "A"
x + x
asS4(x) + x


###################################################
### chunk: exNotRun eval=FALSE
###################################################
## setMethod("foo", "myclass", myS3Method)
## setMethod("foo", "myclass", function(x, y, ...) myS3Method(x, y, ...))


###################################################
### chunk: testG
###################################################
 testG  = function(x, ...) UseMethod("testG")
 setGeneric("testG")
 getMethod("testG", signature="ANY")


###################################################
### chunk: loadpacks
###################################################
library("graph")
library("Rgraphviz")
library("RBGL")


###################################################
### chunk: showGetClass
###################################################
graphClasses = getClasses("package:graph")
head(graphClasses)


###################################################
### chunk: classGraph
###################################################
graphClassgraph = classList2Graph(graphClasses)


###################################################
### chunk: classGraphComps
###################################################
ccomp = connectedComp(graphClassgraph)
complens = sapply(ccomp, length)
length(ccomp)
table(complens)


###################################################
### chunk: printSingletons
###################################################
unlist(ccomp[complens==1], use.names=FALSE)


###################################################
### chunk: renderG
###################################################
subGnodes = ccomp[[which.max(complens)]]
subG = subGraph(subGnodes, graphClassgraph)
nodeRenderInfo(subG) <- list(shape="ellipse")
attrs = list(node=list(fixedsize = FALSE))
x = layoutGraph(subG, attrs = attrs)
renderGraph(x)


###################################################
### chunk: secondlargestCC
###################################################
whichcc = as.integer(names(sort(complens,decr=T))[2])
subGN = ccomp[[whichcc]]
subG = subGraph(subGN, graphClassgraph)
x = layoutGraph(subG, attrs=list(node=list(shape="ellipse",
    fixedsize=FALSE)))
renderGraph(x)


