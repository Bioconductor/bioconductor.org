###################################################
### chunk number 1: startUP
###################################################
## This rather nasty chunk of code is only necessary for reproducing
## the random division into training and validation sets
library("RbcBook1")
library("kidpack")
data(eset)
set.seed(32)
x        <- t(exprs(eset)) 
y        <- pData(phenoData(eset))$type
y        <- as.factor(y)
y.train  <- y[-c(1:10) ]
x.train  <- x[-c(1:10),]  
t.size   <- length(y.train)
l.size   <- round((2/3)*t.size)
v.size   <- t.size-l.size
l.samp   <- NULL
K        <- nlevels(y.train)
props    <- round(l.size/t.size * table(y.train))
props[1] <- l.size - sum(props[2:K])
for (k in 1:K)
  {
    y.num  <- as.numeric(y.train)
    l.samp <- c(l.samp, sample(which(y.num==k))[1:props[k]])
  }
v.samp     <- (1:t.size)[-l.samp]
library("multtest")
yl.num     <- as.numeric(y.train[l.samp])-1
xl.mtt     <- t(x.train[l.samp,])
f.stat     <- mt.teststat(xl.mtt, yl.num, test = "f")
best.genes <- rev(order(f.stat))[1:200]
x.learn    <- x.train[l.samp, best.genes]
x.valid    <- x.train[v.samp, best.genes]
y.learn    <- y.train[l.samp]
y.valid    <- y.train[v.samp]
library("sma")
yl.num      <- as.numeric(y.learn)
yv.num      <- as.numeric(y.valid)-1
dlda.predic <- stat.diag.da(x.learn, yl.num, x.valid)
dlda.error  <- mean(dlda.predic$pred-1 != yv.num)
library("class")
knn.error <- numeric(3)
for(k in c(1,3,5))
  {
    i              <- ((k-1)/2)+1
    knn.predic     <- knn(x.learn, x.valid, y.learn, k = k)
    knn.error[i]   <- mean(knn.predic!=y.valid)
  }
library("e1071")
svm.error <- matrix(0, nrow = 3, ncol = 3)
for(cost in 0:2)
  {
    for(gamma in (-1):1)
      {
        i              <- cost+1
        j              <- gamma+2
        svm.fit        <- svm(x.learn, y.learn, cost = 2^cost,    
                              gamma = 2^gamma/ncol(x.learn))
        svm.predic     <- predict(svm.fit, newdata = x.valid)
        svm.error[i,j] <- mean(svm.predic != y.valid)
      }
  }
randiv.repr <- function(x,y)
  {
    t.size   <- length(y)
    l.size   <- round((2/3)*t.size)
    v.size   <- t.size-l.size
    l.samp   <- NULL
    K        <- nlevels(y)
    props    <- round(l.size/t.size * table(y))
    props[1] <- l.size - sum(props[2:K])
    for (k in 1:K)
      {
        y.num   <- as.numeric(y)
        l.samp  <- c(l.samp, sample(which(y.num==k))[1:props[k]])
      }
    v.samp      <- (1:t.size)[-l.samp]
    yl.num      <- as.numeric(y[l.samp])-1
    f.stat      <- mt.teststat(t(x[l.samp,]), yl.num, test = "f")
    best.genes  <- rev(order(f.stat))[1:200]
    x.learn     <- x[l.samp, best.genes]
    x.valid     <- x[v.samp, best.genes]
    y.learn     <- y[l.samp]
    y.valid     <- y[v.samp]
    yl.num      <- as.numeric(y.learn)
    yv.num      <- as.numeric(y.valid)-1
    dlda.predic <- stat.diag.da(x.learn, as.numeric(y.learn), x.valid)
    dlda.error  <- mean(dlda.predic$pred-1 != yv.num)
    knn.error   <- numeric(3)
    for(k in c(1,3,5))
      {
        i              <- ((k-1)/2)+1
        knn.predic     <- knn(x.learn, x.valid, y.learn, k = k)
        knn.error[i]   <- mean(knn.predic!=y.valid)
      }
    svm.error <- matrix(0, nrow = 3, ncol = 3)
    for(cost in 0:2)
      {
        for(gamma in (-1):1)
          {
            i              <- cost+1
            j              <- gamma+2
            svm.fit        <- svm(x.learn, y.learn, cost = 2^cost,    
                                  gamma = 2^gamma/ncol(x.learn))
            svm.predic     <- predict(svm.fit, newdata = x.valid)
            svm.error[i,j] <- mean(svm.predic != y.valid)
          }
      } 
    list(dlda = dlda.error, knn = knn.error, svm = svm.error,
         lsamp = l.samp, vsamp = v.samp)
  }
runs         <- 50
dlda.errors  <- numeric(runs)
knn.errors   <- matrix(0, nrow = 3, ncol = runs)
svm.errors   <- array(0, dim = c(3, 3, runs))
vsamp        <- matrix(0, runs, length(v.samp))
lsamp        <- matrix(0, runs, length(l.samp))
for(r in 1:runs)
  {
    results         <- randiv.repr(x.train, y.train)
    dlda.errors[ r] <- results$dlda
    knn.errors[, r] <- results$knn
    svm.errors[,,r] <- results$svm
    vsamp[r,]       <- results$vsamp
    lsamp[r,]       <- results$lsamp
    cat("This was run", r, "of", runs, "\n")
  }


###################################################
### chunk number 2: loadkidney
###################################################
set.seed(32)
library("kidpack")
data(eset)


###################################################
### chunk number 3: wegzehn
###################################################
test   <- (1:10)
train  <- (1:length(eset$type))[-test]  # needed!
trEset <- eset[,train]


###################################################
### chunk number 4: 
###################################################
t.size <- length(trEset$type)
l.size <- round((2/3)*t.size)
v.size <- t.size-l.size


###################################################
### chunk number 5: 
###################################################
K        <- nlevels(factor(trEset$type))
l.samp   <- NULL
props    <- round(l.size/t.size * table(trEset$type))
props[1] <- l.size - sum(props[2:K])
for (k in 1:K)
  {
    y.num  <- as.numeric(factor(trEset$type))
    l.samp <- c(l.samp, sample(which(y.num==k))[1:props[k]])
  }
v.samp      <- (1:t.size)[-l.samp]


###################################################
### chunk number 6: 
###################################################
table(trEset$type[l.samp])
table(trEset$type[v.samp])


###################################################
### chunk number 7: 
###################################################
library("multtest")
yl.num     <- as.numeric(factor(trEset$type[l.samp]))-1
xl.mtt     <- exprs(trEset[,l.samp])
f.stat     <- mt.teststat(xl.mtt, yl.num, test = "f")
best.genes <- rev(order(f.stat))[1:200]
trselEset  <- trEset[best.genes,]


###################################################
### chunk number 8: 
###################################################
library("sma")
library("MLInterfaces")
l.samp      <- as.integer(l.samp)
dlda.predic <- stat.diag.daB(trselEset, "type", l.samp)
conf.matrix <- confuMat(dlda.predic)
error.rate  <- function(cm) 1-sum(diag(cm))/sum(cm)
dlda.error  <- error.rate(conf.matrix)


###################################################
### chunk number 9: 
###################################################
library("class")
knn.error <- numeric(3)
for(k in c(1,3,5))
 {
   i            <- ((k-1)/2)+1
   knn.predic   <- knnB(trselEset, "type", l.samp, k = k, prob = FALSE)
   knn.error[i] <- error.rate(confuMat(knn.predic))
 }


###################################################
### chunk number 10: 
###################################################
library("e1071")
svm.error <- matrix(0, nrow = 3, ncol = 3)
for(cost in 0:2)
  {
    for(gamma in (-1):1)
      {
        i       <- cost+1
        j       <- gamma+2
        svm.fit <- svmB(trselEset, "type", l.samp, cost=2^cost, 
                        gamma = 2^gamma/nrow(exprs(trselEset)),
                        type="C-classification")
        svm.error[i,j] <- error.rate(confuMat(svm.fit))
      }
  }


###################################################
### chunk number 11: 
###################################################
dlda.error
knn.error
svm.error


###################################################
### chunk number 12: 
###################################################
comment <- function(x) x


###################################################
### chunk number 13: 
###################################################
randiv <- function(trEset)
  {
    "\#body suppressed; repeat all the code from Section 24.3"
    list(dlda=dlda.error, knn=knn.error, svm=svm.error)
  }


###################################################
### chunk number 14: 
###################################################
randiv <- function(trEset, run)
  {
    ## Defining the size of learning and validation sets
    t.size <- length(trEset$type)
    l.size <- round((2/3)*t.size)
    v.size <- t.size-l.size
    
    ## Balanced sampling
    v.samp <- vsamp[run,]
    l.samp <- lsamp[run,]

    ## Variable selection
    yl.num     <- as.numeric(factor(trEset$type[l.samp]))-1
    xl.mtt     <- exprs(trEset[,l.samp])
    f.stat     <- mt.teststat(xl.mtt, yl.num, test = "f")
    best.genes <- rev(order(f.stat))[1:200]
    trselEset  <- trEset[best.genes,]
    
    ## Classification with DLDA
    l.samp      <- as.integer(l.samp)
    dlda.predic <- stat.diag.daB(trselEset, "type", l.samp)
    conf.matrix <- confuMat(dlda.predic)
    dlda.error  <- error.rate(conf.matrix)
                                        
    ## Classification with Nearest Neighbors 
    knn.error <- numeric(3)
    for(k in c(1,3,5))
      {
        i            <- ((k-1)/2)+1
        knn.predic   <- knnB(trselEset, "type", l.samp, k = k, prob = FALSE)
        knn.error[i] <- error.rate(confuMat(knn.predic))
      }
    
    ## Classification with Support Vector Machines
    svm.error <- matrix(0, nrow = 3, ncol = 3)
    for(cost in 0:2)
      {
        for(gamma in (-1):1)
          {
            i       <- cost+1
            j       <- gamma+2
            svm.fit <- svmB(trselEset, "type", l.samp, 
                            cost=2^cost, 
                            gamma = 2^gamma/nrow(exprs(trselEset)),
                            type="C-classification")
            svm.error[i,j] <- error.rate(confuMat(svm.fit))
          }
      }

    ## Output
    list(dlda = dlda.error, knn = knn.error, svm = svm.error)
  }


###################################################
### chunk number 15: 
###################################################
runs        <- 50
dlda.errors <- numeric(runs)
knn.errors  <- matrix(0, nrow = 3, ncol = runs)
svm.errors  <- array(0, dim = c(3, 3, runs))
for(r in 1:runs)
 {
   results         <- randiv(trEset, r)
   dlda.errors[ r] <- results$dlda
   knn.errors[, r] <- results$knn
   svm.errors[,,r] <- results$svm
   cat("This was run", r, "of", runs, "\n")
 }


###################################################
### chunk number 16: 
###################################################
mean(dlda.errors)
apply(knn.errors, 1, mean)
apply(svm.errors, c(1,2), mean)


###################################################
### chunk number 17: 
###################################################
## These functions are for drawing boxplots and barplots of error rates

## Loading the required package
library("grid")

## Check the settings
if (!exists("pushViewport")) pushViewport <- push.viewport
if (!exists("popViewport"))  popViewport  <- pop.viewport

## The function for the barplots
grid.barplot <- function(x, xlab = FALSE, main = FALSE, xscale = NULL, yticks,
                         ylabels)
{
  ## Customizing the data
  y      <- table(x)
  y.alt  <- table(x)
  x.alt  <- x
  x      <- sort(unique(x))
  if(!is.null(xscale))
    {
      index <- which(x > xscale[1] & x < xscale[2])
      x <- x[index]
      y <- y[index]
    }

  ## Horizontal bars
  grid.segments(x0 = unit(x, "native"), y0 = unit(0, "native"),
                x1 = unit(x, "native"), y1 = unit(y.alt, "native"))

  ## Drawing the lines
  grid.lines(unit(x, "native"), unit(y, "native"))

  ## Plotting the mean
  brks   <- as.numeric(attributes(table(x))$dimnames$x)
  start  <- max(which(brks<mean(x.alt)))
  ende   <- start+1
  xinc   <- (mean(x.alt)-brks[start])/(brks[ende]-brks[start])
  yh     <- ((y.alt[ende]-y.alt[start])*xinc)+(y.alt[start])             
  farbe  <- gpar(col="red")
  grid.segments(x0=unit(mean(x.alt),"native"), y0=unit(0, "native"),
                x1=unit(mean(x.alt),"native"), y1=unit(yh,"native"), gp=farbe)

  ## Annotation
  if(xlab) grid.xaxis()
  grid.yaxis(main = FALSE, at=yticks, name=ylabels)
  grid.rect()
}


## The function for the boxplots
grid.boxplot <- function(x, xlab=FALSE, main=FALSE, ylab=NULL, xscale=NULL)
{
  ## Generic boxplot
  x    <- boxplot(x, plot = FALSE)
  farb <- gpar(col="red")

  ## Transformation to the grid-system
  grid.lines(unit(x$stats[1],  "native"), unit(c(0.3, 0.7), "npc"))
  grid.lines(unit(x$stats[1:2],"native"), unit(0.5, "npc"), gp = gpar(lty = 2))
  grid.rect(unit(x$stats[2], "native"), y = unit(0.5, "npc"),
            height = unit(0.6, "npc"), width=unit(diff(x$stats[2:3]),"native"),
            just = c("left", "centre"))
  grid.rect(unit(x$stats[3], "native"), y = unit(0.5, "npc"),
            height = unit(0.6, "npc"), width=unit(diff(x$stats[3:4]),"native"),
            just = c("left", "centre"))
  grid.lines(unit(x$stats[4:5],"native"), unit(0.5, "npc"), gp = gpar(lty = 2))
  grid.lines(unit(x$stats[5], "native"), unit(c(0.3, 0.7), "npc"))
  grid.lines(unit(x$stats[3], "native"), unit(c(0.2, 0.8), "npc"), gp = farb)

  ## Outliers
  n <- length(x$out)
  if(n > 0) {
    if(!is.null(xscale)) index <- which(x$out > xscale[1] & x$out < xscale[2])
      else index <- 1:n
    if(length(index) > 0)
      {
        grid.points(unit(x$out[index],"native"),unit(rep(0.5,length(index)),
                    "npc"), size = unit(0.5, "char"))
      }
  }

  ## Annotation
  if(xlab) grid.xaxis()
  if(!is.null(ylab))
    {
      grid.text(ylab, x = unit(0, "npc") - unit(0.5, "lines"),
                gp = gpar(fontface = "bold"), just = c("right", "centre"))
    }
}


## Plotting the results
grid.boxNbar <- function(daten, mname)
  {
    grid.newpage()
    pushViewport(plotViewport(c(4, 6, 4, 4)))
    #pushViewport(viewport(layout=lyt))
    n         <- ncol(daten)

    ## Calibration of the x-axis
    xsc     <- range(daten) + c(-1, 1) * max(daten)/15
    xsc2    <- xsc          

    ## Calibration of the y-axis
    ysc     <- c(0,max(sapply(apply(daten,2,function(x) table(x)),max)))
    ysc     <- ysc * 1.1
    ylabels <- c("0", "10", "20", "30", "40", "50")
    ynumb   <- floor((ysc[2]+3)/10)-1
    yticks  <- (0:ynumb)*10
    ylabels <- ylabels[1:(ynumb+1)]
          
    ## The plot contains 3 pieces
    lyt <- grid.layout(1,3,widths=unit(c(0.5-0.03/2,0.03,0.5-0.03/2),"npc"))
    pushViewport(viewport(layout=lyt))

    ## Left side: boxplots
    lyt <- grid.layout(n,1,heights=unit(rep(1/n, n), "npc"))
    pushViewport(viewport(layout.pos.row=1, layout.pos.col=1,xsc=xsc))
    pushViewport(viewport(layout=lyt))
    grid.rect()
    for(i in n:1)
      {
        pushViewport(viewport(layout.pos.r=i,layout.pos.c=1,xsc=xsc))
        main  <- (i < 2)
        grid.boxplot(daten[,i],m=main,xl=(i>(n-1)),yl=mname[i],xsc=xsc)
        popViewport(1)
      }
    popViewport(1)
    grid.text("Boxplot", y = unit(-2.5, "lines"))

    ## Middle: an empty frame  
    popViewport(1)
    #schrift <- gpar(fontface="bold")
    #grid.text(namen[j], y=unit(1,"npc")+unit(2,"lines"), gp=schrift)

    ## Right side: density plots
    lyt <- grid.layout(n,1,hei=unit(rep(1/n,n),"npc"))
    pushViewport(viewport(layout.pos.row=1,layout.pos.col=3,xscale=xsc))
    pushViewport(viewport(layout=lyt))
    grid.rect()
    for(i in n:1)
      {
        pushViewport(viewport(layout.pos.r=i,layout.pos.c=1,xs=xsc2,ys=ysc))
        main  <- (i < 2)
        grid.barplot(daten[,i],main,xsc2,yticks,ylabels,xlab=(i>(n-1)))
        popViewport(1)
      }
    popViewport(1)
    grid.text("Density", y=unit(-2.5,"lines"))
    popViewport(3)
  }


###################################################
### chunk number 18: 
###################################################
daten <- cbind(dlda.errors, knn.errors[2,], svm.errors[3,1,], svm.errors[2,1,])
mname <- c("DLDA", "kNN", "SVM1", "SVM2")
grid.boxNbar(daten, mname)


###################################################
### chunk number 19: 
###################################################
yt.num      <- as.numeric(factor(trEset$type))-1
xt.mtt      <- exprs(trEset)
f.stat      <- mt.teststat(xt.mtt, yt.num, test = "f")
best.genes  <- rev(order(f.stat))[1:200]
selEset     <- eset[best.genes,]
dlda.predic <- stat.diag.daB(selEset, "type", train)
dlda.predic@RObject$pred-1


###################################################
### chunk number 20: 
###################################################
confuMat(dlda.predic)


