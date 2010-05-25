###################################################
### chunk number 1: loadlibs
###################################################
 library("RbcBook1")
 library("ALL")
  
 library("arrayMagic")  ## for the image plots
 stopifnot(package.version("arrayMagic") >= package_version("1.5.4"))


###################################################
### chunk number 2: loadALL
###################################################
library("genefilter")
data(ALL)
Bsub <- (ALL$mol == "BCR/ABL")
Bs <- ALL[,Bsub]


###################################################
### chunk number 3: genefilterALL
###################################################
f1 <- pOverA(.25, log2(100))
f2 <- function(x)(IQR(x)>0.5)
f3 <- function(x) (median(2^x) > 300)
ff <- filterfun(f1,f2,f3)
selected <- genefilter(Bs, ff)
sum(selected)
BSub <- Bs[selected,]
eS <- exprs(BSub)
mads <- apply(eS, 1, mad)
meds <- apply(eS, 1, median)

e1 <- sweep(eS, 1, meds)
e2 <- sweep(e1, 1, mads, FUN="/")

BSubStd <- BSub
exprs(BSubStd) <- e2


###################################################
### chunk number 4: selGenes
###################################################
library("GO")
library("annotate")
GOTERM$"GO:0006917"
library("hgu95av2")

apop <- hgu95av2GO2ALLPROBES$"GO:0006917"
inboth <- apop %in% row.names(e2)
whsel <- apop[inboth]

exprApop <- e2[whsel,]

unlist(mget(whsel, hgu95av2LOCUSID))


###################################################
### chunk number 5: compDists
###################################################
library("bioDist")
man   <- dist(exprApop, "manhattan")
MI       <- MIdist(exprApop)
KLsmooth <- KLD.matrix(exprApop)
KLbin    <- KLdist.matrix(exprApop)


###################################################
### chunk number 6: imageMI
###################################################
library("arrayMagic")
par(mai=c(1.2,1.2,.7,0.1))## to avoid gene names running outside the image
plot.imageMatrix(as.matrix(MI), xlab="", ylab="", main="MI")


###################################################
### chunk number 7: imageKLbin
###################################################
par(mai=c(1.2,1.2,.7,0.1))
plot.imageMatrix(as.matrix(KLbin), xlab="", ylab="", main="KLbin")


###################################################
### chunk number 8: imageKLsmooth
###################################################
par(mai=c(1.2,1.2,.7,0.1)) 
plot.imageMatrix(as.matrix(KLsmooth), xlab="", ylab="", main="KLsmooth")


###################################################
### chunk number 9: compdists
###################################################
v1df <- data.frame(as.vector(MI), as.vector(KLsmooth), as.vector(KLbin))
names(v1df) <- c("MI", "KLsmooth", "KLbin")
pairs(v1df)


###################################################
### chunk number 10: ALLabl
###################################################
abl.locusid <- "25"
hgu95av2.locusid <- unlist(as.list(hgu95av2LOCUSID))
hgu95av2LOCUSID2probes <- split(names(hgu95av2.locusid), hgu95av2.locusid)
abl.probesets <- intersect(hgu95av2LOCUSID2probes[[abl.locusid]], rownames(e2))
pairs(t(e2[abl.probesets,]))


###################################################
### chunk number 11: ALLdist1
###################################################
library("bioDist")
##replace this next call with those below to do all computations
##they are quite extensive though
data(ALL.dist)
##bcr.cor <- cor.dist(e2)
##bcr.spear <- spearman.dist(e2)
##bcr.tau <- tau.dist(e2)
##bcr.euc <- euc(e2)
##bcr.man <- man(e2)
##bcr.kldist <- KLdist.matrix(e2)
##bcr.kld.locfit <- KLD.matrix(e2, method="locfit")
##bcr.kld.density <- KLD.matrix(e2, method="density")
#bcr.kld.normal <- kl.normal(e2, e2.var)
#bcr.mi <- MIdist(e2)

##save(bcr.cor, bcr.spear, bcr.tau, bcr.euc, bcr.man, bcr.kldist, 
##  bcr.mi, file="ALL.dist.rda")


bcr.probesets.top <- list()
top <- 100

for (probeset in abl.probesets) {
  bcr.probesets.top[[probeset]] <- matrix(0,ncol=7,nrow=top)
  colnames(bcr.probesets.top[[probeset]]) <-
      c("cor","spear","tau","euc","man","kld","mi")

bcr.probesets.top[[probeset]][,"cor"] <- closest.top(probeset, bcr.cor, top)
bcr.probesets.top[[probeset]][,"spear"] <- closest.top(probeset, bcr.spear,top)
bcr.probesets.top[[probeset]][,"tau"] <- closest.top(probeset, bcr.tau,top)
bcr.probesets.top[[probeset]][,"euc"] <- closest.top(probeset,bcr.euc,top)
bcr.probesets.top[[probeset]][,"man"] <- closest.top(probeset,bcr.man,top)
bcr.probesets.top[[probeset]][,"kld"] <- closest.top(probeset,bcr.kldist,top)
bcr.probesets.top[[probeset]][,"mi"] <- closest.top(probeset,bcr.mi,top)
}

## compare agreement
## method 1: percentage of agreement
temp <- lapply(lapply(bcr.probesets.top, function(x)
       apply(x,2,function(y) apply(x, 2, function(z)
                     length(intersect(y,z))/top))),function(w) as.dist(w))

sum <- 0
for ( i in 1:length(temp))
{
  sum <- sum + temp[[i]]
}

print(temp)
##cat("Averaging over multiple probesets","\n")
##print(round(sum/length(temp),2))


###################################################
### chunk number 12: ALLdist2
###################################################

  ## method 2: rank of the multiple probesets for the same gene as the
  ## target gene for each distance
for ( i in names(bcr.probesets.top))
{
  cat(i,"\n")
  print(apply(bcr.probesets.top[[i]],2, function(x)
    match(names(bcr.probesets.top)[names(bcr.probesets.top) != i], x)))
}



###################################################
### chunk number 13: ALLadj1
###################################################
library("humanLLMappings")

v1 <- as.list(humanLLMappingsLL2PMID)
v1 <- v1[sapply(v1, function(x) !all(is.na(x)))]

humanLLMappingsLL2PMID <- l2e(v1)

probesets.pubmed <-
    rownames(e2)[unlist(mget(rownames(e2),hgu95av2LOCUSID)) %in%
                 ls(humanLLMappingsLL2PMID)]

## choose top probesets using distances
bcr.probesets.pubmed.top <- list()
top <- 100

for (probeset in abl.probesets)
{

 bcr.probesets.pubmed.top[[probeset]] <- matrix(0,ncol=7,nrow=top)
 colnames(bcr.probesets.pubmed.top[[probeset]]) <-
 c("cor","spear","tau","euc","man","kld","mi")

bcr.probesets.pubmed.top[[probeset]][,"cor"] <- closest.top(probeset, as.dist(as.matrix(bcr.cor)[probesets.pubmed,probesets.pubmed]), top)
bcr.probesets.pubmed.top[[probeset]][,"spear"] <- closest.top(probeset, as.dist(as.matrix(bcr.spear)[probesets.pubmed,probesets.pubmed]),top)
bcr.probesets.pubmed.top[[probeset]][,"tau"] <- closest.top(probeset, as.dist(as.matrix(bcr.tau)[probesets.pubmed,probesets.pubmed]),top)
bcr.probesets.pubmed.top[[probeset]][,"euc"] <- closest.top(probeset,as.dist(as.matrix(bcr.euc)[probesets.pubmed,probesets.pubmed]),top)
bcr.probesets.pubmed.top[[probeset]][,"man"] <- closest.top(probeset,as.dist(as.matrix(bcr.man)[probesets.pubmed,probesets.pubmed]),top)
bcr.probesets.pubmed.top[[probeset]][,"kld"] <- closest.top(probeset,as.dist(as.matrix(bcr.kldist)[probesets.pubmed,probesets.pubmed]),top)
bcr.probesets.pubmed.top[[probeset]][,"mi"] <- closest.top(probeset,as.dist(as.matrix(bcr.mi)[probesets.pubmed,probesets.pubmed]),top)
}



## get adjacency list for ABL1
abl.adj <- unique(unlist(mget(v1[[as.character(abl.locusid)]],
                       humanLLMappingsPMID2LL)))
abl.adj <- abl.adj[abl.adj != abl.locusid]
abl.adj.e2 <- intersect(unique(unlist(mget(rownames(e2), hgu95av2LOCUSID))),
       abl.adj)

temp <- lapply(bcr.probesets.pubmed.top,function(x) {
   apply(x, 2, function(y) 
        intersect(as.character(unlist(mget(y,hgu95av2LOCUSID))), abl.adj.e2))
  })

temp1 <- lapply(temp, function(x) {
     if(is.list(x))
         sapply(x, function(y) length(y))
     else
         apply(x, 2, function(y) length(y))
   })


###################################################
### chunk number 14: ALLadj2
###################################################
 print(sapply(temp1, function(x) x))


###################################################
### chunk number 15: pvcomp
###################################################
nms<-paste("P(X >= ", 1:3, ")", sep="")
pvs <- phyper(0:2, m=3, n=624, k=100, lower.tail = FALSE)
names(pvs) <- nms
print(round(pvs,4))


