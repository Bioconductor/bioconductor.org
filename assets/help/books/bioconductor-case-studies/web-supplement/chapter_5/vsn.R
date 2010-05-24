###################################################
### chunk: exrn1
###################################################
r1 = rnorm(10000, mean=2000, sd=50)
g1 = rnorm(10000, mean=1000, sd=50)
hist(r1/g1, breaks=33, col="azure")


###################################################
### chunk: exrn2
###################################################
r2 = rnorm(10000, mean=200, sd=50)
g2 = rnorm(10000, mean=100, sd=50)
ratio = r2/g2 
ratio = ratio[(ratio>0)&(ratio<6)]
hist(ratio, breaks=33, col="azure", 
     main="Histogram of r2/g2", xlab="r2/g2",
     sub="restricted to [0,6]")


###################################################
### chunk: exrn3
###################################################
hist(log2(r1/g1), breaks=33, col="azure")


###################################################
### chunk: exrn4
###################################################
hist(log2(r2/g2), breaks=33, col="azure")


###################################################
### chunk: loadkidney
###################################################
library("vsn")
data("kidney")


###################################################
### chunk: loadCCl4
###################################################
library("CCl4")
data("CCl4")


###################################################
### chunk: helpCCl4 eval=FALSE
###################################################
## ? CCl4


###################################################
### chunk: CCl4s
###################################################
selArrays = with(pData(CCl4), 
    (Cy3 == "CCl4" & RIN.Cy3 > 9) | 
    (Cy5 == "CCl4" & RIN.Cy5 > 9))
selFeatures = !is.na(featureData(CCl4)$ID) 
CCl4s = CCl4[selFeatures, selArrays]


###################################################
### chunk: denskid
###################################################
library("geneplotter")
pcol = c("green3", "red1")
plty = 1:2
multidensity(exprs(kidney), xlim=c(-200, 1000),
    main = "kidney", xlab="Intensity",
    lty = plty, col = pcol, lwd = 2)
legend("topright", c("green", "red"), 
    lty = plty, col = pcol, lwd = 2)


###################################################
### chunk: densccl4
###################################################
multidensity(cbind(assayData(CCl4s)$G[,1], 
    assayData(CCl4s)$R[,1]), xlim=c(0, 200),
    main = expression(CCl[4]), lwd=2, xlab="Intensity",
    col  = pcol, lty  = plty)


###################################################
### chunk: figgraphs
###################################################
px = seq(-100, 500, length=50)
f  = function(x, b) log2(x+b)
h  = function(x, a) log2((x+sqrt(x^2+a^2))/2)
matplot(px, y=cbind(h(px, a=50), f(px, b=50)), 
  type="l", lty=1:2, xlab="x", ylab="f, h")


###################################################
### chunk: exgraphs
###################################################
px = seq(0, 1e8, length=50)


###################################################
### chunk: figaffine
###################################################
axl = c(30, 300)
plot(assayData(CCl4s)$R[,1], 
     assayData(CCl4s)$G[,1],
     xlim=axl, ylim=axl, pch=".", col="grey", 
     asp=1)
abline(a=0, b=1, col="blue", lty=2, lwd=3)
abline(a=18, b=1.2, col="red", lty=3, lwd=3)


###################################################
### chunk: vsn2
###################################################
CCl4sn = justvsn(CCl4s, backgroundsubtract=TRUE)
class(CCl4sn)


###################################################
### chunk: asd
###################################################
asd = assayData(CCl4sn)
A = (asd$R+asd$G)/2
M = (asd$R-asd$G) 


###################################################
### chunk: MAplot1
###################################################
plot(A[,6], M[,6], pch='.', asp=1, xlab="A", ylab="M")
abline(h=0, col="blue")


###################################################
### chunk: MAplot2
###################################################
smoothScatter(A[,6], M[,6], nrpoints=300, asp=1, 
    xlab="A", ylab="M")
abline(h=0, col="blue")


###################################################
### chunk: Mdensity
###################################################
multidensity(M, xlim=c(-2,2), bw=0.1)
abline(v=0, col="grey")


###################################################
### chunk: meanSdCCl4
###################################################
meanSdPlot(cbind(assayData(CCl4sn)$R, assayData(CCl4sn)$G),
    ylim=c(0, 1.4))


###################################################
### chunk: kidspike
###################################################
kidspike = kidney

## The selection of spike-in features 'sel' is slightly different below
## from how it was done in Edition 1 of the book, namely, the dimmest
## 20% of the features are not selected. This is more realistic.
m = rowMeans(exprs(kidspike))
sel = (runif(length(m)) < 1/3) & (m > quantile(m, 0.2))

delta = 100 + 0.4*abs(m[sel])
exprs(kidspike)[sel,] = exprs(kidspike)[sel,] + 
    cbind(-delta,+delta)


###################################################
### chunk: vkid
###################################################
ltsq = c(1, 0.8, 0.5)
vkid = lapply(ltsq, function(p) vsn2(kidspike, 
    lts.quantile=p))


###################################################
### chunk: getMA
###################################################
getMA = function(x) 
    data.frame(A = rowSums(exprs(x))/2,
        M = as.vector(diff(t(exprs(x)))))
ma = lapply(vkid, getMA)


###################################################
### chunk: rbind
###################################################
for(i in seq(along=ma)) {
    ma[[i]]$group = factor(ifelse(sel, "up", "unchanged"))
    ma[[i]]$lts.quantile = factor(ltsq[i])
}
ma = do.call("rbind", ma)


###################################################
### chunk: vkidplot
###################################################
library("lattice")
lp = xyplot(M ~ A | lts.quantile, group=group, data=ma, 
    layout = c(1,3), pch=".", ylim=c(-3,3), xlim=c(6,16),
    auto.key=TRUE,
    panel = function(...){ 
        panel.xyplot(...) 
        panel.abline(h=0)},
    strip = function(...) strip.default(..., 
        strip.names=TRUE, strip.levels=TRUE))
print(lp)


###################################################
### chunk: normctrl
###################################################
normctrl = sample(which(!sel), 100)
fit = vsn2(kidspike[normctrl, ], lts.quantile=1)
vkidctrl = predict(fit, newdata=kidspike)

ma = getMA(vkidctrl)
ma$group = factor("other", levels=c("other", 
    "normalization control"))
ma$group[normctrl] = "normalization control"


###################################################
### chunk: normctrlplot
###################################################
lp = xyplot(M ~ A , group=group, data=ma,  
    pch=".", ylim=c(-3,3), xlim=c(4,16), auto.key=TRUE,
    panel = function(...){ 
        panel.xyplot(...)
        panel.abline(h=0)})
print(lp)


###################################################
### chunk: esG
###################################################
esG  = channel(CCl4s, "G")
exprs(esG) = exprs(esG) - exprs(channel(CCl4s, "Gb"))


###################################################
### chunk: esGvsn
###################################################
nesG = justvsn(esG)


###################################################
### chunk: esGecdf1
###################################################
library("RColorBrewer")
startedlog = function(x) log2(x+5)
colors = brewer.pal(6, "Dark2")
multiecdf(startedlog(exprs(esG)), col=colors, 
    main="Before normalization")


###################################################
### chunk: esGecdf2
###################################################
multiecdf(exprs(nesG), col=colors, 
    main="After normalization")


###################################################
### chunk: esGscp
###################################################
smoothScatter(getMA(nesG[, 1:2]))
abline(h=0, col="red")


###################################################
### chunk: shrink1
###################################################
g = exprs(kidney)[,1]
r = exprs(kidney)[,2]


###################################################
### chunk: shrink2
###################################################
Aspike = 2^seq(2, 16, by=0.5)
sel = sample(nrow(kidney), length(Aspike))
r[sel] = Aspike*sqrt(2)
g[sel] = Aspike/sqrt(2)


###################################################
### chunk: shrink3
###################################################
A = (log2(r) + log2(g))/2
M.naive =  log2(r) - log2(g)
fit = vsn2(cbind(g, r))
M.vsn = exprs(fit)[,2] - exprs(fit)[,1]


###################################################
### chunk: figshrink
###################################################
plot(A, M.naive, pch=".", ylab="M", 
  xlim=c(2,16), ylim=c(-3,3), col="grey")
points(A, M.vsn, pch=".", col="mistyrose")

sel = sel[order(A[sel])]
lines(A[sel], M.vsn[sel], col="red", lwd=2)
lines(A[sel], M.naive[sel], lwd=2, lty=2)


###################################################
### chunk: ref1
###################################################
ref = vsn2(CCl4s[, 1:5], backgroundsubtract=TRUE)


###################################################
### chunk: ref2
###################################################
x6 = justvsn(CCl4s[, 6], reference = ref, 
    backgroundsubtract=TRUE)


###################################################
### chunk: figref
###################################################
plot(assayData(x6)$G, assayData(CCl4sn)$G[,6], pch=".", 
    asp=1)
abline(a=0, b=1, col="red")


