###################################################
### chunk number 1: setup
###################################################
library("RbcBook1")


###################################################
### chunk number 2: limmaSetup
###################################################
library("marray")
library("limma")
home.dir <- getwd()
generateMA <- function(narrays=2) {
	MA <- new("MAList")
	gn <- paste("Gene",1:4,sep="")
	an <- paste("Array",1:narrays,sep="")
	G <- 2^matrix(rnorm(4*narrays,mean=8,sd=2),4,narrays,dimnames=list(gn,an))
	R <- 2^matrix(rnorm(4*narrays,mean=8,sd=2),4,narrays,dimnames=list(gn,an))
	MA$M <- log2(R)-log2(G)
	MA$A <- (log2(R)+log2(G))/2
	MA$genes <- data.frame(ID=I(paste("Gene",1:4,sep="")))
	MA
}


###################################################
### chunk number 3: 
###################################################
MA <- generateMA(3)


###################################################
### chunk number 4: ReplicatedArrays
###################################################
fit <- lmFit(MA)
fit <- eBayes(fit)
topTable(fit, adjust="fdr")


###################################################
### chunk number 5: 
###################################################
targets <- data.frame(FileName=I(c("File1","File2","File3")),Cy3=I(c("wt","mu","wt")),Cy5=I(c("mu","wt","mu")))
row.names(targets) <- c("File1","File2","File3")


###################################################
### chunk number 6: DyeSwap
###################################################
design <- c(1,-1,1)
fit <- lmFit(MA, design)
fit <- eBayes(fit)
topTable(fit, adjust="fdr")


###################################################
### chunk number 7: DyeSwapDesign
###################################################
design <- modelMatrix(targets, ref="wt")


###################################################
### chunk number 8: Beta7ReadGenePix
###################################################
beta7.dir <- system.file("beta7",package="beta7")
setwd(beta7.dir)
TargetInfo <- read.marrayInfo("TargetBeta7.txt")
mraw <- cache("mraw",read.GenePix(targets = TargetInfo))
normdata <- cache("normdata",maNorm(mraw))
setwd(home.dir)


###################################################
### chunk number 9: beta7DyeEffectFit
###################################################
design <- cbind(Dye=1,Beta7=c(1,-1,-1,1,1,-1))
fit <- lmFit(normdata, design, weights=NULL)
fit <- eBayes(fit)


###################################################
### chunk number 10: Beta7DyeEffectTopTable
###################################################
topTable(fit,coef="Dye",adjust="fdr")


###################################################
### chunk number 11: Beta7Weights
###################################################
w <- 0+(maW(normdata) >= 0)
fit <- lmFit(normdata, design, weights=w)
fit <- eBayes(fit)
tab <- topTable(fit,coef="Beta7",adjust="fdr")
tab$Name <- substring(tab$Name,1,20)
tab


###################################################
### chunk number 12: 
###################################################
MA <- generateMA(4)


###################################################
### chunk number 13: TechRep
###################################################
biolrep <- c(1,1,2,2)
corfit <- duplicateCorrelation(MA, ndups=1, block=biolrep)
fit <- lmFit(MA, block=biolrep, cor=corfit$consensus)
fit <- eBayes(fit)
topTable(fit, adjust="fdr")


###################################################
### chunk number 14: TechRepDyeSwap
###################################################
design <- c(1,-1,1,-1)
corfit <- duplicateCorrelation(MA,design,ndups=1,block=biolrep)
fit <- lmFit(MA, design, block=biolrep, cor=corfit$consensus)
fit <- eBayes(fit)
topTable(fit, adjust="fdr")


###################################################
### chunk number 15: OddSwapTargets
###################################################
targets <- data.frame(FileName=I(c("File1","File2","File3","File4")),Cy3=I(c("wt1","mu1","wt2","wt2")),Cy5=I(c("mu1","wt1","mu2","mu2")))
row.names(targets) <- c("File1","File2","File3","File4")
MA <- generateMA(4)


###################################################
### chunk number 16: OddSwapFit
###################################################
design <- cbind(MU1vsWT1=c(1,-1,0,0),MU2vsWT2=c(0,0,1,1))
fit <- lmFit(MA, design)


###################################################
### chunk number 17: OddSwapFit
###################################################
cont.matrix <- makeContrasts(MUvsWT=(MU1vsWT1+MU2vsWT2)/2, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, adjust="fdr")


###################################################
### chunk number 18:  eval=FALSE
###################################################
## targets


###################################################
### chunk number 19:  eval=FALSE
###################################################
## design <- modelMatrix(targets, ref="wt1")
## design <- cbind(Dye=1,design)
## colnames(design)


###################################################
### chunk number 20:  eval=FALSE
###################################################
## fit <- lmFit(MA, design)
## cont.matrix <- makeContrasts(muvswt=(mu1+mu2+mu3-wt2-wt3)/3,
##                                levels=design)
## fit2 <- contrasts.fit(fit, cont.matrix)
## fit2 <- eBayes(fit2)
## topTable(fit2, adjust="fdr")


###################################################
### chunk number 21: AUTargets
###################################################
targets <- data.frame(FileName=I(c("File1","File2","File3","File4")),Cy3=I(c("U1","A1","U2","A2")),Cy5=I(c("A1","U2","A2","U1")))
row.names(targets) <- c("File1","File2","File3","File4")
MA <- generateMA(4)


###################################################
### chunk number 22: AUFit
###################################################
design <- modelMatrix(targets, ref="U1")
fit <- lmFit(MA, design)
cont.matrix <- makeContrasts(AvsU=(A1+A2-U2)/2, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, adjust="fdr")


###################################################
### chunk number 23:  eval=FALSE
###################################################
## corfit <- duplicateCorrelation(MA,design,ndups=2,spacing="columns")


###################################################
### chunk number 24:  eval=FALSE
###################################################
## fit <- lmFit(MA, design,
##              ndups=2, spacing=1, cor=corfit$consensus)


###################################################
### chunk number 25:  eval=FALSE
###################################################
## design


###################################################
### chunk number 26:  eval=FALSE
###################################################
## fit <- lmFit(MA, design)
## fit <- eBayes(fit)
## topTable(fit, coef="MUvsWT", adjust="fdr")


###################################################
### chunk number 27:  eval=FALSE
###################################################
## fit <- lmFit(MA, design)
## fit2 <- contrasts.fit(fit, c(-1,1))
## fit2 <- eBayes(fit2)
## topTable(fit2, adjust="fdr")


###################################################
### chunk number 28: TwoGrpGroup
###################################################
Group <- factor(c("WT","WT","Mu","Mu","Mu"), levels=c("WT","Mu"))


###################################################
### chunk number 29: TwoGrpDesign11
###################################################
design <- cbind(WTvsRef=1,MUvsWT=c(0,0,1,1,1))


###################################################
### chunk number 30: TwoGrpDesign12
###################################################
design <- model.matrix(~Group)
colnames(design) <- c("WTvsRef","MUvsWT")


###################################################
### chunk number 31: TwoGrpDesign21
###################################################
design <- cbind(WT=c(1,1,0,0,0),MU=c(0,0,1,1,1))


###################################################
### chunk number 32: TwoGrpDesign22
###################################################
design <- model.matrix(~0+Group)
colnames(design) <- c("WT","Mu")


###################################################
### chunk number 33:  eval=FALSE
###################################################
## design


###################################################
### chunk number 34: SeveralGrpTargets
###################################################
targets <- data.frame(Targets=I(c("RNA1","RNA2","RNA3","RNA1","RNA2","RNA3")))
eset <- matrix(rnorm(40*6),40,6)


###################################################
### chunk number 35: SeveralGrpDesign
###################################################
f <- factor(targets$Target, levels=c("RNA1","RNA2","RNA3"))
design <- model.matrix(~0+f)
colnames(design) <- c("RNA1","RNA2","RNA3")


###################################################
### chunk number 36: SeveralGrpFit
###################################################
fit <- lmFit(eset, design)
contrast.matrix <- makeContrasts(RNA2-RNA1, RNA3-RNA2, RNA3-RNA1,
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


###################################################
### chunk number 37: SeveralGrpTopTable1
###################################################
topTable(fit2, coef=1, adjust="fdr")


###################################################
### chunk number 38: SeveralGrpDecideTests
###################################################
results <- decideTests(fit2)


###################################################
### chunk number 39: SeveralGrpVennDiag eval=FALSE
###################################################
## vennDiagram(results)


###################################################
### chunk number 40: SeveralGrpF
###################################################
o <- order(fit2$F.p.value)
fit2$genes[o[1:30],]


###################################################
### chunk number 41:  eval=FALSE
###################################################
## design <- modelMatrix(targets, ref="Ref")


###################################################
### chunk number 42:  eval=FALSE
###################################################
## design <- modelMatrix(targets, ref="CD4")
## design


###################################################
### chunk number 43:  eval=FALSE
###################################################
## fit <- lmFit(MA, design)


###################################################
### chunk number 44:  eval=FALSE
###################################################
## cont.matrix <- cbind("CD8-CD4"=c(1,0), "DN-CD4"=c(0,1),
##      "CD8-DN"=c(1,-1))
## fit2 <- contrasts.fit(fit, cont.matrix)
## fit2 <- eBayes(fit2)


###################################################
### chunk number 45: FactorialTargets
###################################################
targets <- data.frame(
FileName=I(c("File1","File2","File3","File4","File5")),
Strain=I(c("WT","WT","Mu","Mu","Mu")),
Treatment=I(c("U","S","U","S","S"))
)
rownames(targets) <- c("File1","File2","File3","File4","File5")
eset <- matrix(rnorm(4*5),4,5)


###################################################
### chunk number 46: TS
###################################################
TS <- paste(targets$Strain, targets$Treatment, sep=".")
TS


###################################################
### chunk number 47: FactorialFit
###################################################
TS <- factor(TS, levels=c("WT.U","WT.S","Mu.U","Mu.S"))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
fit <- lmFit(eset, design)


###################################################
### chunk number 48: FactorialContrasts
###################################################
cont.matrix <- makeContrasts(
    WT.SvsU=WT.S-WT.U,
    Mu.SvsU=Mu.S-Mu.U,
    Diff=(Mu.S-Mu.U)-(WT.S-WT.U),
    levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


###################################################
### chunk number 49: FactorialResults
###################################################
results <- decideTests(fit2)
vennDiagram(results)


###################################################
### chunk number 50: FactorialModel.Matrix
###################################################
Strain <- factor(targets$Strain, levels=c("WT","Mu"))
Treatment <- factor(targets$Treatment, levels=c("U","S"))
design <- model.matrix(~Strain*Treatment)


###################################################
### chunk number 51: FactorialFit2
###################################################
fit <- lmFit(eset, design)
cont.matrix <- cbind(WT.SvsU=c(0,0,1,0),
                     Mu.SvsU=c(0,0,1,1),
                     Diff   =c(0,0,0,1))
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


###################################################
### chunk number 52: TimeCourseTargets
###################################################
targets <- data.frame(
FileName=I(paste("File",1:8,sep="")),
Target=I(c("wt.0hr","wt.0hr","wt.6hr","wt.24hr","mu.0hr","mu.0hr","mu.6hr","mu.24hr"))
)
rownames(targets) <- targets$FileName
eset <- matrix(rnorm(4*8),4,8)


###################################################
### chunk number 53: TimeCourseFit
###################################################
lev <- c("wt.0hr","wt.6hr","wt.24hr","mu.0hr","mu.6hr","mu.24hr")
f <- factor(targets$Target, levels=lev)
design <- model.matrix(~0+f)
colnames(design) <- lev
fit <- lmFit(eset, design)


###################################################
### chunk number 54: TimeCourseWT
###################################################
cont.wt <- makeContrasts(
     "wt.6hr-wt.0hr",
     "wt.24hr-wt.6hr",
levels=design)
fit2 <- contrasts.fit(fit, cont.wt)
fit2 <- eBayes(fit2)


###################################################
### chunk number 55: TimeCourseSelectWT
###################################################
sel.wt <- p.adjust(fit2$F.p.value, method="fdr") < 0.05


###################################################
### chunk number 56: TimeCourseSelectMu
###################################################
cont.mu <- makeContrasts(
     "mu.6hr-mu.0hr",
     "mu.24hr-mu.6hr",
levels=design)
fit2 <- contrasts.fit(fit, cont.mu)
fit2 <- eBayes(fit2)
sel.mu <- p.adjust(fit2$F.p.value, method="fdr") < 0.05


###################################################
### chunk number 57: TimeCourseSelectDiff
###################################################
cont.dif <- makeContrasts(
    Dif6hr =(mu.6hr-mu.0hr)-(wt.6hr-wt.0hr),
    Dif24hr=(mu.24hr-mu.6hr)-(wt.24hr-wt.6hr),
levels=design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
sel.dif <- p.adjust(fit2$F.p.value, method="fdr") < 0.05


###################################################
### chunk number 58:  eval=FALSE
###################################################
## tstat.ord <- fit$coef / fit$stdev.unscaled / fit$sigma


###################################################
### chunk number 59: Beta7Read
###################################################
beta7.dir <- system.file("beta7",package="beta7")
targets <- readTargets("TargetBeta7.txt",path=beta7.dir)
f <- function(x) as.numeric(x$Flags > -75)
RG <- read.maimages(targets$FileName,source="genepix",path=beta7.dir,wt.fun=f)
RG$printer <- getLayout(RG$genes)


###################################################
### chunk number 60: Beta7Background
###################################################
RGne <- backgroundCorrect(RG,method="normexp",offset=25)


###################################################
### chunk number 61: Beta7Normalize
###################################################
MA <- normalizeWithinArrays(RGne)
design <- cbind(Dye=1,Beta7=c(1,-1,-1,1,1,-1))


###################################################
### chunk number 62: Beta7Fit
###################################################
isGene <- maControls(mraw)=="probes"
fit <- lmFit(MA[isGene,], design)
fit <- eBayes(fit)
tab <- topTable(fit,coef="Beta7",adjust="fdr")
tab$Name <- substring(tab$Name,1,20)
tab[,-(1:3)]


###################################################
### chunk number 63: Beta7vsn eval=FALSE
###################################################
## MA <- normalizeBetweenArrays(RG,method="vsn")


###################################################
### chunk number 64: RestoreOptions
###################################################
setwd(home.dir)


