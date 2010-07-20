###################################################
### chunk number 1: loadPacks
###################################################
library("RbcBook1")
library("Biobase")
library("multtest")
library("genefilter")


###################################################
### chunk number 2: classMTP
###################################################
slotNames("MTP")


###################################################
### chunk number 3: ALL
###################################################
library("ALL")
library("hgu95av2")
data(ALL)


###################################################
### chunk number 4: genefilter
###################################################
ffun <- filterfun(pOverA(p=0.2, A=100), cv(a=0.7, b=10))
filt <- genefilter(2^exprs(ALL), ffun)
filtALL <- ALL[filt,]


###################################################
### chunk number 5: ALLDetails
###################################################
filtX <- exprs(filtALL)
pheno <- pData(filtALL)


###################################################
### chunk number 6: BT
###################################################
table(pData(ALL)$BT)
Bcell<-rep(0,length(pData(ALL)$BT))
Bcell[grep("B",as.character(pData(ALL)$BT))]<-1


###################################################
### chunk number 7: BTBoot
###################################################
seed <- 99
BT.boot <- cache("BT.boot", MTP(X=filtX, Y=Bcell, 
alternative="greater", B=100, method="sd.minP", seed=seed))


###################################################
### chunk number 8: BTOut
###################################################
summary(BT.boot)


###################################################
### chunk number 9: BTGenes
###################################################
BT.diff <- BT.boot@adjp<=0.05
BT.AffyID <- geneNames(filtALL)[BT.diff]
mget(BT.AffyID[1:2], env=hgu95av2GENENAME)


###################################################
### chunk number 10: BTPlot
###################################################
par(mfrow=c(2,2))
plot(BT.boot)


###################################################
### chunk number 11: BTtppfp
###################################################
q <- c(0.05,0.1,0.25)
BT.tppfp <- fwer2tppfp(adjp=BT.boot@adjp, q=q)
comp.tppfp <- cbind(BT.boot@adjp, BT.tppfp)
mtps <- c("FWER",paste("TPPFP(",q,")", sep=""))
mt.plot(adjp=comp.tppfp, teststat=BT.boot@statistic, proc=mtps, 
leg=c(0.1,430), col=1:4, lty=1:4, lwd=3)
title("Comparison of TPPFP(q)-controlling AMTPs\n based on SD minP MTP")


###################################################
### chunk number 12: BTfdr
###################################################
BT.fdr <- fwer2fdr(adjp=BT.boot@adjp, method="both")$adjp
BT.marg.fdr <- mt.rawp2adjp(rawp=BT.boot@rawp, proc=c("BY","BH"))
comp.fdr <- cbind(BT.fdr, BT.marg.fdr$adjp[order(BT.marg.fdr$index),-1])
mtps <- c("AMTP Cons", "AMTP Rest", "BY", "BH")
mt.plot(adjp=comp.fdr, teststat=BT.boot@statistic, proc=mtps, 
leg=c(0.1,430), col=c(2,2,3,3), lty=rep(1:2,2), lwd=3)
title("Comparison of FDR-controlling MTPs")


###################################################
### chunk number 13: mbBoot
###################################################
BX<-filtX[,Bcell==1]
Bpheno<-pheno[Bcell==1,]
mb <- as.character(Bpheno$mol.biol)
table(mb)
other <- c("E2A/PBX1", "p15/p16")
mb.boot <- cache("mb.boot", MTP(X=BX[,!(mb%in%other)], Y=mb[!(mb%in%other)], 
test="f", alpha=c(0.01,0.1), B=100, get.cutoff=TRUE, seed=seed))


###################################################
### chunk number 14: mbOut
###################################################
mb.rej<-summary(mb.boot)$rejections


###################################################
### chunk number 15: mbOut
###################################################
mb.rej


###################################################
### chunk number 16: coxphPrep
###################################################
# Patients with original complete remission and who were followed up
cr.ind <- (Bpheno$remission=="CR")
cr.pheno <- Bpheno[cr.ind,]
times <- strptime(cr.pheno$"date last seen", "%m/%d/%Y")-strptime(cr.pheno$date.cr, "%m/%d/%Y")
time.ind <- !is.na(times)
times <- times[time.ind]
# Patients who haven't relapsed are treated as censored
cens <- ((1:length(times))%in%grep("CR", cr.pheno[time.ind,"f.u"]))
# Time to relapse
rel.times <- Surv(times, !cens)
patients <- (1: ncol(BX))[cr.ind][time.ind]
# Prepare data for MTP
relX <- BX[, patients]
relZ <- Bpheno[patients,]


###################################################
### chunk number 17: coxphBoot
###################################################
cox.boot <- cache("cox.boot", MTP(X=relX, Y=rel.times, Z=relZ,
       Z.incl="sex", Z.test=NULL, test="coxph.YvsXZ", B=100, get.cr=TRUE,
       seed=seed))


###################################################
### chunk number 18: coxphOut
###################################################
cox.diff <- cox.boot@adjp<=0.05
sum(cox.diff)
cox.AffyID <- geneNames(filtALL)[cox.diff]
mget(cox.AffyID, env=hgu95av2GENENAME)


###################################################
### chunk number 19: coxphPlot
###################################################
plot(cox.boot, which=5, top=5, sub.caption=NULL)
abline(h=0,col="red")


