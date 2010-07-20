###################################################
### chunk number 1: RbcBook1
###################################################
library("RbcBook1")


###################################################
### chunk number 2: bagging
###################################################
simple_bagging <- function(x, lsample, M = 100, nu = 0.5) {
    I <- nrow(lsample)
    bsample <- rmultinom(M, I, rep(1, I)/I)
    pred <- rep(0, nrow(x))
    rpc <- rpart.control(xval = 0, cp = 0.01)
    for (m in 1:M) {
        weaktree <- rpart(y ~ ., data = lsample, weights = bsample[,m],
                          control = rpc)
        prtree <- predict(weaktree, newdata = x, type = "class")
        pred <- pred + (prtree == levels(lsample$y)[2])
    }
    factor(pred / M > nu, levels = levels(lsample$y))
}


###################################################
### chunk number 3: ALLpkgs
###################################################
pkgload <- require(Biobase, quietly = TRUE)
pkgload <- require(ALL, quietly = TRUE)
pkgload <- require(MLInterfaces, quietly = TRUE)


###################################################
### chunk number 4: ALLbegin
###################################################
data(ALL)


###################################################
### chunk number 5: ALLBcell
###################################################
cvv <- apply(exprs(ALL), 1, function(x) sd(log(x)))  
ok <- cvv > .08 & cvv < .18
BStagelev <- paste("B", 1:4, sep="")
BALL <- ALL[ok, ALL$BT %in% BStagelev]
pData(phenoData(BALL))$BStage <- factor(BALL$BT, levels = BStagelev)


###################################################
### chunk number 6: ALLtable
###################################################
table(BALL$BStage)


###################################################
### chunk number 7: ALLfilter
###################################################
response <- BALL$BStage
I <- ncol(exprs(BALL))
expressions <- t(apply(exprs(BALL), 1, rank))
Iindx <- 1:I

var_selection <- function(indx, expressions, response, p = 100) {

    y <- switch(class(response),
        "factor" = { model.matrix(~ response - 1)[indx, ,drop = FALSE] },
        "Surv" = { matrix(cscores(response[indx]), ncol = 1) },
        "numeric" = { matrix(rank(response[indx]), ncol = 1) }
    )

    x <- expressions[,indx, drop = FALSE]
    n <- nrow(y)
    linstat <- x %*% y
    Ey <- matrix(colMeans(y), nrow = 1)
    Vy <- matrix(rowMeans((t(y) - as.vector(Ey))^2), nrow = 1)

    rSx <- matrix(rowSums(x), ncol = 1)   
    rSx2 <- matrix(rowSums(x^2), ncol = 1)
    E <- rSx %*% Ey
    V <- n / (n - 1) * kronecker(Vy, rSx2)
    V <- V - 1 / (n - 1) * kronecker(Vy, rSx^2)

    stats <- abs(linstat - E) / sqrt(V)
    stats <- do.call("pmax", as.data.frame(stats))
    return(which(stats > sort(stats)[length(stats) - p]))
}
selected <- var_selection(Iindx, expressions, response)


###################################################
### chunk number 8: ALLrandomForest
###################################################
rf <- randomForestB(BALL[selected,], "BStage", Iindx[-1],
                    sampsize = I - 1)
print(rf)


###################################################
### chunk number 9: ALLbenchmark eval=FALSE
###################################################
## set.seed(290875)
## B <- 100
## performance <- as.data.frame(matrix(0, nrow = B, ncol = 4))
## colnames(performance) <- c("RF", "Bagg", "LBoost", "Guess")
## for (b in 1:B) {
##      bsample <- sample(Iindx, I, replace = TRUE)
##      selected <- var_selection(bsample, expressions, response)
##      rf3 <- randomForestB(BALL[selected,], "BStage", bsample, 
##                           mtry = 3, sampsize = I)
##      predicted3 <- factor(rf3@predLabels, levels = levels(response))
##      performance[b, 1] <- mean(response[-bsample] != predicted3)
## 
##      rfBagg <- randomForestB(BALL[selected,], "BStage", bsample,
##                              mtry = length(selected), 
##                              sampsize = I)
##      predictedBagg <- factor(rfBagg@predLabels, levels = levels(response))
##      performance[b, 2] <- mean(response[-bsample] != predictedBagg)
## 
##      lb <- logitboostB(BALL[selected,], "BStage", bsample, 100)
##      predictedlb <- factor(lb@predLabels, levels = levels(response))
##      performance[b, 3] <- mean(response[-bsample] != predictedlb)   
## 
##      performance[b, 4] <- mean(response[-bsample] !=
##          levels(response)[which.max(tabulate(response[bsample]))])
## 
##      #save(performance, file = "ALLperformance.Rda")
## }


###################################################
### chunk number 10: ALLload
###################################################
data(ALLperformance)
B <- nrow(performance)


###################################################
### chunk number 11: ALL-KW
###################################################
friedman.test(as.matrix(performance))


###################################################
### chunk number 12: ALLinference-plot1
###################################################
matplot(1:ncol(performance), t(performance), 
        xlab="", ylab = "Misclassification error",
        type = "l", col = "#377EB8", lty = 1, axes = FALSE)
axis(1, at = 1:ncol(performance), labels = colnames(performance))
axis(2)
box()


###################################################
### chunk number 13: ALLinference-plot2
###################################################
boxplot(performance, col="#4DAF4A")


###################################################
### chunk number 14: ALLci
###################################################
t.test(performance$RF, performance$Bagg, 
       paired = TRUE, conf.int = TRUE)$conf.int


###################################################
### chunk number 15: pkgs
###################################################
pkgload <- require(randomForest, quietly = TRUE)
pkgload <- require(MLInterfaces, quietly = TRUE)
pkgload <- require(MASS, quietly = TRUE)
# pkgload <- require(exactRankTests, quietly = TRUE)
pkgload <- require(Biobase, quietly = TRUE)
pkgload <- require(kidpack, quietly = TRUE)


###################################################
### chunk number 16: KPbegin
###################################################
set.seed(290875)
data(eset)
pData(phenoData(eset))$type <- as.factor(eset$type)


###################################################
### chunk number 17: KPclass
###################################################
table(pData(phenoData(eset))$type)


###################################################
### chunk number 18: KPfilterKW
###################################################
response <- eset$type
expressions <- t(apply(exprs(eset), 1, rank))
I <- ncol(exprs(eset))
Iindx <- 1:I
selected <- var_selection(Iindx, expressions, response)


###################################################
### chunk number 19: KPrandomForest
###################################################
rf <- randomForestB(eset[selected,], "type", Iindx[-1],
                    sampsize = I - 1)
rf


###################################################
### chunk number 20: KPbenchmark eval=FALSE
###################################################
## B <- 100
## performance <- as.data.frame(matrix(0, nrow = B, ncol = 4))
## colnames(performance) <- c("RF", "Bagg", "LBoost", "Guess")
## for (b in 1:B) {
##      bsample <- sample(Iindx, I, replace = TRUE)
##      selected <- var_selection(bsample, expressions, response)
##      rf3 <- randomForestB(eset[selected,], "type", bsample,   
##                           mtry = 3, sampsize = I)
##      predicted3 <- factor(rf3@predLabels, levels = levels(response))
##      performance[b, 1] <- mean(response[-bsample] != predicted3)
## 
##      rfBagg <- randomForestB(eset[selected,], "type", bsample,
##                              mtry = length(selected), 
##                              sampsize = I)
##      predictedBagg <- factor(rfBagg@predLabels, levels = levels(response))
##      performance[b, 2] <- mean(response[-bsample] != predictedBagg)
## 
##      lb <- logitboostB(eset[selected,], "type", bsample, 100)
##      predictedlb <- factor(lb@predLabels, levels = levels(response))
##      performance[b, 3] <- mean(response[-bsample] != predictedlb)   
## 
##      performance[b, 4] <- mean(response[-bsample] !=
##          levels(response)[which.max(tabulate(response[bsample]))])
## }


###################################################
### chunk number 21: KPload
###################################################
data(kidpackperformance)
B <- nrow(performance)


###################################################
### chunk number 22: KPinference
###################################################
friedman.test(as.matrix(performance))


###################################################
### chunk number 23: KPclass-plot
###################################################
par(mfrow = c(1, 2))
matplot(1:ncol(performance), t(performance),
        type = "l", col = "#377EB8", lty = 1, xlab="", ylab = "Misclassification error",
        axes = FALSE)
        axis(1, at = 1:ncol(performance), labels = colnames(performance))
        axis(2)
        box()
boxplot(performance, col="#4DAF4A")


###################################################
### chunk number 24: KPclass-wilcox
###################################################
t.test(performance$RF, performance$Bagg,
             paired = TRUE, conf.int = TRUE)$conf.int
t.test(performance$RF, performance$LBoost,
             paired = TRUE, conf.int = TRUE)$conf.int


###################################################
### chunk number 25: pkgload
###################################################
pkgload <- require(exactRankTests, quietly = TRUE)
pkgload <- require(kidpack, quietly = TRUE)
pkgload <- require(ipred, quietly = TRUE)
pkgload <- require(survival, quietly = TRUE)
data(eset)
set.seed(290875)


###################################################
### chunk number 26: kidpackSurv
###################################################
remove <- is.na(eset$survival.time)
seset <- eset[,!remove]
response <- Surv(seset$survival.time, seset$died)
response[response[,1] == 0] <- 1
expressions <- t(apply(exprs(seset), 1, rank))
exprDF <- as.data.frame(t(expressions))

I <- nrow(exprDF)
Iindx <- 1:I


###################################################
### chunk number 27: kidpackSurvSelect
###################################################
selected <- var_selection(Iindx, expressions, response)


###################################################
### chunk number 28: kidpackSurvbagg
###################################################
bagg <- bagging(response ~., data = exprDF[,selected], 
                ntrees = 100)
prKM <- predict(bagg)
sbrier(response, prKM)
sbrier(response, survfit(response))


###################################################
### chunk number 29: kidpackSurvplot
###################################################
plot(survfit(response), lwd = 4, conf.int = FALSE, xlab =
     "Survival time in month", ylab = "Probability")
col <- c("lightgray", "darkblue", "red3")
type <- factor(seset$type)
table(type)
for (i in 1:length(prKM))
    lines(prKM[[i]], lty = 2, col = col[as.numeric(type)[i]])


###################################################
### chunk number 30: KSbenchmark eval=FALSE
###################################################
## set.seed(290875)
## B <- 100 
## performance <- as.data.frame(matrix(0, nrow = B, ncol = 2))
## colnames(performance) <- c("Bagging", "Kaplan-Meier") 
## for (b in 1:B) {
##     bsample <- sample(Iindx, I, replace = TRUE)
##     selected <- var_selection(bsample, expressions, response)
##     bagg <- bagging(response ~., data = exprDF[,selected],   
##                     subset = bsample, ntrees = 100)
##     pr <- predict(bagg, newdata = exprDF[-bsample,])
##     KM <- survfit(response[bsample])
##     performance[b, 1] <- sbrier(response[-bsample], pr)
##     performance[b, 2] <- sbrier(response[-bsample], KM)
## }


###################################################
### chunk number 31: KSload
###################################################
data(Survperformance)
B <- nrow(performance)


###################################################
### chunk number 32: kidpackSurv-inference-plot
###################################################
par(mfrow = c(1, 2))
matplot(1:ncol(performance), t(performance),
        type = "l", col = "#377EB8", lty = 1, xlab="", ylab = "Brier scores",
        axes = FALSE)
        axis(1, at = 1:ncol(performance), labels = colnames(performance))
        axis(2)
        box()
boxplot(performance, col="#4DAF4A")


