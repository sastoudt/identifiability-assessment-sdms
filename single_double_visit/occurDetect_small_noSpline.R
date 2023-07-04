#### double visit ####

load("single_double_visit/raw_data/svdvOccurData_smallNoSaturateLine.RData")


numParam <- 4
resultsLogReg <- matrix(NA, nrow = numIter, ncol = numParam)
negLL <- counts <- value <- isConverged <- rstart <- vector("list", numIter)

for (i in 1:25) {

  ## get data in right format
  dataMat <- matrix(NA, n, 2)
  dataMat[, 1] <- obsData[[i]]
  dataMat[, 2] <- obsData2[[i]]

  umf <- unmarkedFrameOccu(y = dataMat)

  siteCovs(umf) <- data.frame(sitevar1 = zData[[i]][1:n])
  obsCovs(umf) <- data.frame(obsCovs = rep(mData[[i]][1:n], each = 2))

  rstart[[i]] <- rnorm(numParam, 0, 2)
  ## fit model
  fm <- occu(~obsCovs ~ sitevar1, umf, starts = rstart[[i]])

  ## optim output
  negLL[[i]] <- attr(fm, "negLogLike")
  counts[[i]] <- attr(fm, "opt")$counts
  value[[i]] <- attr(fm, "opt")$value
  isConverged[[i]] <- attr(fm, "opt")$convergence

  ## estimates
  resultsLogReg[i, ] <- coef(fm)

  print(i)
}

save(resultsLogReg, rstart, negLL, counts, value, isConverged, file = "single_double_visit/results/dvOccurResults_smallNoSaturateLine.RData")

#### Figure 4e ####
load(file = "single_double_visit/results/dvOccurResults_smallNoSaturateLine.RData")

thin <- 1:n ## no need to thin with small data, helps later on when data gets larger

# can sort in the non-spline case instead of order
jpeg("paper_figures/fig4e.jpeg", width = 824, height = 634, units = "px")
plot(sort(zData[[1]][thin]), sort(psiData[[1]][thin]), ylim = c(0, 1), type = "l", lwd = .5, xlim = c(-3, 8), xlab = "x", ylab = "", cex.axis = 2, cex.lab = 2)
## true values
for (i in 2:numIter) {
  lines(sort(zData[[i]][thin]), sort(psiData[[i]][thin]), lwd = .5)
}

## predicted values
for (i in 1:numIter) {
  idx <- order(zData[[i]])
  lines(zData[[i]][idx], exp(resultsLogReg[i, 1] + resultsLogReg[i, 2] * zData[[i]][idx]) / (1 + exp(resultsLogReg[i, 1] + resultsLogReg[i, 2] * zData[[i]][idx])), col = "red", cex = .5)
}
lines(sort(zData[[1]][thin]), sort(psiData[[1]][thin]), lwd = .5)

mtext(
  text = expression(paste(psi, "(x)")),
  side = 2,
  line = 2.5, cex = 2
)
dev.off()

#### Figure 4m ####
jpeg("paper_figures/fig4m.jpeg", width = 824, height = 634, units = "px")
plot(sort(mData[[1]][thin]), sort(pData[[1]][thin]), ylim = c(0, 1), type = "l", lwd = .5, xlim = c(-3, 8), xlab = "z", ylab = "", cex.axis = 2, cex.lab = 2)
## truth
for (i in 2:numIter) {
  lines(sort(mData[[i]][thin]), sort(pData[[i]][thin]), lwd = .5)
}

## predicted
for (i in 1:numIter) {
  idx <- order(mData[[i]])
  lines(mData[[i]][idx], exp(resultsLogReg[i, 3] + resultsLogReg[i, 4] * mData[[i]][idx]) / (1 + exp(resultsLogReg[i, 3] + resultsLogReg[i, 4] * mData[[i]][idx])), col = "red", cex = .5)
}
lines(sort(mData[[1]][thin]), sort(pData[[1]][thin]), lwd = .5)

mtext(
  text = "p(z)",
  side = 2,
  line = 2.5, cex = 2
)
dev.off()

#### single visit ####

## had to modify source code to get extra information from the optimization
## more information about modifications in individual files
source("modified_source_code/svocc.fit.R")
source("modified_source_code/svocc.R")
source("modified_source_code/svisitFormula.R")

# small dependency
# https://github.com/psolymos/detect/blob/master/R/solveneardot.R
.solvenear <- function(x) {
  xinv <- try(solve(x), silent = TRUE)
  if (inherits(xinv, "try-error")) {
    xinv <- as.matrix(solve(Matrix::nearPD(x)$mat))
  }
  xinv
}

resultsLele <- matrix(NA, nrow = numIter, ncol = numParam)
counts <- loglik <- isConverged <- rstart <- vector("list", numIter)

for (i in 1:numIter) {
  rstart[[i]] <- rnorm(numParam, 0, 2)
  mod <- svocc2(obsDataM[[i]] ~ zData[[i]] | mData[[i]], penalized = F, method = "optim", inits = rstart[[i]])
  counts[[i]] <- mod$count
  loglik[[i]] <- mod$loglik
  isConverged[[i]] <- mod$converged[1]

  resultsLele[i, ] <- coef(mod)


  print(i)
}

save(resultsLele, counts, loglik, isConverged, rstart, file = "single_double_visit/results/svOccurResults_smallNoSaturateLine.RData")

#### Figure 4a ####

load(file = "single_double_visit/results/svOccurResults_smallNoSaturateLine.RData")


idx <- order(zData[[1]])
jpeg("paper_figures/fig4a.jpeg", width = 824, height = 634, units = "px")

plot(zData[[1]][idx], psiData[[1]][idx], ylim = c(0, 1), lwd = .5, xlim = c(-3, 8), xlab = "x", ylab = "", cex.axis = 2, cex.lab = 2, type = "l")
## true
for (i in 2:numIter) {
  idx <- order(zData[[i]])
  lines(zData[[i]][idx], psiData[[i]][idx], lwd = .5)
}

## predicted
for (i in 1:numIter) {
  idx <- order(zData[[i]])
  lines(zData[[i]][idx], exp(resultsLele[i, 1] + resultsLele[i, 2] * zData[[i]][idx]) / (1 + exp(resultsLele[i, 1] + resultsLele[i, 2] * zData[[i]][idx])), col = "red", cex = .5)
}



mtext(
  text = expression(paste(psi, "(x)")),
  side = 2,
  line = 2.5, cex = 2
)
dev.off()

#### Figure 4i ####

idx <- order(mData[[1]])
jpeg("paper_figures/fig4i.jpeg", width = 824, height = 634, units = "px")

plot(mData[[1]][idx], pData[[1]][idx], ylim = c(0, 1), lwd = .5, xlim = c(-3, 8), xlab = "z", ylab = "", cex.axis = 2, cex.lab = 2, type = "l") # ,main="Single Visit \n Misspecified Logistic Regressions",cex.main=1.5)
## true
for (i in 2:numIter) {
  idx <- order(mData[[i]])
  lines(mData[[i]][idx], pData[[i]][idx], lwd = .5)
}


## predicted
for (i in 1:numIter) {
  idx <- order(mData[[i]])
  lines(mData[[i]][idx], exp(resultsLele[i, 3] + resultsLele[i, 4] * mData[[i]][idx]) / (1 + exp(resultsLele[i, 3] + resultsLele[i, 4] * mData[[i]][idx])), col = "red", cex = .5)
}

mtext(
  text = "p(z)",
  side = 2,
  line = 2.5, cex = 2
)
dev.off()