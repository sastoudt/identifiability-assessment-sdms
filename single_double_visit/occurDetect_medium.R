load(file = "single_double_visit/raw_data/svdvOccurData_mediumNoSaturateLine.RData")

#### double visit ####
numParam <- 22 ## fix
resultsLogRegSpline <- matrix(NA, nrow = numIter, ncol = numParam)


basisPsi <- basisP <- rstart <- negLL <- counts <- value <- isConverged <- vector("list", numIter)

for (i in 1:numIter) {
  m <- 1:7
  -2.5 + m * 10.5 / 8

  basisPsi[[i]] <- basisPsiUse <- bs(zData[[i]][1:n], knots = -2.5 + m * 10.5 / 8)


  basisP[[i]] <- basisPUse <- bs(mData[[i]][1:n], knots = -2.5 + m * 10.5 / 8)



  ## organize data
  dataMat <- matrix(NA, n, 2)
  dataMat[, 1] <- obsData[[i]]
  dataMat[, 2] <- obsData2[[i]]

  umf <- unmarkedFrameOccu(y = dataMat)


  siteCovs(umf) <- data.frame(sitevar1 = basisPsi[[i]][, 1], sitevar2 = basisPsi[[i]][, 2], sitevar3 = basisPsi[[i]][, 3], sitevar4 = basisPsi[[i]][, 4], sitevar5 = basisPsi[[i]][, 5], sitevar6 = basisPsi[[i]][, 6], sitevar7 = basisPsi[[i]][, 7], sitevar8 = basisPsi[[i]][, 8], sitevar9 = basisPsi[[i]][, 9], sitevar10 = basisPsi[[i]][, 10])

  obsCovs(umf) <- data.frame(obsCov1 = rep(basisP[[i]][, 1], each = 2), obsCov2 = rep(basisP[[i]][, 2], each = 2), obsCov3 = rep(basisP[[i]][, 3], each = 2), obsCov4 = rep(basisP[[i]][, 4], each = 2), obsCov5 = rep(basisP[[i]][, 5], each = 2), obsCov6 = rep(basisP[[i]][, 6], each = 2), obsCov7 = rep(basisP[[i]][, 7], each = 2), obsCov8 = rep(basisP[[i]][, 8], each = 2), obsCov9 = rep(basisP[[i]][, 9], each = 2), obsCov10 = rep(basisP[[i]][, 10], each = 2))

  rstart[[i]] <- rnorm(numParam, 0, 2)
  ## fit
  fm <- occu(~ obsCov1 + obsCov2 + obsCov3 + obsCov4 + obsCov5 + obsCov6 + obsCov7 + obsCov8 + obsCov9 + obsCov10 ~ sitevar1 + sitevar2 + sitevar3 + sitevar4 + sitevar5 + sitevar6 + sitevar7 + sitevar8 + sitevar9 + sitevar10, umf, starts = rstart[[i]], control = list(maxit = 50000), se = F)




  resultsLogRegSpline[i, ] <- coef(fm)

  negLL[[i]] <- attr(fm, "negLogLike")
  counts[[i]] <- attr(fm, "opt")$counts
  value[[i]] <- attr(fm, "opt")$value
  isConverged[[i]] <- attr(fm, "opt")$convergence



  print(i)
}

save(resultsLogRegSpline, negLL, counts, value, isConverged, basisPsi, basisP, file = "single_double_visit/results/dvOccurSplineResults_mediumNoSaturateLine.RData")

#### Figure 4g ####

load(file = "single_double_visit/results/dvOccurSplineResults_mediumNoSaturateLine.RData")


idx <- order(zData[[1]][1:n])
jpeg("paper_figures/fig4g.jpeg", width = 824, height = 634, units = "px")

plot((zData[[1]][idx])[thin], (psiData[[1]][idx])[thin], ylim = c(0, 1), type = "l", lwd = .5, xlim = c(-3, 8), xlab = "x", ylab = "", cex.axis = 2, cex.lab = 2) # ,main="Double Visit \n Misspecified Logistic Regressions \n Spline",cex.main=1.5,sub="n = 100,000")
## true
for (i in 2:numIter) {
  idx <- order(zData[[i]][1:n])
  lines((zData[[i]][idx])[thin], (psiData[[i]][idx])[thin], lwd = .5)
}

## predicted var
for (i in 1:numIter) {
  fitted.val <- exp(resultsLogRegSpline[i, 1] +
    resultsLogRegSpline[i, 2] * basisPsi[[i]][, 1] +
    resultsLogRegSpline[i, 3] * basisPsi[[i]][, 2] +
    resultsLogRegSpline[i, 4] * basisPsi[[i]][, 3] +
    resultsLogRegSpline[i, 5] * basisPsi[[i]][, 4] +
    resultsLogRegSpline[i, 6] * basisPsi[[i]][, 5] +
    resultsLogRegSpline[i, 7] * basisPsi[[i]][, 6]
    + resultsLogRegSpline[i, 8] * basisPsi[[i]][, 7]
    + resultsLogRegSpline[i, 9] * basisPsi[[i]][, 8]
    + resultsLogRegSpline[i, 10] * basisPsi[[i]][, 9]
    + resultsLogRegSpline[i, 11] * basisPsi[[i]][, 10]) / (1 + exp(resultsLogRegSpline[i, 1] +
    resultsLogRegSpline[i, 2] * basisPsi[[i]][, 1] +
    resultsLogRegSpline[i, 3] * basisPsi[[i]][, 2] +
    resultsLogRegSpline[i, 4] * basisPsi[[i]][, 3] +
    resultsLogRegSpline[i, 5] * basisPsi[[i]][, 4] +
    resultsLogRegSpline[i, 6] * basisPsi[[i]][, 5] +
    resultsLogRegSpline[i, 7] * basisPsi[[i]][, 6] +
    resultsLogRegSpline[i, 8] * basisPsi[[i]][, 7]
    + resultsLogRegSpline[i, 9] * basisPsi[[i]][, 8]
    + resultsLogRegSpline[i, 10] * basisPsi[[i]][, 9]
    + resultsLogRegSpline[i, 11] * basisPsi[[i]][, 10]))
  idx <- order(zData[[i]][1:n])
  lines((zData[[i]][idx]), (fitted.val[idx]), col = "red", cex = .5)
}

idx <- order(zData[[1]][1:n])
lines(zData[[1]][idx], psiData[[1]][idx])


mtext(
  text = expression(paste(psi, "(x)")),
  side = 2,
  line = 2.5, cex = 2
)

dev.off()

#### Figure 4o ####

idx <- order(mData[[1]][1:n])
jpeg("paper_figures/fig4o.jpeg", width = 824, height = 634, units = "px")

plot((mData[[1]][idx])[thin], (pData[[1]][idx])[thin], ylim = c(0, 1), type = "l", lwd = .5, xlim = c(-3, 8), xlab = "z", ylab = "", cex.axis = 2, cex.lab = 2) # ,main="Double Visit \n Misspecified Logistic Regressions \n Spline",cex.main=1.5,sub="n = 100,000")
## true
for (i in 2:numIter) {
  idx <- order(mData[[i]][1:n])
  lines((mData[[i]][idx])[thin], (pData[[i]][idx])[thin], lwd = .5)
}
## predicted value
for (i in 1:numIter) {
  fitted.val <- exp(resultsLogRegSpline[i, 12] +
    resultsLogRegSpline[i, 13] * basisP[[i]][, 1] +
    resultsLogRegSpline[i, 14] * basisP[[i]][, 2] +
    resultsLogRegSpline[i, 15] * basisP[[i]][, 3] +
    resultsLogRegSpline[i, 16] * basisP[[i]][, 4] +
    resultsLogRegSpline[i, 17] * basisP[[i]][, 5] +
    resultsLogRegSpline[i, 18] * basisP[[i]][, 6]

    + resultsLogRegSpline[i, 19] * basisP[[i]][, 7]
    + resultsLogRegSpline[i, 20] * basisP[[i]][, 8]
    + resultsLogRegSpline[i, 21] * basisP[[i]][, 9]
    + resultsLogRegSpline[i, 22] * basisP[[i]][, 10]) / (1 +
    exp(resultsLogRegSpline[i, 12] +
      resultsLogRegSpline[i, 13] * basisP[[i]][, 1] +
      resultsLogRegSpline[i, 14] * basisP[[i]][, 2] +
      resultsLogRegSpline[i, 15] * basisP[[i]][, 3] +
      resultsLogRegSpline[i, 16] * basisP[[i]][, 4] +
      resultsLogRegSpline[i, 17] * basisP[[i]][, 5] +
      resultsLogRegSpline[i, 18] * basisP[[i]][, 6]

      + resultsLogRegSpline[i, 19] * basisP[[i]][, 7]
      + resultsLogRegSpline[i, 20] * basisP[[i]][, 8]
      + resultsLogRegSpline[i, 21] * basisP[[i]][, 9]
      + resultsLogRegSpline[i, 22] * basisP[[i]][, 10]))
  idx <- order(mData[[i]][1:n])
  lines((mData[[i]][idx])[thin], (fitted.val[idx])[thin], col = "red", cex = .5)
}


lines(sort(mData[[1]][thin]), sort(pData[[1]][thin]))

mtext(
  text = "p(z)",
  side = 2,
  line = 2.5, cex = 2
)
dev.off()


#### single visit ####

resultsLeleSpline <- matrix(NA, nrow = numIter, ncol = numParam)

basisPsi <- basisP <- rstart <- counts <- loglik <- isConverged <- vector("list", numParam)
ptm <- proc.time()

for (i in 1:25) {
  m <- 1:7
  -2.5 + m * 10.5 / 8

  basisPsi[[i]] <- basisPsiUse <- bs(zData[[i]], knots = -2.5 + m * 10.5 / 8)


  basisP[[i]] <- basisPUse <- bs(mData[[i]], knots = -2.5 + m * 10.5 / 8)



  rstart[[i]] <- rnorm(numParam, 0, 2)
  mod <- svocc2(obsDataM[[i]] ~ basisPsi[[i]][, 1] + basisPsi[[i]][, 2] + basisPsi[[i]][, 3] + basisPsi[[i]][, 4] + basisPsi[[i]][, 5] + basisPsi[[i]][, 6] + basisPsi[[i]][, 7] + basisPsi[[i]][, 8] + basisPsi[[i]][, 9] + basisPsi[[i]][, 10] | basisP[[i]][, 1] + basisP[[i]][, 2] + basisP[[i]][, 3] + basisP[[i]][, 4] + basisP[[i]][, 5] + basisP[[i]][, 6] + basisP[[i]][, 7] + basisP[[i]][, 8] + basisP[[i]][, 9] + basisP[[i]][, 10], penalized = F, method = "optim", inits = rstart[[i]])




  resultsLeleSpline[i, ] <- coef(mod)
  counts[[i]] <- mod$count
  loglik[[i]] <- mod$loglik
  isConverged[[i]] <- mod$converged[1]



  print(i)
}

save(resultsLeleSpline, rstart, basisPsi, basisP, counts, loglik, isConverged, file = "single_double_visit/results/svOccurSplineResults_mediumNoSaturateLine.RData")

#### Figure 4c ####
load(file = "single_double_visit/results/svOccurSplineResults_mediumNoSaturateLine.RData")

idx <- order(zData[[1]])
jpeg("paper_figures/fig4c.jpeg", width = 824, height = 634, units = "px")

plot((zData[[1]][idx])[thin], (psiData[[1]][idx])[thin], ylim = c(0, 1), lwd = .5, xlim = c(-3, 8), xlab = "x", ylab = "", cex.axis = 2, cex.lab = 2, type = "l") # ,main="Single Visit \n Misspecified Logistic Regressions \n Spline",cex.main=1.5,sub="n = 100,000")
## true
for (i in 2:numIter) {
  idx <- order(zData[[i]])
  lines((zData[[i]][idx])[thin], (psiData[[i]][idx])[thin], lwd = .5)
}

## predicted value
for (i in 1:numIter) {
  fitted.val <- exp(resultsLeleSpline[i, 1] + resultsLeleSpline[i, 2] * basisPsi[[i]][, 1] + resultsLeleSpline[i, 3] * basisPsi[[i]][, 2] + resultsLeleSpline[i, 4] * basisPsi[[i]][, 3] +
    resultsLeleSpline[i, 5] * basisPsi[[i]][, 4] + resultsLeleSpline[i, 6] * basisPsi[[i]][, 5] + resultsLeleSpline[i, 7] * basisPsi[[i]][, 6]
    + resultsLeleSpline[i, 8] * basisPsi[[i]][, 7]
    + resultsLeleSpline[i, 9] * basisPsi[[i]][, 8]
    + resultsLeleSpline[i, 10] * basisPsi[[i]][, 9]
    + resultsLeleSpline[i, 11] * basisPsi[[i]][, 10]) / (1 + exp(resultsLeleSpline[i, 1] + resultsLeleSpline[i, 2] * basisPsi[[i]][, 1] + resultsLeleSpline[i, 3] * basisPsi[[i]][, 2] + resultsLeleSpline[i, 4] * basisPsi[[i]][, 3] +
    resultsLeleSpline[i, 5] * basisPsi[[i]][, 4] + resultsLeleSpline[i, 6] * basisPsi[[i]][, 5] + resultsLeleSpline[i, 7] * basisPsi[[i]][, 6] + resultsLeleSpline[i, 8] * basisPsi[[i]][, 7]
    + resultsLeleSpline[i, 9] * basisPsi[[i]][, 8]
    + resultsLeleSpline[i, 10] * basisPsi[[i]][, 9]
    + resultsLeleSpline[i, 11] * basisPsi[[i]][, 10]))
  idx <- order(zData[[i]])
  lines((zData[[i]][idx])[thin], (fitted.val[idx])[thin], col = "red", cex = .5)
}



lines(sort(zData[[1]][thin]), sort(psiData[[1]][thin]))

mtext(
  text = expression(paste(psi, "(x)")),
  side = 2,
  line = 2.5, cex = 2
)
dev.off()

#### Figure 4k ####

idx <- order(mData[[1]])
jpeg("paper_figures/fig4k.jpeg", width = 824, height = 634, units = "px")

plot((mData[[1]][idx])[thin], (pData[[1]][idx])[thin], ylim = c(0, 1), lwd = .5, xlim = c(-3, 8), xlab = "z", ylab = "", cex.axis = 2, cex.lab = 2, type = "l") # ,main="Single Visit \n Misspecified Logistic Regressions \n Spline",cex.main=1.5,sub="n = 100,000")
## true
for (i in 2:numIter) {
  idx <- order(mData[[i]])
  lines((mData[[i]][idx])[thin], (pData[[i]][idx])[thin], lwd = .5)
}

## predicted value
for (i in 1:numIter) {
  fitted.val <- exp(resultsLeleSpline[i, 12] + resultsLeleSpline[i, 13] * basisP[[i]][, 1] + resultsLeleSpline[i, 14] * basisP[[i]][, 2] + resultsLeleSpline[i, 15] * basisP[[i]][, 3] +
    resultsLeleSpline[i, 16] * basisP[[i]][, 4] + resultsLeleSpline[i, 17] * basisP[[i]][, 5] + resultsLeleSpline[i, 18] * basisP[[i]][, 6]
    + resultsLeleSpline[i, 19] * basisP[[i]][, 7]
    + resultsLeleSpline[i, 20] * basisP[[i]][, 8]
    + resultsLeleSpline[i, 21] * basisP[[i]][, 9]
    + resultsLeleSpline[i, 22] * basisP[[i]][, 10]) / (1 + exp(resultsLeleSpline[i, 12] + resultsLeleSpline[i, 13] * basisP[[i]][, 1] + resultsLeleSpline[i, 14] * basisP[[i]][, 2] + resultsLeleSpline[i, 15] * basisP[[i]][, 3] +
    resultsLeleSpline[i, 16] * basisP[[i]][, 4] + resultsLeleSpline[i, 17] * basisP[[i]][, 5] + resultsLeleSpline[i, 18] * basisP[[i]][, 6]
    + resultsLeleSpline[i, 19] * basisP[[i]][, 7]
    + resultsLeleSpline[i, 20] * basisP[[i]][, 8]
    + resultsLeleSpline[i, 21] * basisP[[i]][, 9]
    + resultsLeleSpline[i, 22] * basisP[[i]][, 10]))
  idx <- order(mData[[i]])
  lines((mData[[i]][idx])[thin], (fitted.val[idx])[thin], col = "red", cex = .5)
}


lines(sort(mData[[1]][thin]), sort(pData[[1]][thin]))

mtext(
  text = "p(z)",
  side = 2,
  line = 2.5, cex = 2
)
dev.off()