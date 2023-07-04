load(file = "single_double_visit_abundance/raw_data/abunData_medium.RData")

#### double visit ####
numParam <- 14
resultsLogRegSpline <- matrix(NA, nrow = numIter, ncol = numParam)
basisPsi <- basisP <- varCovar <- vector("list", numIter)

for (i in 1:numIter) {
  m <- 1:7

  basisPsi[[i]] <- bs(zData[[i]][1:n], knots = -2.5 + m * 10.5 / 8)
  basisPsiUse <- basisPsi[[i]]
  basisP[[i]] <- bs(mData[[i]][1:n], knots = -2.5 + m * 10.5 / 8)
  basisPUse <- basisP[[i]]

  dataMat <- matrix(NA, n, 2)
  dataMat[, 1] <- obsData[[i]][1:n]
  dataMat[, 2] <- obsData2[[i]]
  visitMat <- visitMat2 <- visitMat3 <- visitMat4 <- visitMat5 <- visitMat6 <- matrix(NA, n, 2)

  visitMat[, 1] <- visitMat[, 2] <- basisP[[i]][, 1]
  visitMat2[, 1] <- visitMat2[, 2] <- basisP[[i]][, 2]
  visitMat3[, 1] <- visitMat3[, 2] <- basisP[[i]][, 3]
  visitMat4[, 1] <- visitMat4[, 2] <- basisP[[i]][, 4]
  visitMat5[, 1] <- visitMat5[, 2] <- basisP[[i]][, 5]
  visitMat6[, 1] <- visitMat6[, 2] <- basisP[[i]][, 6]
  umf <- unmarkedFramePCount(y = dataMat, siteCovs = data.frame(x1 = basisPsi[[i]][, 1], x2 = basisPsi[[i]][, 2], x3 = basisPsi[[i]][, 3], x4 = basisPsi[[i]][, 4], x5 = basisPsi[[i]][, 5], x6 = basisPsi[[i]][, 6]), obsCovs = list(visit1 = visitMat, visit2 = visitMat2, visit3 = visitMat3, visit4 = visitMat4, visit5 = visitMat5, visit6 = visitMat6))
  test <- pcount(~ visit1 + visit2 + visit3 + visit4 + visit5 + visit6 ~ x1 + x2 + x3 + x4 + x5 + x6, umf, K = 50)

  resultsLogRegSpline[i, ] <- coef(test)

  print(i)
}

save(resultsLogRegSpline, basisP, basisPsi, file = "single_double_visit_abundance/results/dvAbundanceMedium.RData")

#### Figures S1g and S2g ####

load(file = "single_double_visit_abundance/raw_data/abunData_medium.RData")
load("single_double_visit_abundance/results/dvAbundanceMedium.RData")

thin <- 1:n
jpeg("paper_figures/figS1g.jpeg", width = 824, height = 634, units = "px")
plot(zData[[1]][thin], lambdaData[[1]][thin], lwd = .5, xlim = c(-3, 8), ylim = c(0, 25), xlab = "x", ylab = "", cex.axis = 2, cex.lab = 2)
for (i in 2:25) {
  points(zData[[i]][thin], lambdaData[[i]][thin], lwd = .5)
}


for (i in 1:25) {
  fitted.val <- exp(resultsLogRegSpline[i, 1] + resultsLogRegSpline[i, 2] * basisPsi[[i]][, 1] + resultsLogRegSpline[i, 3] * basisPsi[[i]][, 2] + resultsLogRegSpline[i, 4] * basisPsi[[i]][, 3] +
    resultsLogRegSpline[i, 5] * basisPsi[[i]][, 4] + resultsLogRegSpline[i, 6] * basisPsi[[i]][, 5] + resultsLogRegSpline[i, 7] * basisPsi[[i]][, 6])
  idx <- order(zData[[i]][1:n])
  points(zData[[i]][idx], fitted.val[idx], col = "red")
}

mtext(
  text = expression(paste(lambda, "(x)")),
  side = 2,
  line = 2, cex = 2
)
dev.off()

jpeg("paper_figures/figS2g.jpeg", width = 824, height = 634, units = "px")
plot(mData[[1]][thin], pData[[1]][thin], ylim = c(0, 1), lwd = .5, xlim = c(-3, 8), xlab = "z", ylab = "", cex.axis = 2, cex.lab = 2)
for (i in 2:25) {
  points(mData[[i]][thin], pData[[i]][thin], lwd = .5)
}


for (i in 1:25) {
  fitted.val <- exp(resultsLogRegSpline[i, 8] + resultsLogRegSpline[i, 9] * basisP[[i]][, 1] + resultsLogRegSpline[i, 10] * basisP[[i]][, 2] + resultsLogRegSpline[i, 11] * basisP[[i]][, 3] +
    resultsLogRegSpline[i, 12] * basisP[[i]][, 4] + resultsLogRegSpline[i, 13] * basisP[[i]][, 5] + resultsLogRegSpline[i, 14] * basisP[[i]][, 6]) / (1 + exp(resultsLogRegSpline[i, 8] + resultsLogRegSpline[i, 9] * basisP[[i]][, 1] + resultsLogRegSpline[i, 10] * basisP[[i]][, 2] + resultsLogRegSpline[i, 11] * basisP[[i]][, 3] +
    resultsLogRegSpline[i, 12] * basisP[[i]][, 4] + resultsLogRegSpline[i, 13] * basisP[[i]][, 5] + resultsLogRegSpline[i, 14] * basisP[[i]][, 6]))


  points(mData[[i]][thin], fitted.val[thin], col = "red")
}

mtext(
  text = "p(z)",
  side = 2,
  line = 2.5, cex = 2
)
dev.off()


#### single visit ####
resultsLeleSpline <- matrix(NA, nrow = numIter, ncol = numParam)

basisPsi <- basisP <- varCovar <- vector("list", numIter)

for (i in 1:numIter) {
  m <- 1:7

  basisPsi[[i]] <- basisPsiUse <- bs(zData[[i]], knots = -2.5 + m * 10.5 / 8)
  basisP[[i]] <- basisPUse <- bs(mData[[i]], knots = -2.5 + m * 10.5 / 8)

  mod <- svabu(obsData[[i]] ~ basisPsiUse[, 1] + basisPsiUse[, 2] + basisPsiUse[, 3] + basisPsiUse[, 4] + basisPsiUse[, 5] + basisPsiUse[, 6] | basisPUse[, 1] + basisPUse[, 2] + basisPUse[, 3] + basisPUse[, 4] + basisPUse[, 5] + basisPUse[, 6], dist = "P", zeroinfl = F)
  resultsLeleSpline[i, ] <- coef(mod)

  print(i)
}

save(resultsLeleSpline, basisPsi, basisP, file = "single_double_visit_abundance/results/svAbundanceMedium.RData")

#### Figures S1c and S2c ####

load(file = "single_double_visit_abundance/raw_data/abunData_medium.RData")
load("single_double_visit_abundance/results/svAbundanceMedium.RData")

idx <- order(zData[[1]])
jpeg("paper_figures/figS1c.jpeg", width = 824, height = 634, units = "px")
plot(zData[[1]][idx], lambdaData[[1]][idx], lwd = .5, xlim = c(-3, 8), xlab = "x", ylab = "", cex.axis = 2, cex.lab = 2, type = "l", ylim = c(0, 25))
for (i in 2:25) {
  idx <- order(zData[[i]])
  lines(zData[[i]][idx], lambdaData[[i]][idx], lwd = .5)
}


for (i in 1:25) {
  fitted.val <- exp(resultsLeleSpline[i, 1] + resultsLeleSpline[i, 2] * basisPsi[[i]][, 1] + resultsLeleSpline[i, 3] * basisPsi[[i]][, 2] + resultsLeleSpline[i, 4] * basisPsi[[i]][, 3] +
    resultsLeleSpline[i, 5] * basisPsi[[i]][, 4] + resultsLeleSpline[i, 6] * basisPsi[[i]][, 5] + resultsLeleSpline[i, 7] * basisPsi[[i]][, 6])
  idx <- order(zData[[i]])
  lines(zData[[i]][idx], fitted.val[idx], col = "red")
}


i <- 1
idx <- order(zData[[i]])
lines(zData[[i]][idx], lambdaData[[i]][idx], lwd = 2)


mtext(
  text = expression(paste(lambda, "(x)")),
  side = 2,
  line = 2, cex = 2
)
dev.off()

jpeg("paper_figures/figS2c.jpeg", width = 824, height = 634, units = "px")
idx <- order(mData[[1]])
plot(mData[[1]][idx], pData[[1]][idx], ylim = c(0, 1), lwd = .5, xlim = c(-3, 8), xlab = "z", ylab = "", cex.axis = 2, cex.lab = 2, type = "l")
for (i in 2:25) {
  idx <- order(mData[[i]])
  lines(mData[[i]][idx], pData[[i]][idx], lwd = .5)
}


for (i in 1:25) {
  fitted.val <- exp(resultsLeleSpline[i, 8] + resultsLeleSpline[i, 9] * basisP[[i]][, 1] + resultsLeleSpline[i, 10] * basisP[[i]][, 2] + resultsLeleSpline[i, 11] * basisP[[i]][, 3] +
    resultsLeleSpline[i, 12] * basisP[[i]][, 4] + resultsLeleSpline[i, 13] * basisP[[i]][, 5] + resultsLeleSpline[i, 14] * basisP[[i]][, 6]) / (1 + exp(resultsLeleSpline[i, 8] + resultsLeleSpline[i, 9] * basisP[[i]][, 1] + resultsLeleSpline[i, 10] * basisP[[i]][, 2] + resultsLeleSpline[i, 11] * basisP[[i]][, 3] +
    resultsLeleSpline[i, 12] * basisP[[i]][, 4] + resultsLeleSpline[i, 13] * basisP[[i]][, 5] + resultsLeleSpline[i, 14] * basisP[[i]][, 6]))

  idx <- order(mData[[i]])
  lines(mData[[i]][idx], fitted.val[idx], col = "red")
}

mtext(
  text = "p(z)",
  side = 2,
  line = 2.5, cex = 2
)
dev.off()
