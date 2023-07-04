load(file = "single_double_visit_abundance/raw_data/abunData_small.RData")


#### double visit ####
numParam <- 14
resultsLogRegSpline <- matrix(NA, nrow = numIter, ncol = numParam)

basisPsi <- basisP <- vector("list", numIter)
for (i in 1:numIter) {
  basisPsi[[i]] <- basisPsiUse <- bs(zData[[i]][1:n], knots = c(-2.5 + 10.5 / 4, -2.5 + 10.5 / 4 + 10.5 / 4, -2.5 + 10.5 / 4 + 10.5 / 4 + 10.5 / 4))
  basisP[[i]] <- basisPUse <- bs(mData[[i]][1:n], knots = c(-2.5 + 10.5 / 4, -2.5 + 10.5 / 4 + 10.5 / 4, -2.5 + 10.5 / 4 + 10.5 / 4 + 10.5 / 4))

  dataMat <- matrix(NA, n, 2)
  dataMat[, 1] <- obsData[[i]][1:n]
  dataMat[, 2] <- obsData2[[i]]
  visitMat1 <- visitMat2 <- visitMat3 <- visitMat4 <- visitMat5 <- visitMat6 <- matrix(NA, n, 2)

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

save(resultsLogRegSpline, basisPsi, basisP, file = "single_double_visit_abundance/results/dvAbundanceSmall.RData")

#### Figures S1f and S2f ####

load(file = "single_double_visit_abundance/raw_data/abunData_small.RData")
load(file = "single_double_visit_abundance/results/dvAbundanceSmall.RData")

thin <- 1:length(zData[[1]][1:n])
jpeg("paper_figures/figS1f.jpeg", width = 824, height = 634, units = "px")
idx <- order(zData[[1]][1:n])
plot(zData[[1]][idx], lambdaData[[1]][idx], lwd = .5, xlim = c(-3, 8), xlab = "x", ylab = "", cex.axis = 2, cex.lab = 2, ylim = c(0, 25), type = "l")
for (i in 2:25) {
  idx <- order(zData[[i]][1:n])
  lines(zData[[i]][idx], lambdaData[[i]][idx], lwd = .5)
}


for (i in 1:25) {
  fitted.val <- exp(resultsLogRegSpline[i, 1] + resultsLogRegSpline[i, 2] * basisPsi[[i]][, 1] + resultsLogRegSpline[i, 3] * basisPsi[[i]][, 2] + resultsLogRegSpline[i, 4] * basisPsi[[i]][, 3] +
    resultsLogRegSpline[i, 5] * basisPsi[[i]][, 4] + resultsLogRegSpline[i, 6] * basisPsi[[i]][, 5] + resultsLogRegSpline[i, 7] * basisPsi[[i]][, 6])
  idx <- order(zData[[i]][1:n])
  lines(zData[[i]][idx], fitted.val[idx], col = "red")
}

i <- 1
idx <- order(zData[[i]][1:n])
lines(zData[[i]][idx], lambdaData[[i]][idx], lwd = 2)

mtext(
  text = expression(paste(lambda, "(x)")),
  side = 2,
  line = 2, cex = 2
)
dev.off()

jpeg("paper_figures/figS2f.jpeg", width = 824, height = 634, units = "px")
idx <- order(mData[[1]][1:n])
plot(mData[[1]][idx], pData[[1]][idx], ylim = c(0, 1), lwd = .5, xlim = c(-3, 8), xlab = "z", ylab = "", cex.axis = 2, cex.lab = 2, type = "l")
for (i in 2:25) {
  idx <- order(mData[[i]][1:n])
  lines(mData[[i]][idx], pData[[i]][idx], lwd = .5)
}


for (i in 1:25) {
  fitted.val <- exp(resultsLogRegSpline[i, 8] + resultsLogRegSpline[i, 9] * basisP[[i]][, 1] + resultsLogRegSpline[i, 10] * basisP[[i]][, 2] + resultsLogRegSpline[i, 11] * basisP[[i]][, 3] +
    resultsLogRegSpline[i, 12] * basisP[[i]][, 4] + resultsLogRegSpline[i, 13] * basisP[[i]][, 5] + resultsLogRegSpline[i, 14] * basisP[[i]][, 6]) / (1 + exp(resultsLogRegSpline[i, 8] + resultsLogRegSpline[i, 9] * basisP[[i]][, 1] + resultsLogRegSpline[i, 10] * basisP[[i]][, 2] + resultsLogRegSpline[i, 11] * basisP[[i]][, 3] +
    resultsLogRegSpline[i, 12] * basisP[[i]][, 4] + resultsLogRegSpline[i, 13] * basisP[[i]][, 5] + resultsLogRegSpline[i, 14] * basisP[[i]][, 6]))

  idx <- order(mData[[i]][1:n])
  lines(mData[[i]][idx], fitted.val[idx], col = "red")
}

mtext(
  text = "p(z)",
  side = 2,
  line = 2.5, cex = 2
)
dev.off()


#### single visit ####

resultsLeleSpline <- matrix(NA, nrow = numIter, ncol = 14)


for (i in 1:numIter) {
  m <- 1:7

  basisPsi <- bs(zData[[i]], knots = -2.5 + m * 10.5 / 8)
  basisP <- bs(mData[[i]], knots = -2.5 + m * 10.5 / 8)

  mod <- svabu(obsData[[i]] ~ basisPsi[, 1] + basisPsi[, 2] + basisPsi[, 3] + basisPsi[, 4] + basisPsi[, 5] + basisPsi[, 6] | basisP[, 1] + basisP[, 2] + basisP[, 3] + basisP[, 4] + basisP[, 5] + basisP[, 6], dist = "P", zeroinfl = F)
  resultsLeleSpline[i, ] <- coef(mod)


  print(i)
}

basisPsi <- basisP <- vector("list", numIter)

m <- 1:7

for (i in 1:numIter) {
  basisPsi[[i]] <- bs(zData[[i]], knots = -2.5 + m * 10.5 / 8)
  basisP[[i]] <- bs(mData[[i]], knots = -2.5 + m * 10.5 / 8)
}


save(resultsLeleSpline, basisPsi, basisP, file = "single_double_visit_abundance/results/svAbundanceSmall.RData")

#### Figures S1b and S2b ####

load(file = "single_double_visit_abundance/raw_data/abunData_small.RData")
load(file = "single_double_visit_abundance/results/svAbundanceSmall.RData")


jpeg("paper_figures/figS1b.jpeg", width = 824, height = 634, units = "px")
idx <- order(zData[[1]])
plot(zData[[1]][idx], lambdaData[[1]][idx], lwd = .5, xlim = c(-3, 8), xlab = "x", ylab = "", cex.axis = 2, cex.lab = 2, ylim = c(0, 25), type = "l")
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

jpeg("paper_figures/figS2b.jpeg", width = 824, height = 634, units = "px")
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
