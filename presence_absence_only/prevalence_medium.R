
load(file = "presence_absence_only/raw_data/popaData_medium.RData")

#### presence absence ####
resultsLogRegSpline <- matrix(NA, nrow = numIter, ncol = 11)
basisMfull <- vector("list", numIter)

for (i in 1:numIter) {
  zShuffle <- covDataM[[i]][1:(length(covDataM[[i]]) / 2)]
  dataShuffle <- obsDataM[[i]][1:(length(covDataM[[i]]) / 2)]

  m <- 1:7

  basisM <- bs(zShuffle, knots = -2.5 + m * 10.5 / 8)
  basisMfull[[i]] <- basisM


  mod <- glm(dataShuffle ~ basisM[, 1] + basisM[, 2] + basisM[, 3] + basisM[, 4] + basisM[, 5] + basisM[, 6] + basisM[, 7] + basisM[, 8] + basisM[, 9] + basisM[, 10], family = "binomial")


  print(i)
  resultsLogRegSpline[i, ] <- coefficients(mod)
}

save(resultsLogRegSpline, basisMfull, file = "presence_absence_only/results/paSplineResults_mediumMoreKnots.RData")

#### Figure 2g ####
load("presence_absence_only/results/paSplineResults_mediumMoreKnots.RData")

thin <- seq(1, n, by = 1)
jpeg("paper_figures/fig2g.jpeg", width = 824, height = 634, units = "px")
plot(sort(covDataM[[1]][thin]), sort(0.5 * psiDataM[[1]][thin]), ylim = c(0, 1), type = "l", lwd = .5, xlab = "x", ylab = "", cex.axis = 2, cex.lab = 2)
for (i in 2:25) {
  lines(sort(covDataM[[i]][thin]), sort(0.5 * psiDataM[[i]][thin]), lwd = .5)
}
for (i in 1:25) {
  fitted.val <- exp(resultsLogRegSpline[i, 1] + resultsLogRegSpline[i, 2] * basisMfull[[i]][, 1] + resultsLogRegSpline[i, 3] * basisMfull[[i]][, 2] + resultsLogRegSpline[i, 4] * basisMfull[[i]][, 3] +
    resultsLogRegSpline[i, 5] * basisMfull[[i]][, 4] + resultsLogRegSpline[i, 6] * basisMfull[[i]][, 5] + resultsLogRegSpline[i, 7] * basisMfull[[i]][, 6]

    + resultsLogRegSpline[i, 8] * basisMfull[[i]][, 7]
    + resultsLogRegSpline[i, 9] * basisMfull[[i]][, 8]
    + resultsLogRegSpline[i, 10] * basisMfull[[i]][, 9]
    + resultsLogRegSpline[i, 11] * basisMfull[[i]][, 10]) / (1 + exp(resultsLogRegSpline[i, 1] + resultsLogRegSpline[i, 2] * basisMfull[[i]][, 1] + resultsLogRegSpline[i, 3] * basisMfull[[i]][, 2] + resultsLogRegSpline[i, 4] * basisMfull[[i]][, 3] +
    resultsLogRegSpline[i, 5] * basisMfull[[i]][, 4] + resultsLogRegSpline[i, 6] * basisMfull[[i]][, 5] + resultsLogRegSpline[i, 7] * basisMfull[[i]][, 6] + resultsLogRegSpline[i, 8] * basisMfull[[i]][, 7]
    + resultsLogRegSpline[i, 9] * basisMfull[[i]][, 8]
    + resultsLogRegSpline[i, 10] * basisMfull[[i]][, 9]
    + resultsLogRegSpline[i, 11] * basisMfull[[i]][, 10]))
  idx <- order(covDataM[[i]][1:n])
  points((covDataM[[i]][idx])[thin], (fitted.val[idx])[thin], col = "red", cex = .5)
}

lines(sort(covDataM[[1]][thin]), sort(0.5 * psiDataM[[1]][thin]))
mtext(
  text = expression(paste(psi, "(x)")),
  side = 2,
  line = 2.5, cex = 2
)
dev.off()

#### Figure 3b ####
jpeg("paper_figures/fig3b.jpeg", width = 824, height = 634, units = "px")

plot(seq(-2.5, 8, length.out = 100), rep(1 / 10.5, 100), ylim = c(0, 0.2), type = "l", ylab = "", xlab = "x", cex.axis = 2, cex.lab = 2, main = "Presence-Absence \n Observable Distribution", cex.main = 2)
lines(c(-2.5, -2.5), c(0, 1 / 10.5), lty = 2)
lines(c(8, 8), c(0, 1 / 10.5), lty = 2)
idx <- order(covDataM[[1]])
lines(covDataM[[1]][idx], ((psiDataM[[1]] * 1 / 10.5) / mean(psiDataM[[1]]))[idx])
for (i in 2:25) {
  idx <- order(covDataM[[i]])
  lines(covDataM[[i]][idx], ((psiDataM[[i]] * 1 / 10.5) / mean(psiDataM[[i]]))[idx])
}
lines(c(-2.5, -2.5), c(0, min((psiDataM[[1]] * 1 / 10.5) / mean(psiDataM[[1]]))), lty = 2)
lines(c(8, 8), c(0, max((psiDataM[[1]] * 1 / 10.5) / mean(psiDataM[[1]]))), lty = 2)
for (i in 1:25) {
  fitted.val <- exp(resultsLogRegSpline[i, 1] + resultsLogRegSpline[i, 2] * basisMfull[[i]][, 1] + resultsLogRegSpline[i, 3] * basisMfull[[i]][, 2] + resultsLogRegSpline[i, 4] * basisMfull[[i]][, 3] +
    resultsLogRegSpline[i, 5] * basisMfull[[i]][, 4] + resultsLogRegSpline[i, 6] * basisMfull[[i]][, 5] + resultsLogRegSpline[i, 7] * basisMfull[[i]][, 6]
    + resultsLogRegSpline[i, 8] * basisMfull[[i]][, 7]
    + resultsLogRegSpline[i, 9] * basisMfull[[i]][, 8]
    + resultsLogRegSpline[i, 10] * basisMfull[[i]][, 9]
    + resultsLogRegSpline[i, 11] * basisMfull[[i]][, 10]) / (1 + exp(resultsLogRegSpline[i, 1] + resultsLogRegSpline[i, 2] * basisMfull[[i]][, 1] + resultsLogRegSpline[i, 3] * basisMfull[[i]][, 2] + resultsLogRegSpline[i, 4] * basisMfull[[i]][, 3] +
    resultsLogRegSpline[i, 5] * basisMfull[[i]][, 4] + resultsLogRegSpline[i, 6] * basisMfull[[i]][, 5] + resultsLogRegSpline[i, 7] * basisMfull[[i]][, 6]
    + resultsLogRegSpline[i, 8] * basisMfull[[i]][, 7]
    + resultsLogRegSpline[i, 9] * basisMfull[[i]][, 8]
    + resultsLogRegSpline[i, 10] * basisMfull[[i]][, 9]
    + resultsLogRegSpline[i, 11] * basisMfull[[i]][, 10]))
  idx <- order(covDataM[[i]][1:n])
  lines(covDataM[[i]][idx], ((fitted.val * (1 / 10.5)) / mean(fitted.val))[idx], col = "red")
}
idx <- order(covDataM[[1]][1:n])
lines(covDataM[[1]][idx], ((psiDataM[[1]] * 1 / 10.5) / mean(psiDataM[[1]]))[idx])
mtext(
  text = "P(x | y = 1)",
  side = 2, # side 2 = left
  line = 2.5, cex = 2
)
dev.off()



#### presence only ####


resultsRoyleSpline <- matrix(NA, nrow = numIter, ncol = 11)
basisMfull <- basisMDfull <- vector("list", numIter)
for (i in 1:numIter) {
  zShuffle <- covDataM[[i]]
  dataShuffle <- covDataM[[i]][obsDataM[[i]] == 1]

  m <- 1:7


  basisM <- bs(zShuffle, knots = -2.5 + m * 10.5 / 8)
  basisMD <- bs(dataShuffle, knots = -2.5 + m * 10.5 / 8)
  basisMfull[[i]] <- basisM
  basisMDfull[[i]] <- basisMD
  lik <- function(parm) {
    beta0 <- parm[1]
    beta1 <- parm[2]
    beta2 <- parm[3]
    beta3 <- parm[4]
    beta4 <- parm[5]
    beta5 <- parm[6]
    beta6 <- parm[7]
    beta7 <- parm[8]
    beta8 <- parm[9]
    beta9 <- parm[10]
    beta10 <- parm[11]

    gridpsiIn <- beta0 + beta1 * basisM[, 1] + beta2 * basisM[, 2] + beta3 * basisM[, 3] +
      beta4 * basisM[, 4] + beta5 * basisM[, 5] + beta6 * basisM[, 6] + beta7 * basisM[, 7] + beta8 * basisM[, 8] + beta9 * basisM[, 9] + beta10 * basisM[, 10]

    gridpsi <-
      exp(gridpsiIn) / (1 + exp(gridpsiIn))

    datapsiIn <- beta0 + beta1 * basisMD[, 1] + beta2 * basisMD[, 2] + beta3 * basisMD[, 3] +
      beta4 * basisMD[, 4] + beta5 * basisMD[, 5] + beta6 * basisMD[, 6] + beta7 * basisMD[, 7] + beta8 * basisMD[, 8] + beta9 * basisMD[, 9] + beta10 * basisMD[, 10]

    datapsi <-
      exp(datapsiIn) / (1 + exp(datapsiIn))

    -1 * sum(log(datapsi / (sum(gridpsi))))
  }

  out <- tryCatch(optim(rep(0.1, 11), lik, hessian = T), error = function(i) {
    return(rep(NA, 7))
  })
  print(i)
  resultsRoyleSpline[i, ] <- out$par
}

save(resultsRoyleSpline, basisMfull, basisMDfull, file = "presence_absence_only/results/poSplineResults_mediumMoreKnots.RData")

#### Figure 2c ####

load("presence_absence_only/results/poSplineResults_mediumMoreKnots.RData")

thin <- seq(1, 2 * n, by = 1)
jpeg("paper_figures/fig2c.jpeg", width = 824, height = 634, units = "px")
plot(sort(covDataM[[1]][thin]), sort(0.5 * psiDataM[[1]][thin]), ylim = c(0, 1), type = "l", lwd = .5, xlab = "x", ylab = "", cex.axis = 2, cex.lab = 2)
for (i in 2:25) {
  lines(sort(covDataM[[i]][thin]), sort(0.5 * psiDataM[[i]][thin]), lwd = .5)
}
for (i in 1:25) {
  thin <- seq(1, length(covDataM[[i]]), by = 100)
  fitted.val <- exp(resultsRoyleSpline[i, 1] + resultsRoyleSpline[i, 2] * basisMfull[[i]][, 1] + resultsRoyleSpline[i, 3] * basisMfull[[i]][, 2] + resultsRoyleSpline[i, 4] * basisMfull[[i]][, 3] +
    resultsRoyleSpline[i, 5] * basisMfull[[i]][, 4] + resultsRoyleSpline[i, 6] * basisMfull[[i]][, 5] + resultsRoyleSpline[i, 7] * basisMfull[[i]][, 6]

    + resultsRoyleSpline[i, 8] * basisMfull[[i]][, 7]
    + resultsRoyleSpline[i, 9] * basisMfull[[i]][, 8]
    + resultsRoyleSpline[i, 10] * basisMfull[[i]][, 9]
    + resultsRoyleSpline[i, 11] * basisMfull[[i]][, 10]) / (1 + exp(resultsRoyleSpline[i, 1] + resultsRoyleSpline[i, 2] * basisMfull[[i]][, 1] + resultsRoyleSpline[i, 3] * basisMfull[[i]][, 2] + resultsRoyleSpline[i, 4] * basisMfull[[i]][, 3] +
    resultsRoyleSpline[i, 5] * basisMfull[[i]][, 4] + resultsRoyleSpline[i, 6] * basisMfull[[i]][, 5] + resultsRoyleSpline[i, 7] * basisMfull[[i]][, 6] + resultsRoyleSpline[i, 8] * basisMfull[[i]][, 7]
    + resultsRoyleSpline[i, 9] * basisMfull[[i]][, 8]
    + resultsRoyleSpline[i, 10] * basisMfull[[i]][, 9]
    + resultsRoyleSpline[i, 11] * basisMfull[[i]][, 10]))
  idx <- order(covDataM[[i]])
  lines((covDataM[[i]][idx])[thin], (fitted.val[idx])[thin], col = "red", cex = .5)
}

mtext(
  text = expression(paste(psi, "(x)")),
  side = 2, # side 2 = left
  line = 2.5, cex = 2
)
dev.off()


#### Figure 3a ####
jpeg("paper_figures/fig3a.jpeg", width = 824, height = 634, units = "px")

plot(seq(-2.5, 8, length.out = 100), rep(1 / 10.5, 100), ylim = c(0, 0.2), type = "l", ylab = "", xlab = "x", cex.axis = 2, cex.lab = 2, main = "Presence-Only Royle \n Observable Distribution", cex.main = 2)
lines(c(-2.5, -2.5), c(0, 1 / 10.5), lty = 2)
lines(c(8, 8), c(0, 1 / 10.5), lty = 2)
idx <- order(covDataM[[1]])
lines(covDataM[[1]][idx], ((psiDataM[[1]] * 1 / 10.5) / mean(psiDataM[[1]]))[idx])
for (i in 2:25) {
  idx <- order(covDataM[[i]])
  lines(covDataM[[i]][idx], ((psiDataM[[i]] * 1 / 10.5) / mean(psiDataM[[i]]))[idx])
}
lines(c(-2.5, -2.5), c(0, min((psiDataM[[1]] * 1 / 10.5) / mean(psiDataM[[1]]))), lty = 2)
lines(c(8, 8), c(0, max((psiDataM[[1]] * 1 / 10.5) / mean(psiDataM[[1]]))), lty = 2)
for (i in 1:25) {
  fitted.val <- exp(resultsRoyleSpline[i, 1] + resultsRoyleSpline[i, 2] * basisMfull[[i]][, 1] + resultsRoyleSpline[i, 3] * basisMfull[[i]][, 2] + resultsRoyleSpline[i, 4] * basisMfull[[i]][, 3] +
    resultsRoyleSpline[i, 5] * basisMfull[[i]][, 4] + resultsRoyleSpline[i, 6] * basisMfull[[i]][, 5] + resultsRoyleSpline[i, 7] * basisMfull[[i]][, 6]
    + resultsRoyleSpline[i, 8] * basisMfull[[i]][, 7]
    + resultsRoyleSpline[i, 9] * basisMfull[[i]][, 8]
    + resultsRoyleSpline[i, 10] * basisMfull[[i]][, 9]
    + resultsRoyleSpline[i, 11] * basisMfull[[i]][, 10]) / (1 + exp(resultsRoyleSpline[i, 1] + resultsRoyleSpline[i, 2] * basisMfull[[i]][, 1] + resultsRoyleSpline[i, 3] * basisMfull[[i]][, 2] + resultsRoyleSpline[i, 4] * basisMfull[[i]][, 3] +
    resultsRoyleSpline[i, 5] * basisMfull[[i]][, 4] + resultsRoyleSpline[i, 6] * basisMfull[[i]][, 5] + resultsRoyleSpline[i, 7] * basisMfull[[i]][, 6]

    + resultsRoyleSpline[i, 8] * basisMfull[[i]][, 7]
    + resultsRoyleSpline[i, 9] * basisMfull[[i]][, 8]
    + resultsRoyleSpline[i, 10] * basisMfull[[i]][, 9]
    + resultsRoyleSpline[i, 11] * basisMfull[[i]][, 10]))
  idx <- order(covDataM[[i]])
  lines((covDataM[[i]][idx])[thin], (((fitted.val * (1 / 10.5)) / mean(fitted.val))[idx])[thin], col = "red")
}
idx <- order(covDataM[[1]])
lines(covDataM[[1]][idx], ((psiDataM[[1]] * 1 / 10.5) / mean(psiDataM[[1]]))[idx])
mtext(
  text = "P(x | y = 1)",
  side = 2,
  line = 2.5, cex = 2
)
dev.off()
