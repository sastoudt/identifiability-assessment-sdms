
load("presence_absence_only/raw_data/popaData_small.RData")

#### presence absence ####

resultsLogReg <- matrix(NA, nrow = numIter, ncol = 2)

for (i in 1:numIter) {
  zShuffle <- covData[[i]][1:(length(covData[[i]]) / 2)]
  dataShuffle <- obsData[[i]][1:(length(covData[[i]]) / 2)]

  mod <- glm(dataShuffle ~ zShuffle, family = "binomial") ## presence absence

  resultsLogReg[i, ] <- coefficients(mod)
}

save(resultsLogReg, file = "presence_absence_only/results/paResults_small.RData")

#### Figure 2e ####

load("presence_absence_only/results/paResults_small.RData")

thin <- 1:n

jpeg("paper_figures/fig2e.jpeg", width = 824, height = 634, units = "px")
plot(sort(covData[[1]][thin]), sort(0.5 * psiData[[1]][thin]), ylim = c(0, 1), type = "l", lwd = .5, xlim = c(-3, 8), xlab = "x", ylab = "", cex.axis = 2, cex.lab = 2)
for (i in 2:25) {
  lines(sort(covData[[i]][thin]), sort(0.5 * psiData[[i]][thin]), lwd = .5)
}
for (i in 1:25) {
  idx <- order(covData[[i]])
  lines(covData[[i]][idx], exp(resultsLogReg[i, 1] + resultsLogReg[i, 2] * covData[[i]][idx]) / (1 + exp(resultsLogReg[i, 1] + resultsLogReg[i, 2] * covData[[i]][idx])), col = "red", cex = .5)
}

mtext(
  text = expression(paste(psi, "(x)")),
  side = 2,
  line = 2.5, cex = 2
)
dev.off()
#### presence only ####

resultsRoyle <- matrix(NA, nrow = numIter, ncol = 2)

for (i in 1:numIter) {
  zShuffle <- covData[[i]]
  dataShuffle <- covData[[i]][obsData[[i]] == 1]
  lik <- function(parm) {
    beta0 <- parm[1]
    beta1 <- parm[2]
    gridpsi <-
      exp(beta0 + beta1 * zShuffle) / (1 + exp(beta0 + beta1 * zShuffle))
    datapsi <-
      exp(beta0 + beta1 * dataShuffle) / (1 + exp(beta0 +
        beta1 * dataShuffle))
    -1 * sum(log(datapsi / (sum(gridpsi))))
  }

  out <- tryCatch(optim(rep(0.5, 2), lik, hessian = T), error = function(i) {
    return(c(NA, NA))
  })

  resultsRoyle[i, ] <- out$par
}

save(resultsRoyle, file = "presence_absence_only/results/poResults_small.RData")

#### Figure 2a ####

load("presence_absence_only/results/poResults_small.RData")


idx <- order(covData[[1]])
jpeg("paper_figures/fig2a.jpeg", width = 824, height = 634, units = "px")
plot(covData[[1]][idx], 0.5 * psiData[[1]][idx], lwd = .5, xlab = "x", ylim = c(0, 1), ylab = "", cex.axis = 2, cex.lab = 2, type = "l") # ,main="Presence-Only \n Misspecified Royle's Method",cex.main=2)#,xlim=c(-3,5))
for (i in 2:25) {
  idx <- order(covData[[i]])
  lines(covData[[i]][idx], 0.5 * psiData[[i]][idx], lwd = .5)
}
for (i in 1:25) {
  idx <- order(covData[[i]])

  lines(covData[[i]][idx], exp(resultsRoyle[i, 1] + resultsRoyle[i, 2] * covData[[i]][idx]) / (1 + exp(resultsRoyle[i, 1] + resultsRoyle[i, 2] * covData[[i]][idx])), col = "red", cex = .5)
}
mtext(
  text = expression(paste(psi, "(x)")),
  side = 2,
  line = 2.5, cex = 2
)
dev.off()