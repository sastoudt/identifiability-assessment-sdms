
#### small data ####
set.seed(11818)
numIter <- 25
covData <- vector("list", numIter)
psiData <- vector("list", numIter)
obsData <- vector("list", numIter)
n <- 100

for (i in 1:numIter) {
  N <- rpois(1, n) ## number of sites visited
  z <- runif(2 * N, -2.5, 8) # simulate a covariate
  lpsi <- -1 + 1 * z
  psi <- exp(lpsi) / (1 + exp(lpsi)) ### occurrence
  y <- rbinom(2 * N, 1, 0.5 * psi) #### observed data
  psiData[[i]] <- psi
  covData[[i]] <- z
  obsData[[i]] <- y
}

save(psiData, covData, obsData, numIter, n, file = "presence_absence_only/raw_data/popaData_small.RData")


#### medium data ####
numIter <- 25
covDataM <- vector("list", numIter)
psiDataM <- vector("list", numIter)
obsDataM <- vector("list", numIter)

set.seed(11818)
n <- 1000

for (i in 1:numIter) {
  N <- rpois(1, n)
  z <- runif(2 * N, -2.5, 8) # simulate a covariate
  lpsi <- -1 + 1 * z

  psi <- exp(lpsi) / (1 + exp(lpsi)) ## inv logit
  y <- rbinom(2 * N, 1, 0.5 * psi)

  covDataM[[i]] <- z
  psiDataM[[i]] <- psi
  obsDataM[[i]] <- y
}
save(psiDataM, covDataM, obsDataM, numIter, n, file = "presence_absence_only/raw_data/popaData_medium.RData")


#### large data ####
set.seed(11818)
numIter <- 25
covDataM <- vector("list", numIter)
psiDataM <- vector("list", numIter)
obsDataM <- vector("list", numIter)
n <- 100000

for (i in 1:numIter) {
  N <- rpois(1, n) ## number of sites visited
  z <- runif(2 * N, -2.5, 8) # simulate a covariate
  lpsi <- -1 + 1 * z
  psi <- exp(lpsi) / (1 + exp(lpsi)) ### occurrence
  y <- rbinom(2 * N, 1, 0.5 * psi) #### observed data
  psiDataM[[i]] <- psi
  covDataM[[i]] <- z
  obsDataM[[i]] <- y
}

save(psiDataM, covDataM, obsDataM, numIter, n, file = "presence_absence_only/raw_data/popaData_large.RData")
