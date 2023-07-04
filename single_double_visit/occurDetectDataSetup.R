#### get small data ####
set.seed(11818)
numIter <- 25 ## number of data sets
n <- 100 ## size of each data set
zData <- vector("list", numIter)
mData <- vector("list", numIter)
psiData <- vector("list", numIter)
pData <- vector("list", numIter)
obsData <- vector("list", numIter)
obsData2 <- vector("list", numIter)
obsDataM <- vector("list", numIter)

for (i in 1:numIter) {
  z <- runif(2 * n, -2.5, 8) # simulate a covariate
  m <- runif(2 * n, -2.5, 8)
  psi <- 1 / 14 * z + 5 / 28
  p <- 1 / 21 * m + 5 / 42

  alpha <- 1 ## if you want to scale at some point

  occurred <- rbinom(2 * n, 1, alpha * psi)
  detected <- rbinom(2 * n, 1, pmin(1 / alpha * p, 1))
  detected2 <- rbinom(n, 1, pmin(1 / alpha * p[1:n], 1))


  y <- occurred[1:n] * detected[1:n]
  y2 <- occurred[1:n] * detected2
  yM <- occurred * detected


  psiData[[i]] <- alpha * psi
  pData[[i]] <- pmin(1 / alpha * p, 1)
  zData[[i]] <- z
  mData[[i]] <- m ## detection covariate same per visit instead of something that could change
  obsData[[i]] <- y
  obsData2[[i]] <- y2
  obsDataM[[i]] <- yM
}

save(psiData, pData, zData, mData, obsData, obsData2, obsDataM, numIter, n, file = "single_double_visit/raw_data/svdvOccurData_smallNoSaturateLine.RData")

#### get medium data ####
numIter <- 25
set.seed(11818)
zData <- vector("list", numIter)
mData <- vector("list", numIter)
psiData <- vector("list", numIter)
pData <- vector("list", numIter)
obsData <- vector("list", numIter)
obsData2 <- vector("list", numIter)
obsData3 <- vector("list", numIter)
obsDataM <- vector("list", numIter)
n <- 1000
for (i in 1:numIter) {
  z <- runif(2 * n, -2.5, 8) # simulate a covariate
  m <- runif(2 * n, -2.5, 8)
  psi <- 1 / 14 * z + 5 / 28
  p <- 1 / 21 * m + 5 / 42

  alpha <- 1

  occurred <- rbinom(2 * n, 1, alpha * psi)
  detected <- rbinom(2 * n, 1, pmin(1 / alpha * p, 1))
  detected2 <- rbinom(n, 1, pmin(1 / alpha * p[1:n], 1))

  y <- occurred[1:n] * detected[1:n]
  y2 <- occurred[1:n] * detected2
  yM <- occurred * detected


  psiData[[i]] <- alpha * psi
  pData[[i]] <- pmin(1 / alpha * p, 1)
  zData[[i]] <- z
  mData[[i]] <- m
  obsData[[i]] <- y
  obsData2[[i]] <- y2
  obsDataM[[i]] <- yM
}

save(psiData, pData, zData, mData, obsData, obsData2, obsDataM, numIter, n, file = "single_double_visit/raw_data/svdvOccurData_mediumNoSaturateLine.RData")

##### get large data ####
n <- 100000
set.seed(11818)
numIter <- 25
zData <- vector("list", numIter)
mData <- vector("list", numIter)
psiData <- vector("list", numIter)
pData <- vector("list", numIter)
obsData <- vector("list", numIter)
obsData2 <- vector("list", numIter)
obsData3 <- vector("list", numIter)
obsDataM <- vector("list", numIter)
for (i in 1:numIter) {
  z <- runif(2 * n, -2.5, 8) # simulate a covariate
  m <- runif(2 * n, -2.5, 8)
  psi <- 1 / 14 * z + 5 / 28
  p <- 1 / 21 * m + 5 / 42

  alpha <- 1

  occurred <- rbinom(2 * n, 1, alpha * psi)
  detected <- rbinom(2 * n, 1, pmin(1 / alpha * p, 1))
  detected2 <- rbinom(n, 1, pmin(1 / alpha * p[1:n], 1))

  y <- occurred[1:n] * detected[1:n]
  y2 <- occurred[1:n] * detected2
  yM <- occurred * detected


  psiData[[i]] <- alpha * psi
  pData[[i]] <- pmin(1 / alpha * p, 1)
  zData[[i]] <- z
  mData[[i]] <- m
  obsData[[i]] <- y
  obsData2[[i]] <- y2
  obsDataM[[i]] <- yM
}

save(psiData, pData, zData, mData, obsData, obsData2, obsDataM, numIter, n, file = "single_double_visit/raw_data/svdvOccurData_largeNoSaturateLine.RData")
