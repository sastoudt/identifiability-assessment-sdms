
#### small data ####
set.seed(11818)
numIter=25
zData=vector("list",numIter)
mData=vector("list",numIter)
lambdaData=vector("list",numIter)
pData=vector("list",numIter)
obsData=vector("list",numIter)
obsData2=vector("list",numIter)
n=100
for(i in 1:numIter){
  
  z<- runif(2*n, -2.5,8) # simulate a covariate
  m <- runif(2*n, -2.5,8)
  
  lambda = (10*(2/441*(z+2.5)^2+0.5))*-1+15
  p = (2/441*(m+2.5)^2+0.5)*-1+1
  
  
  occurred=rpois(2*n,lambda)
  y=mapply(function(x,y){rbinom(1,x,y)},occurred,p)
  y2=mapply(function(x,y){rbinom(1,x,y)},occurred[1:n],p[1:n])
  y3=mapply(function(x,y){rbinom(1,x,y)},occurred,p)
  
  lambdaData[[i]]=lambda
  pData[[i]]=p
  zData[[i]]=z
  mData[[i]]=m
  obsData[[i]]=y
  obsData2[[i]]=y2
  print(i)
}

save(lambdaData, pData,  zData, mData, obsData, obsData2, numIter, n, file = "single_double_visit_abundance/raw_data/abunData_small.RData")

#### medium data ####

set.seed(11818)
numIter=25
zData=vector("list",numIter)
mData=vector("list",numIter)
lambdaData=vector("list",numIter)
pData=vector("list",numIter)
obsData=vector("list",numIter)
obsData2=vector("list",numIter)
n=1000
for(i in 1:numIter){
  
  z<- runif(2*n, -2.5,8) # simulate a covariate

  m <- runif(2*n, -2.5,8)

  lambda = (10*(2/441*(z+2.5)^2+0.5))*-1+15
  p =  (2/441*(m+2.5)^2+0.5)*-1+1
  
  occurred=rpois(2*n,lambda)
  y=mapply(function(x,y){rbinom(1,x,y)},occurred,p)
  y2=mapply(function(x,y){rbinom(1,x,y)},occurred[1:n],p[1:n])

  lambdaData[[i]]=lambda
  pData[[i]]=p
  zData[[i]]=z
  mData[[i]]=m
  obsData[[i]]=y
  obsData2[[i]]=y2
  print(i)
}

save(lambdaData, pData,  zData, mData, obsData, obsData2, numIter, n, file = "single_double_visit_abundance/raw_data/abunData_medium.RData")

#### large data ####

set.seed(11818)
numIter <- 25
zData <- vector("list", numIter)
mData <- vector("list", numIter)
lambdaData <- vector("list", numIter)
pData <- vector("list", numIter)
obsData <- vector("list", numIter)
obsData2 <- vector("list", numIter)
n <- 100000
for (i in 1:numIter) {
  z <- runif(2 * n, -2.5, 8) # simulate a covariate

  m <- runif(2 * n, -2.5, 8)

  lambda <- (10 * (2 / 441 * (z + 2.5)^2 + 0.5)) * -1 + 15
  p <- (2 / 441 * (m + 2.5)^2 + 0.5) * -1 + 1

  occurred <- rpois(2 * n, lambda)
  y <- mapply(function(x, y) {
    rbinom(1, x, y)
  }, occurred, p)
  y2 <- mapply(function(x, y) {
    rbinom(1, x, y)
  }, occurred[1:n], p[1:n])

  lambdaData[[i]] <- lambda
  pData[[i]] <- p
  zData[[i]] <- z
  mData[[i]] <- m
  obsData[[i]] <- y
  obsData2[[i]] <- y2
  print(i)
}

save(lambdaData, pData, zData, mData, obsData, obsData2, numIter, n, file = "single_double_visit_abundance/raw_data/abunData_large.RData")
