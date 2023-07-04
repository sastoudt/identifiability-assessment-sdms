load(file = "single_double_visit_abundance/raw_data/abunData_small.RData")

#### double visit ####

numParam <- 4
resultsLogReg <- matrix(NA, nrow = numIter, ncol = numParam)

for (i in 1:numIter) {
  dataMat <- matrix(NA, n, 2)
  dataMat[, 1] <- obsData[[i]][1:n]
  dataMat[, 2] <- obsData2[[i]]
  visitMat <- matrix(NA, n, 2)
  visitMat[, 1] <- visitMat[, 2] <- mData[[i]][1:n]
  umf <- unmarkedFramePCount(y = dataMat, siteCovs = data.frame(x = zData[[i]][1:n]), obsCovs = list(visit = visitMat))
  test <- pcount(~visit ~ x, umf, K = 100)
  resultsLogReg[i, ] <- coef(test)

  print(i)
}

save(resultsLogReg, file = "single_double_visit_abundance/results/dvAbundanceSmall_noSpline.RData")

#### Figure S1e and S2e ####

load("single_double_visit_abundance/raw_data/abunData_small.RData")
load("single_double_visit_abundance/results/dvAbundanceSmall_noSpline.RData")


thin=1:n
jpeg("paper_figures/figS1e.jpeg",width=824, height=634,units='px')
idx = order(zData[[1]][1:n])
plot(zData[[1]][idx],lambdaData[[1]][idx],lwd=.5,xlim=c(-3,8),xlab="x",ylab="",cex.axis=2,cex.lab=2,ylim=c(0,25),type="l")
for(i in 2:25){
  idx = order(zData[[i]][1:n])
  lines(zData[[i]][idx],lambdaData[[i]][idx],lwd=.5)
}


for(i in 1:25){
  orderIdx=order(zData[[i]][1:n])
  lines(zData[[i]][orderIdx],exp(resultsLogReg[i,1]+resultsLogReg[i,2]*zData[[i]][orderIdx]),col="red")
}

i=1
idx = order(zData[[i]][1:n])
lines(zData[[i]][idx],lambdaData[[i]][idx],lwd=2)
mtext(text =expression(paste(lambda,"(x)")),
      side = 2,
      line = 2,cex=2)

dev.off()

jpeg("paper_figures/figS2e.jpeg",width=824, height=634,units='px')
idx = order(mData[[1]][1:n])
plot(mData[[1]][idx],pData[[1]][idx],ylim=c(0,1),lwd=.5,xlim=c(-3,8),xlab="z",ylab="",cex.axis=2,cex.lab=2,type="l")
for(i in 2:25){
  idx = order(mData[[i]][1:n])
  lines(mData[[i]][idx],pData[[i]][idx],lwd=.5)
}

for(i in 1:25){
  orderIdx=order(mData[[i]][1:n])
  
  lines(mData[[i]][orderIdx],exp(resultsLogReg[i,3]+resultsLogReg[i,4]*mData[[i]][orderIdx])/(1+exp(resultsLogReg[i,3]+resultsLogReg[i,4]*mData[[i]][orderIdx])),col="red")
}

mtext(text ="p(z)",
      side = 2, 
      line = 2.5,cex=2)
dev.off()


#### single visit ####
numParam <- 4
resultsLele <- matrix(NA, nrow = numIter, ncol = numParam)
for (i in 1:numIter) {
  mod <- svabu(obsData[[i]] ~ zData[[i]] | mData[[i]], dist = "P", zeroinfl = F)
  resultsLele[i, ] <- coef(mod)

  print(i)
}

save(resultsLele, file = "single_double_visit_abundance/results/svAbundanceSmall_noSpline.RData")

#### Figure S1a and S2a ####

load("single_double_visit_abundance/raw_data/abunData_small.RData")
load("single_double_visit_abundance/results/svAbundanceSmall_noSpline.RData")

jpeg("paper_figures/figS1a.jpeg",width=824, height=634,units='px')
idx = order(zData[[1]])
plot(zData[[1]][idx],lambdaData[[1]][idx],lwd=.5,xlim=c(-3,8),ylim=c(0,25),xlab="x",ylab="",cex.axis=2,cex.lab=2,type="l")
for(i in 2:25){
  idx = order(zData[[i]])
  lines(zData[[i]][idx],lambdaData[[i]][idx],lwd=.5)
}


for(i in 1:25){
  idx = order(zData[[i]])
  lines(zData[[i]][idx],exp(resultsLele[i,1]+resultsLele[i,2]*zData[[i]][idx]),col="red")
}


mtext(text =expression(paste(lambda,"(x)")),
      side = 2, 
      line = 2,cex=2)
dev.off()

jpeg("paper_figures/figS2a.jpeg",width=824, height=634,units='px')
idx = order(mData[[1]])
plot(mData[[1]][idx],pData[[1]][idx],ylim=c(0,1),lwd=.5,xlim=c(-3,8),xlab="z",ylab="",cex.axis=2,cex.lab=2,type="l")
for(i in 2:25){
  idx = order(mData[[i]])
  lines(mData[[i]][idx],pData[[i]][idx],lwd=.5)
}

for(i in 1:25){
  idx = order(mData[[i]])
  lines(mData[[i]][idx],exp(resultsLele[i,3]+resultsLele[i,4]*mData[[i]][idx])/(1+exp(resultsLele[i,3]+resultsLele[i,4]*mData[[i]][idx])),col="red")
}

mtext(text ="p(z)",
      side = 2, 
      line = 2.5,cex=2)
dev.off()
