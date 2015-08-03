setwd("~/Desktop/Research/")

require(Amelia)
require(numDeriv)

#linear model comparison
data = read.csv("2008.csv",header = T)

lmcomparison <- function(data)
{
  n = nrow(data)
  lmfit.time <- system.time(lmfit <- lm(ArrDelay~DepDelay+Distance , data=data))
  #Prepare subset data
  idx = sample(n, 100000, replace = FALSE)
  y = data$ArrDelay[idx]
  x = as.data.frame(cbind(data$DepDelay[idx],data$Distance[idx]))
  names(x) = c("DepDelay", "Distance")
  #amelia2
  amdata = cbind(x,y)
  amfit.time<- system.time({amX <- amelia(x = amdata,p2s = 0,m = 5)$imputations
                          amMat <- Reduce("+", amX) / length(amX)
                          amfit <- lm(y ~ . , data = amMat)})
  x = cbind(1,x)
  mvfit.time <- system.time(mvfit <- lmmv(x = x, y = y))
  rest <- list(lmfit.coeffecoents = lmfit$coefficients,
               amfit.coeffecients = amfit$coefficients,
               mvfit.coeffecients = mvfit$beta,
               lmfit.time = lmfit.time,
               amfit.time = amfit.time,
               mvfit.time = mvfit.time)
  rest
}
#######################
#> lmcomparison(data)
#$lmfit.coeffecoents
#(Intercept)     DepDelay     Distance 
#-1.061369464  1.019154260 -0.001213193 
#
#$amfit.coeffecients
#(Intercept)     DepDelay     Distance 
#-1.164011527  1.021171995 -0.001028826 
#
#$mvfit.coeffecients
#[1] -1.109610463  1.011967440 -0.001056084
#
#$lmfit.time # Alldata
#user  system elapsed 
#9.363   1.799  15.158 
#
#$amfit.time
#user  system elapsed 
#21.673   1.319  26.937 
#
#$mvfit.time
#user  system elapsed 
#0.544   0.288   1.009 


#PCA comparison
#prepare data
idx = sample(n,100000,FALSE)
pcdata = data[idx ,c (12:16, 19:21)]

#Amelia does not work here
tryCatch(amX <- amelia(x = pcdata, p2s = 0, m = 5),finally = print("Amelia can not work here!"))


#normal PCA

pcacomparison <- function(pcdata)
{
  pca.time <- system.time({cc <- na.omit(pcdata)
                         ccpca <- prcomp(x = cc)$sdev[1:3]})
  mvpca.time1 <- system.time(corPcaMV <- PCAmv(data = pcdata,retX = FALSE,Cor = TRUE)$eigenvalues[1:3])
  mvpca.time2 <- system.time(covPcaMV <- PCAmv(data = pcdata,retX = FALSE,Cor = FALSE)$eigenvalues[1:3])
  res <- list(pca.time = pca.time,
              CorrelationMat.time = mvpca.time1,
              CovarianceMat.time = mvpca.time2,
              ccpca = ccpca,
              corPcaMV =sqrt(corPcaMV),
              covPcaMV =sqrt(covPcaMV))
  res
}
pcacomparison(pcdata)  
#ccpca
#mvpca

#pca.time



#######################

#log-linear comparison


