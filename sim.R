require(Amelia)
require(numDeriv)

# n: number of observations
# missp : probility of missing data

makeNA <- function(n,missp) {
  nmiss <- rbinom(1,n,missp)
  sample(1:n,nmiss,replace=FALSE)
}

# nrep: number of replicate
# n: number of rows
# p: number of parameters
# NAprob: probility of NA

simlm <- function(nrep, n, p, NAprob)
{
  onerep <- function(n, p, NAprob)
  {
    np <- n*p
    x <- matrix(runif(np), ncol = p)
    y <- rowSums(x) + runif(n)
    NAidx <- makeNA(np, NAprob)
    NAidy <- makeNA(n, NAprob)
    x[NAidx] <- NA
    y[NAidy] <- NA
    #using elapsed time
    cctime <- system.time(lmfit <- lm(y ~ x)$coefficients[2])[3]
    amtime <- system.time({am <- amelia(x = x, m = 5,p2s = 0)$imputations
    amX <- Reduce("+", am)/length(am)
    amfit <- lm(y ~ amX)$coefficients[2]})[3]
    x <- cbind(1,x)
    actime <- system.time(acfit <- lmmv(x = x, y =y)$beta[2])[3]
    c(lmfit, amfit, acfit,cctime,amtime,actime)
  }
  tmp <- t(replicate(nrep, onerep(n,p,NAprob)))
  tmp <- as.data.frame(tmp)
  cat("Mean of normal linear model:",mean(tmp[,1]),"\n")
  cat("Mean of Amelia 2:",mean(tmp[,2]),"\n")
  cat("Mean of Available Cases:",mean(tmp[,3]),"\n")
  cat("Variance of normal linear model:",var(tmp[,1]),"\n")
  cat("Variance of Amelia 2:",var(tmp[,2]),"\n")
  cat("Variance of Available Cases:",var(tmp[,3]),"\n")
  cat("Time for complete obeservation is", sum(tmp[,4]),"\n")
  cat("Time for Amelia is", sum(tmp[,5]),"\n")
  cat("Time for Available Cases is", sum(tmp[,6]),"\n")
}

#simlm(100, 10000, 3, 0.1)

# Cor: Use correlation or covariance matrix

simeig <- function(nreps,n,p,NAprob,Cor = FALSE) {
  onerep <- function(n,p,NAprob) {
    np <- n * p
    x <- matrix(runif(np),ncol=p)
    x[,p] <- x[,p] + x[,-p] %*% rep(1,p-1)
    idxs <- makeNA(n = np, missp = NAprob)
    x[idxs] <- NA
    if(!Cor)
    {
      mvtime <- system.time(eig1mv <- sqrt(PCAmv(data = x)$eigenvalues[1]))[3]
      cctime <- system.time({corcc <- cov(x = x,use='complete')
                            eig1cc <- sqrt(eigen(corcc)$values[1])})[3]
      amtime <- system.time({am <- amelia(x = x, p2s = 0 , m = 5)$imputations
                            amX <- Reduce("+", am)/length(am)
                            #using complete because there may still be complete missing pairs after imputations
                            amMat <- cov(x = amX,use='complete')
                            eig1am <- sqrt(eigen(amMat)$values[1])})[3]
      c(eig1cc,eig1am,eig1mv,cctime,amtime,mvtime)
    }
    else
    {
      mvtime <- system.time(eig1mv <- sqrt(PCAmv(data = x,Cor = TRUE)$eigenvalues[1]))[3]
      cctime <- system.time({corcc <- cor(x = x,use='complete')
                            eig1cc <- sqrt(eigen(corcc)$values[1])})[3]
      amtime <- system.time({am <- amelia(x = x, p2s = 0 , m = 5)$imputations
                            amX <- Reduce("+", am)/length(am)
                            #using complete because there may still be complete missing pairs after imputations
                            amMat <- cor(x = amX,use='complete')
                            eig1am <- sqrt(eigen(amMat)$values[1])})[3]
      c(eig1cc,eig1am,eig1mv,cctime,amtime,mvtime)
    }
  }
  tmp <- t(replicate(nreps, onerep(n,p,NAprob)))
  tmp <- as.data.frame(tmp)
  cat("The largest Eigenvalue simulation for PCA","\n")
  cat("Mean of complete obeservations",mean(tmp[,1]),"\n")
  cat("Mean of Amelia 2:",mean(tmp[,2]),"\n")
  cat("Mean of Available Cases:",mean(tmp[,3]),"\n")
  cat("Variance of complete obeservations:",var(tmp[,1]),"\n")
  cat("Variance of Amelia 2:",var(tmp[,2]),"\n")
  cat("Variance of Available Cases:",var(tmp[,3]),"\n")
  cat("Time for complete obeservation is", sum(tmp[,4]),"\n")
  cat("Time for Amelia is", sum(tmp[,5]),"\n")
  cat("Time for Available Cases is", sum(tmp[,6]),"\n")
}

#simeig(100,1000,3,0.1,Cor = FALSE)
#simeig(100,1000,3,0.1,Cor = TRUE)


Eig1mv <- function(x, Cor = FALSE)
{
  sqrt(PCAmv(data = x,Cor = TRUE)$eigenvalues[1])
}

Eig1cc <- function(x)
{
  corcc <- cor(x = x,use='complete')
  sqrt(eigen(corcc)$values[1])
}

Eig1am <- function(x)
{
  am <- amelia(x = x, p2s = 0 , m = 5)$imputations
  amX <- Reduce("+", am)/length(am)
  #using complete because there may still be complete missing pairs after imputations
  amMat <- cor(x = amX,use='complete')
  sqrt(eigen(amMat)$values[1])
}