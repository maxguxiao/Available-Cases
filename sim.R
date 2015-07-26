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
    lmfit <- lm(y ~ x)$coefficients[2]
    am <- amelia(x = x, m = 5,p2s = 0)$imputations
    amX <- Reduce("+", am)/length(am)
    amfit <- lm(y ~ amX)$coefficients[2]
    x <- cbind(1,x)
    acfit <- lmmv(x = x, y =y)$beta[2]
    c(lmfit, amfit, acfit)
  }
  tmp <- t(replicate(nrep, onerep(n,p,NAprob)))
  tmp <- as.data.frame(tmp)
  cat("Mean of normal linear model:",mean(tmp[,1]),"\n")
  cat("Mean of Amelia 2:",mean(tmp[,2]),"\n")
  cat("Mean of Available Cases:",mean(tmp[,3]),"\n")
  cat("Variance of normal linear model:",var(tmp[,1]),"\n")
  cat("Variance of Amelia 2:",var(tmp[,2]),"\n")
  cat("Variance of Available Cases:",var(tmp[,3]),"\n")
}

#simlm(500, 10000, 3, 0.1)


simeig <- function(nreps,n,p,NAprob) {
  onerep <- function(n,p,NAprob) {
    np <- n * p
    x <- matrix(runif(np),ncol=p)
    x[,p] <- x[,p] + x[,-p] %*% rep(1,p-1)
    idxs <- makeNA(n = np, missp = NAprob)
    x[idxs] <- NA
    eig1mv <- sqrt(PCAmv(data = x)$eigenvalues[1])
    corcc <- cov(x = x,use='complete')
    eig1cc <- sqrt(eigen(corcc)$values[1])
    am <- amelia(x = x, p2s = 0 , m = 5)$imputations
    amX <- Reduce("+", am)/length(am)
    #using complete because there may still be complete missing pairs after imputations
    amMat <- cov(x = amX,use='complete')
    eig1am <- sqrt(eigen(amMat)$values[1])
    c(eig1cc,eig1am,eig1mv)
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
}

#simeig(100,1000,3,0.1)


