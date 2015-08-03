#' lmmv
#'
#' This function uses Available Cases method to conduct linear regression for numeric data with missing value.
#' @param x numeric matrix or numeric data.frame. The first column must be 1, the constant.
#' @param y response for linear regression
#' @return A list contains estimates of coefficient and variances of estimates.
#' @export
#' @examples
#' n = 1000
#' x1 <- runif(n)
#' x2 <- runif(n)
#' x3 <- runif(n)
#' y <- x1 + x2 + x3 + runif(n)
#' y <- x1 + x2 + x3 + runif(n)
#' makeNA <- function(k,missp) {
#'    nmiss <- rbinom(1,k,missp)
#'    sample(1:k,nmiss,replace=FALSE)
#' }
#' 
#' idx <- makeNA(3*n, 0.1)
#' x <- cbind(x1,x2,x3)
#' x[idx] <- NA
#' idy <- makeNA(n, 0.1)
#' y[idy] <- NA
#' #The first column of the Data matrix is 1.
#' #By Now, we only provide the model includes the intersect.
#' x <- cbind(1,x)
#' lmmv(x,y)
require(numDeriv)
lmmv <- function(x,y)
{ 
  x <- as.matrix(x)
  n <- dim(x)[1]
  mcol <- dim(x)[2]
  # length of parameter length
  parlen <- mcol*(mcol+1)/2 - 1 + mcol
  # compute the estimate of E(X1), E(X1X2) etc. 
  xest <- unlist(sapply(1:(mcol-1),function(i) colMeans(x[,i]* x[,i:mcol],na.rm = TRUE)))
  xest <- c(xest[-1], mean(x[,mcol] * x[,mcol],na.rm = TRUE))
  yest <- colMeans(y*x, na.rm = T)
  ests <- c(xest, yest)
  
  # theta: parameters(ests)
  # mcol: the column of x
  # p: the pth beta in beta vector 1 <= p <= m
  betahat = function(theta, mcol , pars , p)
  {
    D = matrix(c(theta[(pars - mcol + 1):pars]), nrow = mcol)
    B = matrix(1,nrow = mcol, ncol = mcol)
    B[1,2:mcol] = theta[1:(mcol-1)]
    B[2:mcol,1] = theta[1:(mcol-1)]
    k = mcol
    for(i in 2:mcol)
    {
      for(j in i:mcol)
        {
          B[j,i] <- B[i,j] <- theta[k]
          k = k + 1
        }
    }
    solve(B)[p,] %*% D
  }  
  # compute the covariance matrix of ests
  
  cvdata = matrix(0, nrow = n, ncol = parlen)
  cvdata[,1:(mcol - 1)] = x[,1] * x[,2:mcol]
  k = mcol
  for(i in 2:mcol)
  {
    cvdata[,k : (k + mcol - i)] = x[,i] * x[,(i:mcol)]
    k = k + (mcol - i) + 1
  }
  cvdata[,k: (k + mcol - 1)] = x * y
  k = k + mcol - 1
  
  cvxy <- cov(cvdata,use = "pairwise.complete.obs")
  
  for(i in 1:parlen)
  {
    for(j in i:parlen)
    {
      #browser()
      ni <- as.double(sum(!is.na(cvdata[,i])))
      nj <- as.double(sum(!is.na(cvdata[,j])))
      nij<- as.double(sum(!is.na(cvdata[,i] * cvdata[,j])))
      
      cvxy[j,i] <- cvxy[i,j] <- cvxy[i,j] * (nij / (ni * nj))
    }
  }
  beta <- sapply(1:mcol, function(i) 
    betahat(theta = ests, mcol = mcol, pars = parlen, p = i))
  
  derivs <- sapply(1:mcol, function(i) 
    genD(betahat, ests, mcol = mcol, pars = parlen, p = i)$D[1:parlen])
  
  vars <- diag(t(derivs) %*% cvxy %*% derivs)
  
  list(beta = beta, vars = vars)
}






