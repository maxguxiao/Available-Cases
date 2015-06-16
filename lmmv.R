# columns of x, the first column of x should should be 1.

lmmv <- function(x,y)
{ 
  n <- dim(x)[1]
  mcol <- dim(x)[2]
  # length of parameter length
  parlen <- mcol*(mcol+1)/2 - 1 + mcol
  # compute the estimate of E(X1), E(X1X2) etc. 
  ests <- numeric(parlen)
  k <- 0
  for(i in 1:mcol)
  {
    for(j in i:mcol)
    {
      ests[k] <- mean(x[,i] * x[,j], na.rm = T)
      k <- k + 1
    }
  }
  for(i in 1:mcol)
  {
    ests[k] <- mean(x[,i] * y, na.rm = TRUE)
    k = k + 1
  }
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
  # compute the covariance matrix
  
  cvdata = matrix(0, nrow = n, ncol = parlen)
  k = 0
  for(i in 1:mcol)
  {
    for(j in i:mcol)
    {
      #browser()
      cvdata[,k] = x[,i] * x[,j]
      k = k + 1
    }
  }
  for(i in 1:mcol)
  {
    #browser()
    cvdata[,k] = x[,i] * y
    k = k + 1
  }
  cvxy <- cov(cvdata,use = "pairwise.complete.obs")
  
  for(i in 1:parlen)
  {
    for(j in i:parlen)
    {
      #browser()
      ni <- sum(!is.na(cvdata[,i]))
      nj <- sum(!is.na(cvdata[,j]))
      nij<- sum(!is.na(cvdata[,i] * cvdata[,j]))
      cvxy[j,i] <- cvxy[i,j] <- cvxy[i,j] * nij / (ni * nj)
    }
  }
  beta <- sapply(1:mcol, function(i) 
    betahat(theta = ests, mcol = mcol, pars = parlen, p = i))
  
  derivs <- sapply(1:mcol, function(i) 
    genD(betahat, ests, mcol = mcol, pars = parlen, p = i)$D[1:parlen])
  
  vars <- diag(t(derivs) %*% cvxy %*% derivs)
  
  list(beta = beta, vars = vars)
}

#try three variables
n = 1000
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y <- runif(n)
idx1 = sample(n, round(0.1 * n), replace = FALSE)
idx2 = sample(n, round(0.1 * n), replace = FALSE)
idx3 = sample(n, round(0.1 * n), replace = FALSE)
idy = sample(n, round(0.1 * n), replace = FALSE)
x1[idx1] = NA
x2[idx2] = NA
x3[idx3] = NA
y[idy] = NA    
x = cbind(1,x1,x2,x3)
lmmv(x,y)

#simulation
simmv <- function(num,n,nreps) {
  #numm is the num of variables 
  res <- replicate(nreps,{
    x <- matrix(1,nrow = n, ncol = num + 1)# the first column is 1
    for(i in 2:(num+1))
    {
      x[,i] = runif(n)
      idx = sample(n, round(0.1 * n), replace = FALSE)
      x[idx,i] = NA
    }
    y <- runif(n)
    idy = sample(n, round(0.1 * n), replace = FALSE)
    y[idy] = NA    
    tmp <- lmmv(x,y)
    beta0 = (tmp$beta[1] - 0.5) / sqrt(tmp$vars[1])
    beta = sapply(2:(num+1), function(i) (tmp$beta[i] - 0) / sqrt(tmp$vars[i]))
    c(beta0,beta)
  })
  sapply(1:(num+1), function(i) mean(res[i,] < 1.28))
}

system.time(simmv(3,1000,200))


