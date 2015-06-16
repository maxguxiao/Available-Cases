# one variable case
# x is the covariate, y is the response

lmmv.simple = function(x, y)
{
  xbar <- mean(x, na.rm = TRUE)
  x2bar <- mean(x^2, na.rm = TRUE)
  ybar <- mean(y, na.rm = TRUE)
  xybar <- mean(x*y, na.rm = TRUE)
  b0hat <- function(theta)
  {
    B = matrix(c(1,theta[1],theta[1],theta[2]),nrow = 2,byrow = TRUE)
    D = matrix(c(theta[3],theta[4]),nrow = 2)
    solve(B)[1,] %*% D
  }
  b1hat <- function(theta)
  {
    B = matrix(c(1,theta[1],theta[1],theta[2]),nrow = 2,byrow = TRUE)
    D = matrix(c(theta[3],theta[4]),nrow = 2)
    solve(B)[2,] %*% D
  }
  derivs0 <- genD(b0hat, c(xbar, x2bar, ybar, xybar))$D[1:4]
  derivs1 <- genD(b1hat, c(xbar, x2bar, ybar, xybar))$D[1:4]
  cvxy <- cov(cbind(x,x^2,y,x*y),use="pairwise.complete.obs") 
  nx <- sum(!is.na(x))
  ny <- sum(!is.na(y))
  nxy <- sum(!is.na(x*y))
  cvxy[1,1] <- cvxy[1,1] / nx
  cvxy[1,2] <- cvxy[1,2] / nx
  cvxy[1,3] <- cvxy[1,3] * nxy / (nx*ny)
  cvxy[1,4] <- cvxy[1,4] * nxy / (nx*nxy)
  cvxy[2,1] <- cvxy[1,2]
  cvxy[2,2] <- cvxy[2,2] / nx
  cvxy[2,3] <- cvxy[2,3] * nxy / (nx*ny)
  cvxy[2,4] <- cvxy[2,4] * nxy / (nx*nxy)
  cvxy[3,1] <- cvxy[1,3]
  cvxy[3,2] <- cvxy[2,3]
  cvxy[3,3] <- cvxy[3,3] / ny
  cvxy[3,4] <- cvxy[3,4] * nxy / (ny*nxy)
  cvxy[4,1] <- cvxy[1,4]
  cvxy[4,2] <- cvxy[2,4]
  cvxy[4,3] <- cvxy[3,4]
  cvxy[4,4] <- cvxy[4,4] /nxy 
  var0 <- t(derivs0) %*% cvxy %*% derivs0
  var1 <- t(derivs1) %*% cvxy %*% derivs1
  beta0 = b0hat(c(xbar, x2bar, ybar, xybar))
  beta1 = b1hat(c(xbar, x2bar, ybar, xybar))
  list(beta = list(beta0 = beta0, beta1 = beta1), vars = list(beta0 = var0, beta1 = var1))
}

#simulation
simmv <- function(n,nreps) {
  res <- replicate(nreps,{
    x <- runif(n)
    y <- runif(n)
    idx = sample(n, round(0.1 * n), replace = FALSE)
    idy = sample(n, round(0.1 * n), replace = FALSE)
    x[idx] = NA
    y[idy] = NA    
    tmp <- lmmv.simple(x,y)
    b0 = (tmp$beta$beta0 - 0.5) / sqrt(tmp$vars$beta0)
    b1 = (tmp$beta$beta1 - 0) / sqrt(tmp$vars$beta1)
    c(b0,b1)
  })
  
  beta0test = mean(res[1,] < 1.28)
  beta1test = mean(res[2,] < 1.28)
  list(beta0test = beta0test, beta1test = beta1test)
}






