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
  beta = bhat(c(xbar, x2bar, ybar, xybar))
  list(beta = beta, vars = list(beta0 = var0, beta1 = var1))
}




