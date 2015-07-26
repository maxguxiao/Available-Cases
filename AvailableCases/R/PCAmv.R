#' PACmv
#'
#' This function uses Available Cases method to conduct Principle Components Analysis for numeric data with missing value.
#' @param data numeric matrix or numeric data.frame.
#' @param retX Whether return eigenvectors or not. Defaults to TRUE.
#' @param Cor use correlation matrix or covariance matrix
#' @return A list contains eigenvlue, eigenvector(optional)
#' @export
#' @examples
#' X = runif(6000)
#' idx = sample(600,60,replace = FALSE)
#' X[idx] = NA
#' data = matrix(X, nrow = 1000)
#' PCAmv(data)

PCAmv <- function(data, retX = TRUE,Cor = FALSE)
{
  # Using Available Case to estimate the covariance matrix.
  if(!is.matrix(data) & !is.data.frame(data))
  {
    stop("The data is not either a matrix nor data.frame!")
  }
  n <- ncol(data)
  covar <- matrix(0,nrow = n,ncol = n)
  for(i in 1:n)
  {
    for(j in i:n)
      {
        ni <- sum(!is.na(data[,i]))
        nj <- sum(!is.na(data[,j]))
        nij <- sum(!is.na(data[,i] * data[,j]))
        covar[j,i] <- covar[i,j] <- mean(data[,i] * data[,j],na.rm = T) -
                                    mean(data[,i], na.rm = T) * mean(data[,j],na.rm = T)
      }
  }
  if(Cor)
  {
    correl <- matrix(0, ncol =n, nrow = n)
    for(i in 1:n)
    {
      for(j in i:n)
      {
        correl[i,j] <- correl[j,i] <- covar[i,j]/(sqrt(covar[i,i])*sqrt(covar[j,j]))
      }
    }
    covar <- correl
  }
  # Apply eigen() to the covariance matrix.
  eigens <- eigen(covar)
  eigenvalues <- eigens$values
  eigenvectors <-eigens$vectors
  if(retX)
  {
    ACpca <- list(eigenvalues = eigenvalues, Rotation = eigenvectors)
  }
  else
  {
    ACpca <- list(eigenvalues = eigenvalues)
  }
  class(ACpca) <- append("ACPCA", class(ACpca))
  ACpca  
}




