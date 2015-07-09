PCAmv <- function(dat)
{
  dimen <- dim(dat)
  m <- dimen[1]
  n <- dimen[2]
  covar <- matrix(0,nrow = n,ncol = n)
  for(i in 1:n)
  {
    for(j in i:n)
      {
        ni <- sum(!is.na(dat[,i]))
        nj <- sum(!is.na(dat[,j]))
        nij <- sum(!is.na(dat[,i] * dat[,j]))
        covar[j,i] <- covar[i,j] <- mean(dat[,i] * dat[,j],na.rm = T) -
                                    mean(dat[,i], na.rm = T) * mean(dat[,j],na.rm = T)
          #sum((dat[,i] - sum(dat[,i],na.rm = TRUE)/ni) * 
          #                          (dat[,j] - sum(dat[,j],na.rm = TRUE)/nj),na.rm = TRUE)/(nij-1)
                                          
                                          
        #covar[j,i] <- covar[i,j] <- sum(dat[,i]*dat[,j],na.rm = T)/nij - 
        #              ((sum(dat[,i],na.rm = T)/ni) * 
        #              (sum(dat[,j],na.rm = T)/nj))
        
      }
  }
  eigens <- eigen(covar)
  eigenvalues <- eigens$values
  eigenvectors <-eigens$vectors
  list(StandardDeviation = eigenvalues, Rotation = eigenvectors)
       
}


#X = runif(6000)
#dat = matrix(X, nrow = 1000)
#PCAmv(dat)
#prcomp(dat)

#idx = sample(600,60,replace = FALSE)
#X[idx] = NA




