#' screen plot BML object
#' 
#' screen plot BML object
#' @param ACPCA ACpac object
#' @param npc Number of principle components to display.
#' @param use defalut to Variance of eigenvalues. 
#' You can also choose "cumvariance", which is the cumulative variance or "eigenvalue"
#' @export

plot.ACPCA <- function(ACPCA, npc = "maximum",use = "variance")
{
  eigenvalues <- ACPCA$eigenvalues
  variance <- eigenvalues * 100/ sum(eigenvalues)
  cumvar <- cumsum(variance)
  df <- data.frame(dimension = 1:length(eigenvalues),
                   eigenvalues = eigenvalues, 
                   variance = variance, 
                   cumvariance = cumvar)
  if(npc == "maximum")
  {
    npc = length(eigenvalues)
    p <- ggplot(df)
  }
  else
  {
    df <- df[1:npc,]
    p <- ggplot(df)
  }
  # screen plot of eigen value
  use = tolower(use)
  method = c("variance","cumvariance","eigenvalues")
  if(use %in% method)
  {p + geom_bar(aes_string(x = "dimension",y = use), stat = "identity", 
               fill = "lightblue") +
    geom_point(aes_string(x = "dimension",y = use), stat = "identity") + 
    geom_line(aes_string(x = "dimension",y = use), stat = "identity",colour = 'red') +
    ggtitle(paste("Screen plot of fisrt ",npc ," Eigenvalues ", sep = "" )) 
  }
  else
  {
    stop(paste("Do not know how to deal with",use,sep = " "))
  }
}


