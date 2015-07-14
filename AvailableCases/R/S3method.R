#' screen plot BML object
#' 
#' screen plot BML object
#' @param ACPCA ACpac object
#' @param npc Number of principle components to display.
#' @export

plot.ACPCA <- function(ACPCA, npc = "maximum")
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
  p + geom_bar(aes(x = dimension, y = variance), stat = "identity", 
               fill = "lightblue") +
    geom_point(aes(x = dimension, y = variance), stat = "identity") + 
    geom_line(aes(x = dimension, y = variance), stat = "identity",colour = 'red') +
    ggtitle(paste("Screen plot of fisrt ",npc ," Eigenvalues ", sep = "" )) 
}


