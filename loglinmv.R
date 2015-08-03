# author: N. Matloff
#
# log-linear model; at present, handles only the 3-factor casea
#
# arguments:
#
#    x: data frame/matrix, one row per observation
#    margin: a list of vectors specifying the model, 
#            as in loglin()
#
# value:  $para component in the value emitted from R's loglin()

loglinmv <- function(x,margin) {
  # find lengths of the elements in the model, do determine what
  # situtation we are in
  termlengths <- Map(length,margin)
  n1 <- sum(termlengths == 1)  # singletons
  n2 <- sum(termlengths == 2)  # 2-way interactions
  # fully independent?
  if (n1 == 3) mdlf <- mindep else
    # X independent of Y and Z?
    if (n1 == 1 && n2 == 1) mdlf <- mxindyz else
      # Y and Z conditionally independent, given X?
      if (n1 == 1 && n2 == 2) mdlf <-myzcondindx else
        # all possible 2-way interactions
        mdlf <- mall2s
  # need an appropriate shell, with the right dimensions, labels etc.;
  # the contents here are irrelevant and will be overwritten
  x <- as.data.frame(na.omit(x))
  tbl <- table(x)
  tbl <- mdlf(x,margin,tbl,termlengths)
  loglin(tbl,margin,param=TRUE)$param
}

mindep <- function(x,margin,tbl,termlengths) {
  nc <- ncol(x)  # 3
  probs <- list()
  nvals <- vector(length=nc)
  for (i in 1:nc) {
    tmp <- table(x[,i])
    probs[[i]] <- tmp / sum(tmp)
    nvals[i] <- length(tmp)
  }
  for (i in 1:nvals[1]) 
    for (j in 1:nvals[2]) 
      for (k in 1:nvals[3]) {
        tbl[i,j,k] <- 
          probs[[1]][i] *
          probs[[2]][j] *
          probs[[3]][k] 
      }
  tbl <- nrow(x) * tbl
}

mxindyz <- function(x,margin,tbl,termlengths) {
  nc <- ncol(x)  # 3
  # which variable is X?
  ix <- which(termlengths == 1)
  # and which are Y and Z?
  iyz <- setdiff((1:nc),ix)
  probs <- list()
  nvals <- vector(length=nc)
  nvals[1] <- length(table(x[,ix]))
  nvals[2] <- length(table(x[,iyz[1]]))
  nvals[3] <- length(table(x[,iyz[2]]))
  tmp <- table(x[,ix])
  probs[[1]] <- tmp / sum(tmp)
  tmp <- table(x[,iyz])
  probs[[2]] <- tmp / sum(tmp)
  for (i in 1:nvals[1]) 
    for (j in 1:nvals[2]) 
      for (k in 1:nvals[3]) {
        tbl[i,j,k] <- 
          probs[[1]][i] *
          probs[[2]][j,k] 
      }
  tbl <- nrow(x) * tbl
}

myzcondindx <- function(x,margin,tbl,termlengths) {
  nc <- ncol(x)  # 3
  # which variable is X?
  ix <- which(termlengths == 1)
  # and which are Y and Z?
  iyz <- setdiff((1:nc),ix)
  iy <- iyz[1]
  iz <- iyz[2]
  probs <- list()
  nvals <- vector(length=nc)
  nvals[1] <- length(table(x[,ix]))
  nvals[2] <- length(table(x[,iy]))
  nvals[3] <- length(table(x[,iz]))
  tmp <- table(x[,ix])
  probs[[1]] <- tmp / sum(tmp)
  tmp <- table(x[,c(ix,iy)])
  probs[[2]] <- tmp / sum(tmp)
  tmp <- table(x[,c(ix,iz)])
  probs[[3]] <- tmp / sum(tmp)
  for (i in 1:nvals[1]) 
    for (j in 1:nvals[2]) 
      for (k in 1:nvals[3]) {
        tbl[i,j,k] <- 
          probs[[3]][i,k] *
          probs[[2]][i,j] /
          probs[[1]][i] 
      }
  tbl <- nrow(x) * tbl
}
makena <- function(k,missp) {
  nmiss <- rbinom(1,k,missp)
  sample(1:k,nmiss,replace=FALSE)
}

# converts an R table to a fake data frame; if a cell has frequency k,
# it will appear k times in the output
tbltofakedf <- function(tbl) {
  adf <- as.data.frame(tbl)
  nc <- ncol(adf)
  onecell <- function(adfrow) {
    freq <- as.numeric(adfrow[nc])
    if (freq == 0) return(NULL)
    remainingrow <- adfrow[-nc]
    matrix(rep(remainingrow,freq),byrow=TRUE,nrow=freq)
  }
  Reduce(rbind,apply(adf,1,onecell))
}

#ucb <- tbltofakedf(UCBAdmissions)
#ucb <- as.data.frame(ucb)
#nr <- nrow(ucb)
#nc <- ncol(ucb)
#idx <- makena(nc*nr,0.1)
#ucb <- as.matrix(ucb)
#ucb[idx] <- NA
#ucb <- as.data.frame(ucb)

#a <- table(ucb)
#b <- table(na.omit(ucb))
#all.equal(target = a,current = b)
##The result is TRUE for all.equal()

##Same result
#loglinmv(x = ucb,margin = list(1,2,3))
#loglin(table(ucb),list(1,2,3),param=TRUE)$param
##same result
#a <- loglinmv(ucb,list(1,2:3))
#b <- loglin(table(ucb),list(1,2:3),param=TRUE)$param
#all.equal(a,b)

#a <- loglinmv(x = ucb,margin = list(1,c(1:2),c(2:3)))
#b <- loglin(table(ucb),list(1,c(1:2),c(2:3)),param=TRUE)$param
#all.equal(a,b)
#[1] "Component “(Intercept)”: Mean relative difference: 0.04317363"
#[2] "Component “V2”: Mean relative difference: 0.5079941"          
#[3] "Component “V3”: Mean relative difference: 1.537664"           
#[4] "Component “V2.V3”: Mean relative difference: 10.08918"  


#a <- loglinmv(ucb,list(2,c(1,3)))
#b <- loglin(table(ucb),list(2,c(1,3)),param=TRUE)$param
#all.equal(a,b)


