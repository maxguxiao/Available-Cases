% Generated by roxygen2 (4.1.1.9000): do not edit by hand
% Please edit documentation in R/PCAmv.R
\name{PCAmv}
\alias{PCAmv}
\title{PACmv}
\usage{
PCAmv(data, retX = TRUE)
}
\arguments{
\item{data}{numeric matrix or numeric data.frame.}

\item{retX}{Whether return eigenvectors or not. Defaults to TRUE.}
}
\value{
A list contains eigenvlue, eigenvector(optional)
}
\description{
This function uses Available Cases method to conduct Principle Components Analysis for numeric data with missing value.
}
\examples{
X = runif(6000)
idx = sample(600,60,replace = FALSE)
X[idx] = NA
data = matrix(X, nrow = 1000)
PCAmv(dat)
}

