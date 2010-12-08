\name{pamr.geneplot}
\alias{pamr.geneplot}
\title{ A function to  plot the genes that surive the thresholding
  from the nearest shrunken centroid classifier}
\description{ A function to plot the genes that survive the thresholding,
  from the nearest shrunken centroid classifier produced by pamr.train}

\usage{
pamr.geneplot(fit, data, threshold)
}

\arguments{
  \item{fit}{The result of a call to pamr.train}
  \item{data}{The input data.  In the same format as the input data for
    pamr.train}
  
  \item{threshold}{The desired threshold value}
}    

\details{
  \code{pamr.geneplot} Plots the raw gene expression for genes
  that survive the specified threshold. Plot is stratified by class. 
  Plot is set up to display only up to about 20 or 25 genes,
  otherwise it gets too crowded. Hence threshold should be chosen to
  yield at most about 20 or 25 genes.
}






\author{ Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{
set.seed(120)
x <- matrix(rnorm(1000*20),ncol=20)
y <- sample(c(1:4),size=20,replace=TRUE)
mydata <- list(x=x,y=y)
mytrain <-   pamr.train(mydata)
pamr.geneplot(mytrain, mydata, threshold=1.6)
 
}

