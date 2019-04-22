\name{pamr.plotcv}
\alias{pamr.plotcv}
\title{ A function to plot the cross-validated  error curves from
  the nearest shrunken centroid
  classifier} 
\description{A function to plot the cross-validated error curves
  the nearest shrunken centroid
  classifier}

\usage{
pamr.plotcv(fit)
}

\arguments{
  \item{fit}{The result of a call to pamr.cv}
}    

\details{
  \code{pamr.plotcv} plots the cross-validated misclassification error
  curves,   from nearest shrunken
  centroid classifier. An overall plot, and a plot by class, are produced.
}



\author{ Trevor Hastie,Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{
suppressWarnings(RNGversion("3.5.0"))
set.seed(120)
x <- matrix(rnorm(1000*20),ncol=20)
y <- sample(c(1:4),size=20,replace=TRUE)
mydata <- list(x=x,y=y)
mytrain <-   pamr.train(mydata)
mycv <-  pamr.cv(mytrain, mydata)
pamr.plotcv(mycv)
}


