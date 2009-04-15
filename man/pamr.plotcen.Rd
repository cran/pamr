\name{pamr.plotcen}
\alias{pamr.plotcen}
\title{ A function to plot the shrunken class centroids, from the nearest shrunken centroid classifier}
\description{ A function to plot the shrunken class centroids, from the
  nearest shrunken centroid classifier produced by pamr.train}


\usage{
pamr.plotcen(fit, data,  threshold)
}

\arguments{
  \item{fit}{The result of a call to pamr.train}
  \item{data}{The input data, in the same form as that used
    by pamr.train},
  \item{threshold}{The desired threshold value}
}    

\details{
  \code{pamr.plotcen} plots the shrunken class centroids for each
  class,
  for genes surviving the threshold for at least once class. If
  genenames
  are included in "data", they are added to the plot. Note: for many
  classes
  and long gene names, this plot may need some manual prettying.
}


\references{
}




\author{ Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{
set.seed(120)
x <- matrix(rnorm(1000*20),ncol=20)
y <- sample(c(1:4),size=20,replace=TRUE)
mydata <- list(x=x,y=y,genenames=as.character(1:1000))
mytrain <-   pamr.train(mydata)
mycv <- pamr.cv(mytrain,mydata)
pamr.plotcen(mytrain, mydata,threshold=1.6)
 
}
\keyword{ }

