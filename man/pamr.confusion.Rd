\name{pamr.confusion}
\alias{pamr.confusion}
\title{ A function giving a table of true versus predicted values,
  from a nearest shrunken centroid fit.}
\description{ A function giving a table of true versus predicted values,
  from a nearest shrunken centroid fit.}

\usage{
pamr.confusion(fit, threshold, extra=TRUE)
}

\arguments{
  
  \item{fit}{The result of a call to pamr.train or pamr.cv}
  \item{threshold}{The desired threshold value}
\item{extra}{Should the  classwise and overall error rates be returned?
Default TRUE}
  
  
}    

\details{
  \code{pamr.confusion} Gives a cross-tabulation of true versus
  predicted classes for the fit returned by pamr.train or pamr.cv,
  at the specified threshold.
}


\references{}




\author{ Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{
set.seed(120)
x <- matrix(rnorm(1000*20),ncol=20)
y <- sample(c(1:4),size=20,replace=TRUE)
mydata <- list(x=x,y=y)
mytrain <-   pamr.train(mydata)
mycv <- pamr.cv(mytrain,mydata)
pamr.confusion(mytrain,  threshold=2)
pamr.confusion(mycv,  threshold=2)
 
}
\keyword{ }

