\name{pamr.predictmany}
\alias{pamr.predictmany}
\title{ A function giving prediction information for many threshold values,
  from a nearest shrunken centroid fit.}
\description{A function giving prediction information for many threshold values,
  from a nearest shrunken centroid fit}


\usage{
pamr.predictmany(fit, newx, threshold=fit$threshold, prior =fit$prior, threshold.scale = fit$
        threshold.scale, ...)
}

\arguments{
  \item{fit}{The result of a call to pamr.train }
\item{newx}{Matrix of features at which predictions are to be made}
  \item{threshold}{The desired threshold values}
  \item{prior}{Prior probabilities for each class. Default is that
    specified in "fit"}
  \item{threshold.scale}{Additional scaling factors to be applied
    to the thresholds. Vector of length equal to the number of
    classes.
    Default is that
    specified in "fit".}
\item{...}{Additional arguments to be passed to pamr.predict}
}


\details{
}


\references{ }




\author{ Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{
set.seed(120)
x <- matrix(rnorm(1000*20),ncol=20)
y <- sample(c(1:4),size=20,replace=TRUE)
mydata <- list(x=x,y=y)
mytrain <-   pamr.train(mydata)

pamr.predictmany(mytrain, mydata$x)
 
}
\keyword{ }

