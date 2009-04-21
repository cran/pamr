\name{pamr.predict}
\alias{pamr.predict}
\title{ A function giving prediction information,
  from a nearest shrunken centroid fit.}
\description{A function giving prediction information,
  from a nearest shrunken centroid fit.}


\usage{
pamr.predict(fit, newx, threshold, type= c("class", "posterior", "centroid", "nonzero"), prior = fit$prior, threshold.scale = fit$threshold.scale)
}

\arguments{
  
  \item{fit}{The result of a call to pamr.train }
\item{newx}{Matrix of features at which predictions are to be made}
  \item{threshold}{The desired threshold value}
  \item{type}{Type of prediction desired: class predictions, posterior
    probabilities, (unshrunken) class centroids, vector of  genes
    surviving the threshold}
  \item{prior}{Prior probabilities for each class. Default is that
    specified in "fit"}
  \item{threshold.scale}{Additional scaling factors to be applied
    to the thresholds. Vector of length equal to the number of
    classes.
    Default is that
    specified in "fit".}
}


\details{
  \code{pamr.predict} Give a cross-tabulation of true versus
  predicted classes for the fit returned by pamr.train or pamr.cv,
  at the specified threshold
}


\references{ }




\author{ Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{
set.seed(120)
x <- matrix(rnorm(1000*20),ncol=20)
y <- sample(c(1:4),size=20,replace=TRUE)
mydata <- list(x=x,y=y)
mytrain <-   pamr.train(mydata)
mycv <- pamr.cv(mytrain,mydata)
pamr.predict(mytrain, mydata$x , threshold=1)
 
}
\keyword{ }

