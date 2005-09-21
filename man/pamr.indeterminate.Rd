\name{pamr.indeterminate}
\alias{pamr.indeterminate}
\title{ A function that takes estimate class probabilities
and produces a class prediction or indeterminate prediction}




\description{
 A function that takes estimate class probabilities
and produces a class prediction or indeterminate prediction}

\usage{
pamr.indeterminate(prob, mingap=0) 
}

\arguments{
\item{prob}{Estimated class probabilities, from pamr.predict with type="posterior")}
\item{mingap}{Minimum difference between highest and second highest probability. If difference is < mingap, prediction is set to  indeterminate (NA)}
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
prob<- pamr.predict(mytrain, mydata$x , threshold=1, type="posterior")
pamr.indeterminate(prob,mingap=.75)
}
\keyword{ }

