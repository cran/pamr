\name{pamr.adaptthresh}
\alias{pamr.adaptthresh}
\title{ A function  to adaptive choose threshold scales, for
  use in pamr.train}
\description{A function  to adaptive choose threshold scales, for
  use in pamr.train}

\usage{
pamr.adaptthresh(object, ntries = 10, reduction.factor = 0.9, full.out = FALSE)
}

\arguments{
  
  \item{object}{The result of a call to pamr.train }
  \item{ntries}{Number of iterations to use in algorithm}
  \item{reduction.factor}{Amount by which a scaling is reduced in one
    step of the algorithm}
\item{full.out}{Should full output be returned? Default FALSE}
}    


\details{
  \code{pamr.adaptthresh} Adaptively searches for  set of
  good threshold scales.  The baseline (default) scale is 1 for
  each class. The idea is that for easy to classify classes,
  the threshold scale can be increased without increasing the error
  rate for that class, and resulting in fewer genes needed for the
  classification rule. The scalings from pamr.adaptthresh are then used
  in pamr.train, and pamr.cv. The results may be better than those obtained
  with the default values of threshold.scale.
}


\references{
  
  Robert Tibshirani, Trevor Hastie, Balasubramanian Narasimhan, and Gilbert Chu. "Diagnosis of multiple cancer types by shrunken centroids of gene expression"  PNAS 2002 99:6567-6572 (May 14). 
  
  Robert Tibshirani,  Trevor Hastie, Balasubramanian Narasimhan, and Gilbert Chu (2002).
  Class prediction by  
  nearest shrunken centroids,with applications
  to DNA microarrays. Stanford tech report.
}




\author{ Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{
suppressWarnings(RNGversion("3.5.0"))
set.seed(120)
x <- matrix(rnorm(1000*20),ncol=20)
y <- sample(c(1:4),size=20,replace=TRUE)
mydata <- list(x=x,y=y)
mytrain <-   pamr.train(mydata)
new.scales <- pamr.adaptthresh(mytrain)

 
mytrain2 <- pamr.train(mydata, threshold.scale=new.scales)

myresults2 <- pamr.cv(mytrain2, mydata)

}

