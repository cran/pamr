\name{pamr.batchadjust}
\alias{pamr.batchadjust}
\title{ A function to mean-adjust microarray data by batches}
\description{
  A function to mean-adjust microarray data by batches
}
\usage{
pamr.batchadjust(data)
}

\arguments{
  \item{data}{The input data. A list with components: x- an expression
    genes in the rows, samples in the columns, and y-  a vector of
    the class labels for each sample, and batchlabels- a vector of batch
    labels for each sample.}
  This object if the same form as that produced by pamr.from.excel.}


\details{
  \code{pamr.batchadjust} 
  does a genewise one-way ANOVA adjustment for expression values.
  Let x(i,j) be the expression for gene i in sample j.
  Suppose sample j in in batch b, and let B be the set of all samples
  in batch b. Then \code{pamr.batchadjust}  adjusts  x(i,j) to
  x(i,j) - mean[x(i,j)]  where the mean is taken over all samples j in B
}

\value{ A data object  of the same form as the input data, with x
  replaced by the adjusted x}




\author{ Trevor Hastie,Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{
set.seed(120)
#generate some data
x <- matrix(rnorm(1000*20),ncol=20)
y <- sample(c(1:4),size=20,replace=TRUE)
batchlabels <- sample(c(1:5),size=20,replace=TRUE)
mydata <- list(x=x,y=factor(y),batchlabels=factor(batchlabels))

mydata2 <- pamr.batchadjust(mydata)
}


