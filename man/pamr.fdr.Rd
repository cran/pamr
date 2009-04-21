\name{pamr.fdr}
\alias{pamr.fdr}
\title{A function to estimate false discovery rates for the nearest shrunken centroid
  classifier}
\description{A function to estimate false discovery rates for the nearest shrunken centroid  classifier}


\usage{
pamr.fdr(trained.obj, data,  nperms=100, 
 xl.mode=c("regular","firsttime","onetime","lasttime"),xl.time=NULL, xl.prevfit=NULL)

}

\arguments{
  \item{trained.obj}{The result of a call to pamr.train}
\item{data}{Data object; same as the one  passed to pamr.train}
  \item{nperms}{Number of permutations for estimation of FDRs.  Default is 100} 
 \item{xl.mode}{Used by Excel interface}
\item{xl.time}{Used by Excel interface}
\item{xl.prevfit}{Used by Excel interface}
}    

\details{
  \code{pamr.fdr} estimates false discovery rates  for a nearest shrunken
  centroid classifier
}

\value{
  A list with components:
\item{results}{Matrix of estimates FDRs for various various threshold values.
Reported are both the median and 90th percentile of the  FDR over permutations}
\item{pi0}{The estimated proportion of genes that are null, i.e. not
significantly different}
}

\references{}


\author{ Trevor Hastie,Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{
set.seed(120)
x <- matrix(rnorm(1000*20),ncol=20)
y <- sample(c(1:4),size=20,replace=TRUE)

mydata <- list(x=x,y=factor(y), geneid=as.character(1:nrow(x)),
               genenames=paste("g", as.character(1:nrow(x)), sep=""))

mytrain <-   pamr.train(mydata)
myfdr <- pamr.fdr(mytrain, mydata)
}

\keyword{ }

