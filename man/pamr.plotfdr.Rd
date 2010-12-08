\name{pamr.plotfdr}
\alias{pamr.plotfdr}
\title{ A function to plot the FDR curve from
  the nearest shrunken centroid
  classifier} 
\description{A function to plot the FDR curve
  the nearest shrunken centroid
  classifier}

\usage{
pamr.plotfdr(fdrfit,  call.win.metafile = FALSE)
}

\arguments{
  \item{fdrfit}{The result of a call to pamr.fdr}
\item{call.win.metafile}{Used by Excel interface}
}    

\details{
  \code{pamr.plotfdr} plots the FDR
  curves   from nearest shrunken
  centroid classifier. The median FDR (solid line) and upper 90 percentile
(broken line) are shown
}



\author{ Trevor Hastie,Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{
set.seed(120)
x <- matrix(rnorm(1000*20),ncol=20)
y <- sample(c(1:2),size=20,replace=TRUE)
x[1:50,y==2]=x[1:50,y==2]+3
mydata <- list(x=x,y=y)
mytrain <-   pamr.train(mydata)
myfdr <-  pamr.fdr(mytrain, mydata)
pamr.plotfdr(myfdr)
}


