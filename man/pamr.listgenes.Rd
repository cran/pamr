\name{pamr.listgenes}
\alias{pamr.listgenes}
\title{ A function to list the genes that survive the thresholding,
  from the nearest shrunken centroid classifier}
\description{ A function to list the genes that survive the thresholding,
  from the nearest shrunken centroid classifier produced by pamr.train}

\usage{
pamr.listgenes(fit, data,  threshold, fitcv=NULL, genenames=FALSE)
}

\arguments{
  \item{fit}{The result of a call to pamr.train}
  \item{data}{The input data.  In the same format as the input data for
    pamr.train}
  
  \item{threshold}{The desired threshold value}
\item{fitcv}{Optional object, result of a call to pamr.cv}
  \item{genenames}{Include genenames in the list? If yes, they
    are taken from "data". Default is false (geneid is always included
    in the list).}
}    

\details{
  \code{pamr.listgenes} List the geneids, and standardized centroids for
  each class, for genes surviving at the given threshold.
If fitcv is provided, the function also reports
the average rank of the gene in the
cross-validation folds, and the proportion of times  that the gene is chosen
(at the given threshold) in the cross-validation folds.
}



\author{ Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{

#generate some data
set.seed(120)
x <- matrix(rnorm(1000*20),ncol=20)
y <- sample(c(1:4),size=20,replace=TRUE)

mydata <- list(x=x,y=factor(y), geneid=as.character(1:nrow(x)), genenames=paste("g",as.character(1:nrow(x)),sep=""))


#train classifier
mytrain<-   pamr.train(mydata)

pamr.listgenes(mytrain, mydata, threshold=1.6)
 
}

