\name{pamr.geneplot}
\alias{pamr.geneplot}
\title{ A function to plot  genes found by
                          the nearest shrunken centroid
                         classifier} 
   
\description{ function to plot  genes found by
                          the nearest shrunken centroid
                         classifier}

\usage{
pamr.geneplot(fit, data, threshold)
}

\arguments{
    \item{fit}{The result of a call to pamr.train}
    \item{data}{A data object, like that passed to pamr.train}
    \item{threshold}{Threshold for selecting which genes to plot}
}    
    
\details{
    \code{pamr.geneplot} 
plots the gene expression from a pamr fit, for genes which survive the "threshold".
The genes are plotted from strongest to weakest.
}
\value{}

   \references{}


\author{ Trevor Hastie,Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }


\keyword{ }

