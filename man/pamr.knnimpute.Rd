\name{pamr.knnimpute}
\alias{pamr.knnimpute}
\title{ A function to impute missing expression data}
\description{ A function to impute missing expression data}

   
\usage{
pamr.knnimpute(data ,k)
}

\arguments{
\item{data}{The input data. A list with components: x- an expression
        genes in the rows, samples in the columns, and y-  a vector of
        the class labels for each sample. Same form as used by pamr.train,
 and same as that produced by pamr.from.excel}

\item{k}{Number of neighbors to be used in the imputation.Default=10 }
 } 
    
\details{
    \code{pamr.knnimpute} 
used k-nearest neighbors in the space of genes to impute missing
expression values.
For large data matrices with lots of missing values, this function can take a while
to run.}

\value{
  \item{data} {The input data list, with x replaced by the 
   imputed version of x}
   }

   \references{
Olga Troyanskaya, Michael Cantor, Gavin Sherlock, Pat Brown, Trevor Hastie, Robert Tibshirani, David Botstein and Russ B. Altman, Missing value estimation methods for DNA microarrays BIOINFORMATICS Vol. 17 no. 6, 2001 Pages 520-525 
 }

   
 

\author{ Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\keyword{ }

