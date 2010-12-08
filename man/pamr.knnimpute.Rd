\name{pamr.knnimpute}
\alias{pamr.knnimpute}
\title{A function to impute missing expression data}
\description{A function to impute missing expression data, using
  nearest neighbor averaging.}


\usage{
pamr.knnimpute(data ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500)
}

\arguments{
  \item{data}{The PAM input data. A list with components: \code{x}, an expression
    matrix with 
    genes in the rows, samples in the columns, and \code{y}, a vector of
    the class labels for each sample. Same form as used by pamr.train,
    and same as that produced by \code{pamr.from.excel}}
  \item{k}{Number of neighbors to be used in the
    imputation (default=10)}
  \item{rowmax}{The maximum percent missing data allowed in any row
    (default 50\%). For any rows with more than \code{rowmax}\% missing
    are imputed using the overall mean per sample.}
  \item{colmax}{The maximum percent missing data allowed in any column
    (default 80\%). If any column has more than \code{colmax}\% missing data,
    the program halts and reports an error.}
  \item{maxp}{The largest block of genes imputed using the knn
    algorithm inside \code{pamr.knnimpute} (default
    1500); larger blocks are divided by two-means clustering
    (recursively) prior to imputation. If \code{maxp=p}, only knn
    imputation is done.}
} 
  
\details{
  \code{pamr.knnimpute} 
  uses k-nearest neighbors in the space of genes to impute missing
  expression values.
  
  For each gene with missing values, we find the k nearest neighbors using
  a Euclidean metric, confined to the columns for which that gene is NOT
  missing. Each candidate neighbor might be missing some of the
  coordinates used to calculate the distance. In this case we average the
  distance from the non-missing coordinates. Having found the k nearest
  neighbors for a gene, we impute the missing elements by averaging those
  (non-missing) elements of its neighbors. This can fail if ALL the
  neighbors are missing in a particular element. In this case we use the
  overall column mean for that block of genes.
  
  Since nearest neighbor imputation costs
  \eqn{O(p\log(p))}{O(p*log(p))} operations per gene, where p is the
  number of rows, the computational time can be excessive for large p and
  a large number of missing rows. Our strategy is to break blocks with
  more than \code{maxp} genes into two smaller blocks using two-mean
  clustering. This is done recursively till all blocks have less than
  \code{maxp} genes. For each block, knn imputation is done separately.
  We have set the default value of \code{maxp} to 1500. Depending on the
  speed of the machine, and number of samples, this number might be
  increased. Making it too small is counter-productive, because the
  number of two-mean clustering algorithms will increase.
}

\value{
  \item{data}{The input data list, with x replaced by the 
    imputed version of x}
}

\references{
  Hastie, T., Tibshirani, R., Sherlock, G., Eisen, M., Brown, P. and
  Botstein, D., Imputing Missing Data for Gene Expression Arrays,
  Stanford University Statistics Department Technical report (1999), http://www-stat.stanford.edu/~hastie/Papers/missing.pdf
  
  Olga Troyanskaya, Michael Cantor, Gavin Sherlock, Pat Brown, Trevor Hastie, Robert Tibshirani, David Botstein and Russ B. Altman, Missing value estimation methods for DNA microarrays BIOINFORMATICS Vol. 17 no. 6, 2001 Pages 520-525 
}

\author{ Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu}


