\name{order.restricted}
\alias{order.restricted}
\title{Compute the order of a vector subject to restrictions. }
\description{
  When data have been clustered there are generally many different
  orderings of the rows or columns of the data matrix that are
  consistent with the clustering information.

  This function takes a set of labels, a set of weights and a set of
  clusters. It returns the labels ordered by increasing weight within
  clusters. The clusters are ordered by average within cluster weight.
}
\usage{
order.restricted(labels, weights, clusters)
}
\arguments{
  \item{labels}{The labels, or indices to be ordered. }
  \item{weights}{A vector of weights to be used for the ordering. }
  \item{clusters}{A vector indicating cluster membership. }
}
\details{
  Currently all three input vectors are assumed to be the same
  length. Eventually we would like to extend the functionality of this
  routine so that it can handle inputs of class \code{\link{dendrogram}}.

  This is a helper function for plotting microarray data. We want to
  reorder the data so that genes with similar levels of expression (ie
  those that cluster together) get plotted together.

  Other distance measures etc may prove useful.
}
\value{
   The input \code{labels} ordered by weights within clusters.
}

\author{R. Gentleman }

\seealso{ \code{\link{plot}} }

\examples{
   library(genefilter)
   library(mva)
   data(tdata)
   tdata <- as.matrix(tdata)
   km1 <- kmeans(t(tdata), 3, iter.max=100)
   km2 <- kmeans(tdata, 6, iter.max=100)
   rmns <- apply(tdata, 1, mean)
   cmns <- apply(tdata, 2, mean)
   dn <- dimnames(tdata)
   ord.cols <- order.restricted(dn[[2]], cmns, km1$clust)
   ord.rows <- order.restricted(dn[[1]], rmns, km2$clust)
   tdata <- tdata[ord.rows, ord.cols]
}
\keyword{ manip }

