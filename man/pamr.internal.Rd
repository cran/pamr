\name{pamr-internal}
\title{Internal pamr functions}
\alias{pamr.pairscore}
\alias{pamr.pvalue.survival}
\alias{pamr.score.to.class1}
\alias{pamr.score.to.class2}
\alias{pamr.knnimpute.old}
\description{Internal pamr functions}
\usage{
pamr.pairscore(x, pair.ind=NULL)
pamr.pvalue.survival(group, survival.time, censoring.status, ngroup.survival)
pamr.score.to.class1(x, scores, cutoff=2, n.class=2)
pamr.score.to.class2(x, scores, cutoff=2, n.pc=1, n.class=2)
pamr.knnimpute.old(data, k = 10) 
}
\author{Balasubramanian Narasimhan
 and Rob Tibshirani}
\details{
  These are not to be called by the user.
}
\keyword{internal}

