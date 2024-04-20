\name{pamr-internal}
\title{Internal pamr functions}
\alias{pamr.cube.root}
\alias{pamr.pairscore}
\alias{pamr.pvalue.survival}
\alias{pamr.score.to.class1}
\alias{pamr.score.to.class2}
\alias{print.nsc}
\alias{print.nsccv}
\alias{print.pamrcved}
\alias{print.pamrtrained}

\description{Internal pamr functions}
\usage{
pamr.pairscore(x, pair.ind=NULL)
pamr.pvalue.survival(group, survival.time, censoring.status, ngroup.survival)
pamr.score.to.class1(x, scores, cutoff=2, n.class=2)
pamr.score.to.class2(x, scores, cutoff=2, n.pc=1, n.class=2)
pamr.cube.root(x)
\method{print}{nsc}(x, ...)
\method{print}{nsccv}(x, ...)
\method{print}{pamrtrained}(x, ...)
\method{print}{pamrcved}(x, ...)
}
\author{Balasubramanian Narasimhan
 and Rob Tibshirani}
\details{
  These are not to be called by the user.
}
\keyword{internal}

