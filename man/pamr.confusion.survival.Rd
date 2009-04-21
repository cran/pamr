\name{pamr.confusion.survival}
\alias{pamr.confusion.survival}
\title{ Compute confusin matrix from pamr survival fit}
\description{computes confusion matrix for (survival.time,censoring) outcome
 based on fit object "fit" and class predictions "yhat"
 soft response probabilities for (survival.time,censoring) are first estimated
  using Kaplan-Meier method applied to training data}


\usage{
pamr.confusion.survival(fit, survival.time, censoring.status, yhat)
}

\arguments{
  
  \item{fit}{The result of a call to pamr.train or pamr.cv}
  \item{survival.time}{Survival time}
  \item{censoring.status}{censoring status}
  \item{yhat}{class predictions}
  
}    

\details{}



\references{}




\author{ Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{
}
\keyword{ }

