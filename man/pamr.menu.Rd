\name{pamr.menu}
\alias{pamr.menu}
\title{A function that interactively leads the user through a PAM analysis}
\description{A function that interactively leads the user through a PAM analysis}

   
\usage{
pamr.menu(data)
}

\arguments{
    \item{data}{A list with at least two components: x- an expression
        genes in the rows, samples in the columns), and y-  a vector of
        the class labels for each sample. Same form as data object
         used by pamr.train.}
}    
    
\details{
    \code{pamr.menu} provides a menu for training, cross-validating
and plotting a nearest shrunken centroid analysis.
}
\value{
   }

   \references{}


\author{ Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{
set.seed(120)
x <- matrix(rnorm(1000*20),ncol=20)
y <- sample(c(1:4),size=20,replace=TRUE)
mydata <- list(x=x,y=y)
#  pamr.menu(mydata)
}

\keyword{ }

