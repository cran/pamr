\name{pamr.to.excel}
\alias{pamr.to.excel}
\title{ A function to write out a data object into a tab-delimited text file}
\description{ A function to write out a data object into a tab-delimited text file}


\usage{
pamr.to.excel(data,file,trace=TRUE)
}

\arguments{
  \item{data}{A data object, of the same form as is read in by pamr.from.excel.
    Must have components x (the matrix of expression values), y- a vector of
    class labels for each sample, geneid- a vector of gene identifiers and
    genenames- a vector of gene names. Optional components: samplelabels and
    batchlabels, both character vectors}
\item{file}{Character name of a text file.}
\item{trace}{Optional argument. If true, progress in writing out file is reported.}
}
\details{
  \code{pamr.to.excel} 
    writes out the data object into a tab-delimited text file, of the same
    form as is read in by pamr.from.excel.
    Useful for writing out data that has been imputed by pamr.knnimpute or
    adjusted by pamr.batchadjust. Note- this function writes the file out
    one line at a time, and hence can take a while for big datasets.
}




\author{ Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

