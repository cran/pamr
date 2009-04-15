\name{pamr.from.excel}
\alias{pamr.from.excel}
\title{ A function to read in a text file saved from Excel}
\description{ A function to read in a text file saved from Excel.
  The spreadsheet is assumed to be of the format used by the SAM program.}

\usage{
pamr.from.excel(file, ncols, sample.labels = FALSE, batch.labels = FALSE) 
}

\arguments{
  \item{file}{Character name of a text file. This is assumed to be a tab-delimited text file saved from an excel spreadsheet from  "SAM". 
    The spreadsheet has one row of expression values per gene. In addition
    there is one information row and two information columns.  The first row
    has class labels for each of the samples. The first column had gene identifiers,
    and the second column has gene names. In the SAM program, for the multiclass
    option, the samples must be labelled 1,2,3 etc. Here we allow general labels,
    like "lymphoma", "colon cancer" etc 
  }
  \item{ncols}{Number of columns in file}
  \item{sample.labels}{Optional argument. If true, "file" is assumed to
    have an additional row at the top, consisting of two blank cells
    followed by a sample labels for each of the columns. If available, these sample labels
    are used by various plotting routines.} 
  \item{batch.labels}{Optional argument. If true, "file" is assumed to
    have an additional row at the top, consisting of two blank cells
    followed by a batch labels for each of the columns.
    If sample.labels=T as well, the row of batch labels are assumed to 
    come after the row of sample labels.
    The batch labels are used by the function pamr.batchadjust.
  } 
}  
\details{
  \code{pamr.from.excel} Reads in the  text file "file", and creates an
    object with components x (the matrix of expression values), y- a vector of
    class labels for each sample, geneid- a vector of gene identifiers and
    genenames- a vector of gene names
}
\value{
  A list with components
  \item{x}{the matrix of expression values}
  \item{y}{a vector of
    class labels for each sample} 
  \item{geneid}{a vector of gene identifiers}
  \item{genenames}{a vector of gene names}
  \item{samplelabels}{a vector of sample labels, if provided in "file"}
  \item{batchlabels}{a vector of batch labels, if provided in "file"}
}

\references{ }




\author{ Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\keyword{ }

