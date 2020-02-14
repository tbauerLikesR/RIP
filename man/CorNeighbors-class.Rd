\name{CorNeighbors-class}
\Rdversion{1.1}
\docType{class}
\alias{CorNeighbors-class}
\alias{CorNeighbors}

\title{Class "CorNeighbors"
}
\description{
This is a class representation for storing objects after calculating the correlation neighbors of a series of ExpressionSets.
}
\section{Objects from the Class}{
Objects can be created by calls to the function \code{getCorNeighbors} or in the form of \code{new("CorNeighbors", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{datasets}:}{Object of class \code{"character"}. The name(s) of the processed ExpressionSet(s).}
    \item{\code{parameters}:}{Object of class \code{"numeric"}. The defined CC and FoC parameters.}
    \item{\code{entrez_annotation}:}{Object of class \code{"character"}. The provided Entrez gene ID annotation.}
    \item{\code{cor_neighbors}:}{Object of class \code{"matrix"}. Contains the correlation neighbors (Entrez gene IDs) - one per row.}
    \item{\code{mean_cor}:}{Object of class \code{"numeric"}. the mean correlation coefficient of all gene pairs from the annotation. The order is equvalent to the \code{upper.tri} of the correlation matrix.}
  }
}
\section{Methods}{
No methods defined with class "CorNeighbors" in the signature.
}
\references{
Bauer T et al. (2011) Bioinformatics 27:2239-47 - RIP: the regulatory interaction predictor--a machine learning-based approach for predicting target genes of transcription factors.
}
\author{Tobias Bauer \email{t.bauer@dkfz.de}}
\seealso{\code{\link{getCorNeighbors}}
}
\keyword{classes}
