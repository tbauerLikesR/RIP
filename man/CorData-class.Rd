\name{CorData-class}
\Rdversion{1.1}
\docType{class}
\alias{CorData-class}
\alias{CorData}

\title{Class "CorData"
}
\description{
This is a class representation to store objects when applying Pearson correlation to a series of ExpressionSet objects.
}
\section{Objects from the Class}{
Objects can be created by calls to the function \code{createCorData} or in the form of \code{new("CorData", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{dataset}:}{Object of class \code{"character"}. The name of the processed ExpressionSet.}
    \item{\code{series}:}{Object of class \code{"character"}. An optional name of a series of ExpressionSets.}
    \item{\code{cor_coef}:}{Object of class \code{"numeric"}. The Pearson correlation coefficients of all possible feature pairs from the ExpressionSet expression values.}
    \item{\code{entrez_annotation}:}{Object of class \code{"character"}. The Entrez gene IDs of the ExpressionSet features.}
  }
}
\section{Methods}{
No methods defined with class "CorData" in the signature.
}
\author{Tobias Bauer \email{t.bauer@dkfz.de}}
\seealso{\code{\link{createCorData}}
}
\keyword{classes}
