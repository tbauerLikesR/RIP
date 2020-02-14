\name{RIfeatures-class}
\Rdversion{1.1}
\docType{class}
\alias{RIfeatures-class}
\alias{RIfeatures}

\title{Class "RIfeatures"
}
\description{
This is a class representation for storing features of potential regulatory interactions (RIs) for the RIP classifier.
}
\section{Objects from the Class}{
Objects can be created by calls to the function \code{calculateRIPfeatures} or in the form of \code{new("RIfeatures", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{features}:}{Object of class \code{"matrix"}. The resulting feature values required by the RIP classifier.}
    \item{\code{pRIs}:}{Object of class \code{"data.frame"}. The potential RIs.}
    \item{\code{CorNeighbors}:}{Object of class \code{"character"}. The name of the CorNeighbors input used for feature calculation.}
    \item{\code{PWMscan}:}{Object of class \code{"character"}. The name of the PWMscan input used for feature calculation.}
  }
}
\section{Methods}{
No methods defined with class "RIfeatures" in the signature.
}
\references{
Bauer et al. (2011) Bioinformatics 27:2239-47 - RIP: the regulatory interaction predictor--a machine learning-based approach for predicting target genes of transcription factors.
}
\author{Tobias Bauer \email{t.bauer@dkfz.de}}
\seealso{\code{\link{calculateRIPfeatures}},\code{\link{createPotentialRIs}}
}
\keyword{classes}
