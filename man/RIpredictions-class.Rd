\name{RIpredictions-class}
\Rdversion{1.1}
\docType{class}
\alias{RIpredictions-class}
\alias{RIpredictions}

\title{Class "RIpredictions"
}
\description{
This is a class representation for the results of potential regulatory interactions (RIs) classification by the RIP classifier.
}
\section{Objects from the Class}{
Objects can be created by calls to the function \code{classifyRIs} or in the form of \code{new("RIpredictions", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{predictions}:}{Object of class \code{"data.frame"}. Contains the positively classifed RIs: TF symbol, target gene Entrez ID and symbol, the percentage of positive votes from the RIP classifier, and the estimated true positive rate taken from the RIP paper (see references).}
    \item{\code{percent.votes}:}{Object of class \code{"numeric"}. The user-defined vote cutoff in percent.}
    \item{\code{CorNeighbors}:}{Object of class \code{"character"}. The name of the CorNeighbors input used for feature calculation.}
    \item{\code{PWMscan}:}{Object of class \code{"character"}. The name of the PWMscan input used for feature calculation.}
  }
}
\section{Methods}{
No methods defined with class "RIpredictions" in the signature.
}
\references{
Bauer et al. (2011) Bioinformatics 27:2239-47 - RIP: the regulatory interaction predictor--a machine learning-based approach for predicting target genes of transcription factors.
}
\author{Tobias Bauer \email{t.bauer@dkfz.de}}
\seealso{\code{\link{classifyRIs}}
}
\keyword{classes}
