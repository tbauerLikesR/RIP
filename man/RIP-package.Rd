\name{RIP-package}
\alias{RIP-package}
\alias{RIP}
\docType{package}
\title{
RIP - the regulatory interaction predictor
}
\description{
Application of the RIP classifier to predict regulatory interactions (RIs) between human transcription factors (TFs) and target genes.
}
\details{
\tabular{ll}{
Package: \tab RIP\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2012-06-13\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
Depends: \tab R (>= 2.12.0), utils (>= 2.12.0), methods (>= 2.12.0), Biobase (>= 2.10.0)\cr
}
The package functions are designed to be applied in the following order:
\tabular{ll}{
createCorData()\tab Wrapper to calculate Pearson correlation coefficients between features of a series of ExpressionSets.\cr
getCorNeighbors()\tab Derive 'correlation neighbors' from the correlation data.\cr
createPotentialRIs()\tab Create RI candidates for defined TFs and target genes.\cr
calculateRIPfeatures()\tab Calculate features used by the RIP classifier for RI candidates.\cr
classifyRIs()\tab Apply the RIP classifier to predict RIs based on the provided features.\cr
}
}
\author{
Tobias Bauer, Roland Eils, and Rainer Koenig

Maintainer: Tobias Bauer <t.bauer@dkfz.de>
}
\references{
Bauer et al. (2011) Bioinformatics 27:2239-47 - RIP: the regulatory interaction predictor--a machine learning-based approach for predicting target genes of transcription factors.
}
\keyword{ package }
\seealso{
\url{http://www.ichip.de/software/RIP.html}
}
