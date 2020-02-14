\name{RIP.PWMscan}
\alias{RIP.PWMscan}
\docType{data}
\title{
Provided RIP PWM scan
}
\description{
Position weigh matrix scan of target gene promoters for the RIP classifier.
}
\usage{data(RIP.PWMscan)
}
\details{
A \code{list} containing the PWM scan results for transcription factors (TFs) available in the RIP universe. Each named item corresponds to one TF and contains a named vector of p-values (with Entrez gene IDs as the vector names). The scan is truncated at p<=0.1.
}
\source{
\url{http://www.ichip.de/software/RIP.html}
}
\references{
Bauer et al. (2011) Bioinformatics 27:2239-47 - RIP: the regulatory interaction predictor--a machine learning-based approach for predicting target genes of transcription factors.
}
\examples{
data(RIP.PWMscan)
}
\keyword{datasets}
