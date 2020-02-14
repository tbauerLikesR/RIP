\name{RIP.universe}
\alias{RIP.universe}
\docType{data}
\title{
Available genes and transcription factors (TFs) in RIP
}
\description{
Provides the available RIP universe of genes and TFs that can be used to predict RIs by the RIP master classifer.
}
\usage{data(RIP.universe)}
\details{
    \code{entrez}:Contains a the Entrez gene IDs and symbols used for the RIP classifier. They are derived from the Affymetrix Human Genome U133 Set annotation data (chip hgu133a) of Bioconductor release 2.4. Therefore, rownames of the data.frame correspond to affymetrix probe IDs.
    \code{TF}:The available TF symbols. A colon ':' is used to separate subunits of TF complexes.
    \code{core.genes}:These are the genes used to estimate the performance with the gold standard.
}
\source{
Bauer et al. (2011) Bioinformatics 27:2239-47 - RIP: the regulatory interaction predictor--a machine learning-based approach for predicting target genes of transcription factors.
\url{http://www.ichip.de/software/RIP.html}
}
\references{
Affymetrix hgu133a annotation used: \url{http://www.bioconductor.org/packages/2.4/AnnotationData.html}
}
\examples{
data(RIP.universe)
}
\keyword{datasets}
