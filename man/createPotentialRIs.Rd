\name{createPotentialRIs}
\alias{createPotentialRIs}
\title{Create regulatory interaction candidate sets from transcription factors and target genes
}
\description{
Creates a set of regulatory interaction (RI) candidates from a list input of transcription factors (TFs) and target genes.
}
\usage{
createPotentialRIs(target.genes=NULL, TFs=NULL, annotation.targets=c("entrezID", "symbol"))
}
\arguments{
  \item{target.genes}{character. A vector of either Entrez IDs or symbols to be used as target genes for RI feature calculation and prediction.}
  \item{TFs}{character. A vector of TF symbols to be used as regulators for RI feature calculation and prediction.}
  \item{annotation.targets}{character. The annotation type of the provided target genes: either 'entrezID' or 'symbol'.}
}
\value{a \code{data.frame} with two columns named 'TF' and 'gene'. Contains all possible combinations of the the input TFs and target genes. The data is checked to match the TFs and genes available for RIP as provided in \code{RIP.universe}.
}
\author{Tobias Bauer \email{t.bauer@dkfz.de}}
\examples{
## using Entrez gene IDs
TFs <- RIP.universe$TF[1:10]
targets <- sample(RIP.universe$entrez$ID, size=100)

potentialRIs <- createPotentialRIs(target.genes=targets, TFs=TFs, annotation.targets="entrezID")

## using Entrez gene symbols
TFs <- RIP.universe$TF[1:10]
targets <- sample(RIP.universe$entrez$symbol, size=100)

potentialRIs <- createPotentialRIs(target.genes=targets, TFs=TFs, annotation.targets="symbol")
}
