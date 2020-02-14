\name{calculateRIPfeatures}
\alias{calculateRIPfeatures}
\alias{calculateRIPfeatures.default}
\title{Calculate ten features for RIP classification
}
\description{
Calculates the ten features required to predict regulatory interactions (RIs) between transcription factors (TFs) and target genes with the RIP classifier.
}
\usage{
calculateRIPfeatures(pRIs=NULL, cor.neighbors=NULL, show.progressbar=TRUE)
}
\arguments{
  \item{pRIs}{a \code{data.frame} with two columns named \code{TF} and \code{gene} containing potential RIs for feature calculation.}
  \item{cor.neighbors}{an object of class \code{CorNeighbors}.}
  \item{show.progressbar}{logical. Show text progress-bars of some features?}
}
\details{
It is recommended to use the function \code{createPotentialRIs} to create \code{pRIs}. Otherwise, make sure that all RI-candidates match the annotation provided in \code{RIP.universe}.

\code{cor.neighbors} is expected to be a \code{CorNeighbors} class object as returned by \code{getCorNeighbors}.

By default, a provided object is used: \code{RIP.PWMscan}. This particular list contains the PWM scan results for each gene and TF in RIP. A general cutoff of p<0.1 was applied. The names of \code{RIP.PWMscan} are that of the TFs as provided in \code{RIP.universe} and each entry contains a named vector of p-values, where the names correspond to the Entrez gene ID. Don't overwrite that variable name. (Why would you do that anyway?)
}
\value{An object of class \code{RIfeatures} which is a \code{list} of four arguments:
  \item{features}{A matrix with the resulting feature values required by the RIP classifier.}
  \item{pRIs}{A data.frame containing the pRIs input.}
  \item{CorNeighbors}{The name of the CorNeighbors input.}
  \item{PWMscan}{The name of the PWMscan input.}
}
\references{
Bauer et al. (2011) Bioinformatics 27:2239-47 - RIP: the regulatory interaction predictor--a machine learning-based approach for predicting target genes of transcription factors.
}
\author{Tobias Bauer \email{t.bauer@dkfz.de}}
\seealso{\code{\linkS4class{RIfeatures}},\code{\link{createPotentialRIs}},\code{\linkS4class{CorNeighbors}},\code{\link{RIP.universe}}
}
\examples{
data(RIP.demo_ExpressionSets)

## Calculate the Pearson correlation objects for each ExpressionSet:
test_eset.filenames <- createCorData(ExpressionSets=test_eset.names, entrez.annotation=NULL, series="RIP_testset", path=".")

## Retreive the correlation neighbors and the overall mean correlation within the datasets:
test_eset.CorNeighbors <- getCorNeighbors(cordata.filenames=test_eset.filenames, path=".", CC = 0.6, FoC = 0.25)

## Create some potential regulatory interactions (RIs) to be classified:
TFs <- sample(RIP.universe$TF, size=25)
targets <- sample(test_eset.CorNeighbors@entrez_annotation, size=100)
test_eset.RIs <- createPotentialRIs(target.genes=targets, TFs=TFs, annotation.targets="entrezID")

## Calculate the required features of the RIs for classification:
test_eset.features <- calculateRIPfeatures(pRIs=test_eset.RIs, cor.neighbors=test_eset.CorNeighbors)
}
