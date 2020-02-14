\name{classifyRIs}
\alias{classifyRIs}
\title{Classify potential RIs by RIP features
}
\description{
Applies the RIP classifier to potential regulatory interaction (RIs) between transcription factors (TFs) and target genes by ten previously calculated features.
}
\usage{
classifyRIs(RIfeatures=NULL, percent.votes=80, filename=NULL, row.names=FALSE, show.progressbar=TRUE, classifier.name=NULL, ...)
}
\arguments{
  \item{RIfeatures}{an object of class \code{RIfeatures}.}
  \item{percent.votes}{a numeric defining the fraction of votes (in percent) required to assign an RI to a TF - target gene pair.}
  \item{filename}{character. An optional file name for csv-export of the results.}
  \item{row.names}{logical. Should row-names be written to csv ouput file? (Row-names correspond to RI indices of input.)}
  \item{show.progressbar}{logical. Show a text progress-bar?}
  \item{classifier.name}{character. Advanced option - only use if you have a trained SVM master-classifier other than that provided by the RIP-package! (Provide variable name then.)}
  \item{\dots}{further optional arguments passed to the \code{write.csv} function, if \code{filename} is provided.}
}
\value{An object of class \code{RIpredictions} which is a \code{list} of four arguments:
  \item{predictions}{A data.frame with the resulting RIs: \code{TF} is the TF symbol, \code{gene} and \code{symbol} are the target gene Entrez ID and symbol respectively, \code{percent.votes} is the percentage of positive votes from the RIP classifier, and \code{confidence} is the estimated true positive rate taken from the RIP paper (see references).}
  \item{percent.votes}{A numeric with the user-defined vote cutoff in percent.}
  \item{CorNeighbors}{The name of the CorNeighbors input used for feature calculation.}
  \item{PWMscan}{The name of the PWMscan input used for feature calculation.}
}
\references{
Bauer et al. (2011) Bioinformatics 27:2239-47 - RIP: the regulatory interaction predictor--a machine learning-based approach for predicting target genes of transcription factors.
}
\author{Tobias Bauer \email{t.bauer@dkfz.de}}
\seealso{\code{\link{RIfeatures-class}},\code{\link{RIpredictions-class}}
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

## Classify the RIs and write results to file:
test_eset.classifedRIs <- classifyRIs(RIfeatures=test_eset.features, percent.votes=80, filename="test_classification_votes80.csv")

## Re-classify the RIs with a different vote cutoff and write results to file:
test_eset.classifedRIs <- classifyRIs(RIfeatures=test_eset.features, percent.votes=60, filename="test_classification_votes60.csv")
}
