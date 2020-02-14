\name{getCorNeighbors}
\alias{getCorNeighbors}
\title{Processing CorData objects to select RIP correlation neighbors and their mean correlation
}
\description{
CorData objects are loaded from file and processed to select RIP correlation neighbors by given parameters CC and FoC. Additionally, their mean correlation is stored.
}
\usage{
getCorNeighbors(cordata.filenames=NULL, path=".", CC=0.6, FoC=0.25)
}
\arguments{
  \item{cordata.filenames}{A character vector of file names containing objects of class \code{CorData}.}
  \item{path}{character. The path where the files are located.}
  \item{CC}{numeric. A value between 0 and 1 defining the CC parameter. see details}
  \item{FoC}{numeric. A value between 0 and 1 defining the FoC parameter. see details}
}
\details{Each file containing a CorData object is loaded and processed by two parameters: CC and FoC.

CC defines the absolute value of the Pearson correlation threshold.

FoC defines the Fraction of Conditions threshold. The number of conditions is equal to the number of datasets that were used to create the CorData objects, where each dataset is expected to represent one class, condition or experimental setup (e.g. tissue/tumor type or treatment).

CC and FoC are applied so that gene pairs are selected if they have a correlation coefficient of >|CC| in more than 'FoC * number of datasets' conditions.

The standard parameters are set to CC = 0.6 and FoC = 0.25, so that in an analysis of 10 datasets, a gene pair is define a correlation neighbor if the absolute value of the Pearson correlation is higher than 0.6 in more than 2.5 datasets (i.e. at least 3).
}
\value{an object of class \code{CorNeighbors} which is a \code{list} of five arguments:
  \item{datasets}{A character vector containing the name(s) of the initial ExpressionSet(s) used to create the CorData object(s).}
  \item{parameters}{A numeric vector with the defined CC and FoC parameters.}
  \item{entrez_annotation}{A character vector with the provided Entrez gene ID annotation.}
  \item{cor_neighbors}{A matrix listing the correlation neighbors (Entrez gene IDs) - one per row.}
  \item{mean_cor}{A numeric vector with the mean correlation coefficient of all gene pairs from the annotation. The order is equvalent to the \code{upper.tri} of the correlation matrix.}
}
\references{
Bauer et al. (2011) Bioinformatics 27:2239-47 - RIP: the regulatory interaction predictor--a machine learning-based approach for predicting target genes of transcription factors.
}
\author{Tobias Bauer \email{t.bauer@dkfz.de}}
\seealso{\code{\linkS4class{CorNeighbors}}
}
\examples{
data(RIP.demo_ExpressionSets)

## Calculate the Pearson correlation objects for each ExpressionSet:
test_eset.filenames <- createCorData(ExpressionSets=test_eset.names, entrez.annotation=NULL, series="RIP_testset", path=".")

## Retreive the correlation neighbors and the overall mean correlation within the datasets:
test_eset.CorNeighbors <- getCorNeighbors(cordata.filenames=test_eset.filenames, path=".", CC = 0.6, FoC = 0.25)
}
