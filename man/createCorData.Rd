\name{createCorData}
\alias{createCorData}
\title{Calculating the Pearson correlation coefficient for all pairs of features of ExpressionSets
}
\description{
Calculates the Pearson correlation coefficient for all pairs of features of given ExpressionSets. Data is adjusted to Entrez gene annotation available in the RIP classifier.
}
\usage{
createCorData(ExpressionSets=NULL, entrez.annotation=NULL, series="", path=".", overwrite.default=FALSE)
}
\arguments{
  \item{ExpressionSets}{a character vector containing one or more objects of class \code{ExpressionSet} available in the workspace. see details}
  \item{entrez.annotation}{a character vector providing the Entrez gene ID annotation for the features in \code{ExpressionSets}. Must be valid for all ExpressionSets. see details}
  \item{series}{character. An optional name for the dataset series.}
  \item{path}{character. The path files should be written to (if so).}
  \item{overwrite.default}{logical. The default used for overwriting an existing file in when running R in batch mode. If FALSE, no file will be overwritten, if TRUE, all will be.}
}
\details{
Currently, all ExpressionSet objects named in \code{ExpressionSets} must have the same form in terms of number and order of features.

If the \code{featureNames} of the ExpressionSets are Entrez gene IDs, \code{entrez.annotation} can be set to NULL. Otherwise, you need to provide the Entrez gene ID for each feature there. The annotation must be valid for all ExpressionSets.

The data is checked to match the genes available for RIP as provided in \code{RIP.universe}.

All ExpressionSets are processed and saved to an individual file in the given \code{path} assembled from the ExpressionSet name and the ending '.RData'. This is done because most machines do not have the RAM capacity (64Gb up) to keep a series of ExpressionSet with thousands of genes in the active memory and process them. Note that the number of gene pairs for n genes equals n * (n-1) / 2. In the future I might write another function to skip the disc-writing (and follow-up loading) for high-end machines, but implementation has not been done so far.
}
\value{a character vector with the file names of each processed ExpressionSet. Each file contains an object of class \code{CorData}
}
\author{Tobias Bauer \email{t.bauer@dkfz.de}}
\seealso{\code{\linkS4class{CorData}}
}
\examples{
data(RIP.demo_ExpressionSets)

## Calculate the Pearson correlation objects for each ExpressionSet:
test_eset.filenames <- createCorData(ExpressionSets=test_eset.names, entrez.annotation=NULL, series="RIP_testset", path=".")
}
