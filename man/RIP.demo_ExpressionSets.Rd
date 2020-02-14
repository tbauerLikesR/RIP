\name{RIP.demo_ExpressionSets}
\alias{RIP.demo_ExpressionSets}
\alias{test_eset.names}
\alias{test_eset_01}
\alias{test_eset_02}
\alias{test_eset_03}
\alias{test_eset_04}
\alias{test_eset_05}
\alias{test_eset_06}
\alias{test_eset_07}
\alias{test_eset_08}
\alias{test_eset_09}
\alias{test_eset_10}
\docType{data}
\title{
Provided ExpressionSets for demonstation of RIP
}
\description{
This data provides 10 datasets for demonstation of RIP-functions. Each ExpressionSet resembles Affymetrix hgu133a-based gene expression data truncated to 1,000 genes after rma-normalization.
}
\usage{data(RIP.demo_ExpressionSets)}
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
\keyword{datasets}
