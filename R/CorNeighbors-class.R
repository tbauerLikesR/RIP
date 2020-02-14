setClass(
    Class="CorNeighbors",
    representation(datasets="character", parameters="numeric", entrez_annotation="character", cor_neighbors="matrix", mean_cor="numeric"),
    prototype=list(datasets=character(), parameters=c(CC=0, FoC=0), entrez_annotation=character(), cor_neighbors=matrix("",0,0), mean_cor=numeric())
)
