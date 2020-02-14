setClass(
    Class="CorData",
    representation(dataset="character", series="character", cor_coef="numeric", entrez_annotation="character"),
    prototype=list(name=character(), datasets=character(), cor_coef=numeric(), entrez_annotation=character())
)
