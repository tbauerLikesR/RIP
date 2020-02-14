setClass(
    Class="RIfeatures",
    representation(features="matrix", pRIs="data.frame", CorNeighbors="character", PWMscan="character"),
    prototype=list(features=matrix(0,10,0), pRIs=data.frame(TF=character(), gene=character()), CorNeighbors=character(), PWMscan=character())
)
