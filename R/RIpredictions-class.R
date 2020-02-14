setClass(
    Class="RIpredictions",
    representation(predictions="data.frame", percent.votes="numeric", CorNeighbors="character", PWMscan="character"),
    prototype=list(predictions=data.frame(TF=character(), gene=character(), symbol=character(), percent.votes=numeric(), confidence=numeric()), CorNeighbors=character(), PWMscan=character())
)
