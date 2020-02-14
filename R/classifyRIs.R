classifyRIs <-
function(RIfeatures=NULL, percent.votes=80, filename=NULL, row.names=FALSE, show.progressbar=TRUE, classifier.name=NULL, ...) {
    if (! is.null(classifier.name)) {
        if (! exists(classifier.name)) {
            warning("Given classifier.name not found in workspace, using classifier 'RIP.classifiers' provided by the RIP-package instead.")
        } else {
            RIP.classifiers <- get(classifier.name)
        }
    }
    if (any(is.na(RIfeatures@features)))
        stop("Features must not contain NAs. This could be caused by missing 'core' genes in the 'CorNeighbors' during feature calculation.")
    percent.votes <- round(percent.votes)
    if (! all((percent.votes >= 0), (percent.votes <= 100)))
        stop("'percent.votes' must be between 0 and 100!")
    on.exit(invisible(gc()))
    pRIs <- RIfeatures@pRIs
    pRIs$symbol <- RIP.universe$entrez$symbol[match(as.character(pRIs$gene), RIP.universe$entrez$ID)]
    classifier.indices <- as.numeric(which(sapply(RIP.classifiers, function(x) "classifier.for.matrix" %in% names(x))))
    cat("Applying RIP classifier to ", ncol(RIfeatures@features), " features.\n", sep="")
    total <- 100
    if (length(classifier.indices) >= total && show.progressbar) {
        pb <- txtProgressBar(min=0, max=total, style=3)
        progress <- 0
        step <- length(classifier.indices) %/% total
    }
    total.votes <- sapply(X=seq_along(classifier.indices), FUN=function(x) {
        if (length(classifier.indices) >= total && show.progressbar) {
            if ((x %% step) == 0) {
                progress <<- progress + 1
                setTxtProgressBar(pb, progress)
            }
        }
        return(as.numeric(RIP.classifiers[[x]]$classifier.for.matrix(t(as.matrix(RIfeatures@features)[RIP.classifiers$feature.order,]))) - 1)
    })
    total.votes <- rowSums(total.votes)
    pRIs$percent.votes <- total.votes / length(classifier.indices) * 100
    if (length(classifier.indices) >= total && show.progressbar)
        close(pb)
    confidence <- pRIs$percent.votes
    confidence[confidence > 0 ] <- RIP.classifiers$precision[round(pRIs$percent.votes)] * 100
    pRIs$confidence <- round(confidence, digits=1)
    pRIs <- pRIs[order(pRIs$symbol),]
    pRIs <- pRIs[order(pRIs$TF),]
    pRIs <- pRIs[order(pRIs$confidence, decreasing=TRUE),]
    pRIs <- pRIs[pRIs$percent.votes >= percent.votes,]
    result <- new("RIpredictions", predictions=pRIs, percent.votes=percent.votes, CorNeighbors=RIfeatures@CorNeighbors, PWMscan=RIfeatures@PWMscan)
    if (! is.null(filename))
        try(write.csv(pRIs, file=filename, row.names=row.names, ...))
    return(invisible(result))
}

