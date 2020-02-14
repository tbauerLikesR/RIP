createPotentialRIs <-
function(target.genes=NULL, TFs=NULL, annotation.targets=c("entrezID", "symbol")) {
    stopifnot(all(c(length(target.genes), length(TFs)) > 0))
    target.genes <- as.character(target.genes)
    annotation.targets <- match.arg(arg=annotation.targets, choices=c("entrezID", "symbol"), several.ok=FALSE)
    if (annotation.targets == "symbol") {
        if (! any(target.genes %in% RIP.universe$entrez$symbol)) stop("No target gene symbol matches the available gene pool (see 'RIP.universe$entrez')!")
        target.genes <- RIP.universe$entrez$ID[match(target.genes, RIP.universe$entrez$symbol, nomatch=0)]
    }
    if(! any(target.genes %in% RIP.universe$entrez$ID)) stop("No target gene ID matches the available gene pool (see 'RIP.universe$entrez')!")
    if (any(duplicated(target.genes))) {
        target.genes <- unique(target.genes)
        cat("Deleted duplicated target gene entries (", length(target.genes), " genes kept).\n", sep="")
    }
    if (! any(TFs %in% RIP.universe$TF)) stop("No TF name matches the available TF pool (see 'RIP.universe$TF')!")
    if (any(duplicated(TFs))) {
        TFs <- unique(TFs)
        cat("Deleted duplicated TF entries (", length(TFs), " TFs kept).\n", sep="")
    }
    pRIs <- data.frame(cbind(
        TF=rep(TFs, each=length(target.genes)),
        gene=rep(target.genes, times=length(TFs))
    ), stringsAsFactors=FALSE)
    cat("Returning ", nrow(pRIs), " potential RIs from ", length(TFs), " TFs and ", length(target.genes), " target genes.\n")
    return(pRIs)
}

