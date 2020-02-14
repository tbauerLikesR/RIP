calculateRIPfeatures <-
function(pRIs=NULL, cor.neighbors=NULL, show.progressbar=TRUE) (UseMethod("calculateRIPfeatures"))

calculateRIPfeatures.default <-
function(pRIs=NULL, cor.neighbors=NULL, show.progressbar=TRUE) {
    features <- list()
    if (! all(is.data.frame(pRIs), c("TF", "gene") %in% names(pRIs)))
        stop("'pRIs' must be a 'data.frame' containing two columns named 'TF' and 'gene'!\n")
    stopifnot(class(cor.neighbors) == "CorNeighbors")
    if (nrow(pRIs) > 10000) cat("\nCalculating features for a large number (>10,000) of candidate RIs can be very time-consuming.\n==> Please consider running several instances simultaneously if this takes too long!\n")
    if (! all(pRIs$gene %in% cor.neighbors@entrez_annotation)) {
        warning("Only ", sum(pRIs$gene %in% cor.neighbors@entrez_annotation), " out of ", length(pRIs$gene), " RIs are in the CorNeighbors genes. The RIs were reduced to those genes!" , sep="")
        pRIs <- pRIs[pRIs$gene %in% cor.neighbors@entrez_annotation,]
    }
    if (! all(pRIs$TF %in% RIP.universe$TF))
        stop("Not all TFs are in the available TF pool (see 'RIP.universe$TF')!")
    if (! all(pRIs$gene %in% RIP.universe$entrez$ID))
        stop("Not all genes are in the available target gene pool (see 'RIP.universe$entrez$ID')!")
    if (! is.list(RIP.PWMscan))
        stop("You need to provide a PWM-scan list object containing . See 'RIP.PWMscan' for example!")
    if (! all(RIP.universe$core.genes %in% cor.neighbors@entrez_annotation))
        warning("Only ", sum(RIP.universe$core.genes %in% cor.neighbors@entrez_annotation), " out of ", length(RIP.universe$core.genes), " genes from the gold standard are included in the 'CorNeighbors' genes.\nThis may lead to biased results for some features!", sep="")
    on.exit(invisible(gc()))
    cat("Preparing look-up tables...\n")
    cor_neighbor_list <- unstack(as.data.frame(rbind(cor.neighbors@cor_neighbors, cor.neighbors@cor_neighbors[,2:1])))
    invisible(gc())
    .gs.md <- sapply(1:400, RIP.classifiers$mdi)
    cat("Calculating Feature 1...\n")
    features$"01" <- sapply(X=pRIs$gene, FUN=function(x) return(length(cor_neighbor_list[[x]])))
    cat("Calculating Feature 2...\n")
    .t <- match(pRIs$TF, RIP.universe$TF, nomatch=0)
    features$"02" <- sapply(X=seq_along(pRIs$TF), FUN=function(x) length(intersect(RIP.universe$entrez$ID[.gs.md[[.t[x]]]], cor_neighbor_list[[pRIs$gene[x]]])))
    cat("Calculating Feature 3...\n")
    total <- 100
    if (nrow(pRIs) >= total && show.progressbar) {
        pb <- txtProgressBar(min=0, max=total, style=3)
        progress <- 0
        step <- nrow(pRIs) %/% total
    }
    features$"03" <- sapply(X=seq_along(pRIs$TF), FUN=function(x) {
        if (nrow(pRIs) >= total && show.progressbar) {
            if ((x %% step) == 1) {
                progress <<- progress + 1
                setTxtProgressBar(pb, progress)
            }
        }
        gene <- pRIs$gene[x]
        is.cor.neighbor <- cor.neighbors@entrez_annotation %in% c(gene, cor_neighbor_list[[gene]])
        is.m <- cor.neighbors@entrez_annotation %in% RIP.universe$entrez$ID[.gs.md[[.t[x]]]]
        if ((sum(is.cor.neighbor) <= 1) | (! any(is.m))) return(0)
        contingency.table <- table(is.cor.neighbor, is.m)
        -log10(fisher.test(contingency.table, alternative="greater")$p.value)
    })
    if (nrow(pRIs) >= total && show.progressbar) {
        close(pb)
        rm(total, pb, progress, step)
    }
    cat("Calculating Feature 4...\n")
    features$"04" <- sapply(X=seq_along(pRIs$TF), FUN=function(x) {
        gene <- pRIs$gene[x]
        neighbors <- cor_neighbor_list[[gene]]
        if (length(neighbors) > 0) return(length(intersect(names(RIP.PWMscan[[pRIs$TF[x]]]), neighbors)))
        return(0)
    })
    cat("Calculating Feature 5...\n")
    total <- 100
    if (nrow(pRIs) >= total && show.progressbar) {
        pb <- txtProgressBar(min=0, max=total, style=3)
        progress <- 0
        step <- nrow(pRIs) %/% total
    }
    features$"05" <- sapply(X=1:nrow(pRIs), FUN=function(x) {
        if (nrow(pRIs) >= total && show.progressbar) {
            if ((x %% step) == 1) {
                progress <<- progress + 1
                setTxtProgressBar(pb, progress)
            }
        }
        TF <- pRIs$TF[x]
        gene <- pRIs$gene[x]
        is.cor.neighbor <- cor.neighbors@entrez_annotation %in% c(gene, cor_neighbor_list[[gene]])
        has.PWM <- cor.neighbors@entrez_annotation %in% names(RIP.PWMscan[[TF]])
        if (!(any(is.cor.neighbor) & any(has.PWM))) return(0)
        contingency.table <- table(is.cor.neighbor, has.PWM)
        -log10(fisher.test(contingency.table, alternative="greater")$p.value)
    })
    if (nrow(pRIs) >= total && show.progressbar) {
        close(pb)
        rm(total, pb, progress, step)
    }
    cat("Calculating Feature 6...\n")
    features$"06" <- sapply(X=seq_along(pRIs$TF), FUN=function(x) {
        corneighbors <- cor_neighbor_list[[pRIs$gene[x]]]
        if (length(corneighbors) > 0) return(sum((corneighbors %in% names(RIP.PWMscan[[pRIs$TF[x]]])) & (corneighbors %in% RIP.universe$entrez$ID[.gs.md[[.t[x]]]])))
        return(0)
    })
    cat("Calculating Feature 7...\n")
    features$"07" <- sapply(X=seq_along(pRIs$TF), FUN=function(x) {
        TF <- pRIs$TF[x]
        gene <- pRIs$gene[x]
        if (pRIs$gene[x] %in% names(RIP.PWMscan[[TF]])) return(-log10(RIP.PWMscan[[TF]][gene]))
        return(0)
    })
    cat("Calculating Feature 8...\n")
    features$"08" <- sapply(X=.gs.md, FUN=length)[.t]
    cat("Creating correlation indices for Features 9 and 10...\n")
    cor.index.matrix <- matrix(NA, nrow=length(cor.neighbors@entrez_annotation), ncol=length(cor.neighbors@entrez_annotation))
    for (i in 2:ncol(cor.index.matrix)) {
        indices <- (choose(n=i-1, k=2) + 1):(choose(n=i, k=2))
        cor.index.matrix[seq_along(indices),i] <- indices
    }
    rm(i)
    get.index <- function(x, i.mat=cor.index.matrix) {
        res <- c(i.mat[,x], i.mat[x,])
        res[!is.na(res)]
    }
    cat("Creating feature indices...\n")
    feature.indices.match <- match(pRIs$gene, cor.neighbors@entrez_annotation, nomatch=0)
    gs.indices.match <- sapply(X=.gs.md, FUN=function(x) {
        res <- match(RIP.universe$entrez$ID[x], cor.neighbors@entrez_annotation)
        res[! is.na(res)]
    })
    names(gs.indices.match) <- RIP.universe$TF
    gs.indices.match <- gs.indices.match[names(gs.indices.match) %in% pRIs$TF]
    if (any(sapply(X=gs.indices.match, FUN=length) == 0)) cat("\nSome modules do not have any genes in common with the genes of the correlation data. Feature 9 (average correlation to module) will be set to 0 for these.\n")
    cat("Calculating Features 9 and 10...\n")
    total <- 100
    if (nrow(pRIs) >= total && show.progressbar) {
        pb <- txtProgressBar(min=0, max=total, style=3)
        progress <- 0
        step <- nrow(pRIs) %/% total
    }
    features.9and10 <- lapply(X=seq_along(pRIs$TF), FUN=function(index) {
        if (nrow(pRIs) >= total && show.progressbar) {
            if ((index %% step) == 1) {
                progress <<- progress + 1
                setTxtProgressBar(pb, progress)
            }
        }

        TF <- pRIs$TF[index]
        gene <- feature.indices.match[index]
        gene.indices <- get.index(gene)
        m <- setdiff(gs.indices.match[[TF]], gene)
        ### Catch problem - no matching target gene (current TF-module is empty because all its genes are missing in the CorNeighbors):
        ### Set the mean correlation with TF-module to 0 (no correlation)
        if (length(gs.indices.match[[TF]]) == 0) {
            meancor.md <- 0
            meancor.nonmd <- mean(cor.neighbors@mean_cor[setdiff(gene.indices, as.vector(unlist(sapply(X=m, FUN=get.index))))])
            return(c(meancor.md, meancor.nonmd))
        }
        ### Catch problem - current gene is itself the only matching target gene (may occur if some RIP core genes are missing in CorNeighbors):
        ### Set the mean correlation with TF-module to 1 (correlates perfect with itself)
        if (setequal(gs.indices.match[[TF]], gene)) {
            meancor.md <- 1
            meancor.nonmd <- mean(cor.neighbors@mean_cor[setdiff(gene.indices, as.vector(unlist(sapply(X=m, FUN=get.index))))])
            return(c(meancor.md, meancor.nonmd))
        }
        meancor.md <- mean(cor.neighbors@mean_cor[intersect(gene.indices, as.vector(unlist(sapply(X=m, FUN=get.index))))])
        meancor.nonmd <- mean(cor.neighbors@mean_cor[setdiff(gene.indices, as.vector(unlist(sapply(X=m, FUN=get.index))))])
        return(c(meancor.md, meancor.nonmd))
    })
    if (nrow(pRIs) >= total && show.progressbar) {
        close(pb)
        rm(total, pb, progress, step)
    }
    features$"09" <- sapply(X=features.9and10, FUN=function(x) x[1])
    features$"10" <- sapply(X=features.9and10, FUN=function(x) x[2])
    features.matrix <- NULL
    for (index in seq_along(features)) {
        cf <- as.vector(features[[index]])
        names(cf) <- NULL
        features.matrix <- rbind(features.matrix, cf)
    }
    rownames(features.matrix) <- NULL
    rm(cor_neighbor_list, features, .gs.md, get.index, features.9and10, cor.index.matrix, feature.indices.match, gs.indices.match)
    NAs.found <- is.na(as.matrix(features.matrix))
    if (any(NAs.found)) {
        cat("\nWarning: Calculated RI-features contain NAs:\n")
        NA.rows <- unique(row(NAs.found)[NAs.found])
        NA.rowcounts <- rowSums(NAs.found)[NA.rows]
        cat("Feature:\t", paste(NA.rows, collapse="\t"), "\n", sep="")
        cat("NA-count:\t", paste(NA.rowcounts, collapse="\t"), "\n\n", sep="")
    }
    return(new("RIfeatures", features=as.matrix(features.matrix), pRIs=pRIs, CorNeighbors=deparse(substitute(cor.neighbors)), PWMscan="RIP.PWMscan"))
}
