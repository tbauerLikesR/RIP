getCorNeighbors <-
function(cordata.filenames=NULL, path=".", CC=0.6, FoC=0.25) {
    if (! all((0 <= c(CC, FoC)) & (c(CC, FoC) <= 1)))
        stop("Both 'CC' and 'FoC' must be between 0 and 1!")
    if (! all(cordata.filenames %in% dir(path)))
        stop("Could not find all files in given path.")
    old.path <- getwd()
    on.exit(try(setwd(old.path)))
    setwd(path)
    cor.list <- list()
    entrez.annotation <- entrez.indices <- mean.cor <- NULL
    cat("\nCombining all correlation data:\n")
    datasets <- NULL
    for (index in seq_along(cordata.filenames)) {
        cat("Loading file ", index, " of ", length(cordata.filenames), ".\n", sep="")
        load(cordata.filenames[index])
        if (exists("cordata")) {
            if (class(cordata) != "CorData") {
                cat("Could not find a CorData object in a varable named 'cordata' for file '", cordata.filenames[index], "'! Processing next filename.\n", sep="")
                next
            }
        } else {
            cat("Could not find a CorData object in a varable named 'cordata' for file '", cordata.filenames[index], "'! Processing next filename.\n", sep="")
            next
        }
        if (is.null(entrez.annotation)) {
            entrez.annotation <- cordata@entrez_annotation
            entrez.indices <- cbind(feature1.index=unlist(sapply(X=2:length(entrez.annotation), FUN=function(x) return(rep(x, each=x-1)))), feature2.index=unlist(sapply(X=1:(length(entrez.annotation) - 1), FUN=function(x) return(1:x))))
        }
        if (! length(cordata@cor_coef) == choose(length(entrez.annotation), 2)) {
            cat("Problem with data in file '", cordata.filenames[index], "':\nCorrelation data does either not fit the provided annotation or is not of the right size! Processing next filename.\n", sep="")
            next
        }
        cat("\tprocessing correlation coefficients...\n", sep="")
        if (is.null(mean.cor)) {
            mean.cor <- cordata@cor_coef / length(cordata.filenames)
        } else {
            mean.cor <- mean.cor + cordata@cor_coef
        }
        cor.coefs <- abs(cordata@cor_coef)
        cor.list[[index]] <- which(cor.coefs > CC)
        datasets <- c(datasets, cordata@dataset)
        rm(cordata, cor.coefs)
        invisible(gc())
    }
    dataset.count <- length(datasets)
    if (dataset.count != length(cordata.filenames))
        warning("Number of processed datasets is different from number of filenames. Correlation means were adjusted to ", dataset.count, " datasets accordingly!")
    cat("\nProcessed", dataset.count, "datasets.\nProcessing correlation indices...\n")
    cor.neighbor.candidates <- cor.neighbors.indices <- cor.neighbors <- NULL
    cor.neighbor.candidates <- table(unlist(cor.list))
    cor.neighbors.indices <- as.numeric(names(cor.neighbor.candidates[((cor.neighbor.candidates / dataset.count) > FoC)]))
    rm(cor.list, cor.neighbor.candidates)
    invisible(gc())
    cat("Resolving ", length(cor.neighbors.indices), " indices...\n", sep="")
    cor.neighbors <- cbind(entrez.annotation[entrez.indices[cor.neighbors.indices,1]], entrez.annotation[entrez.indices[cor.neighbors.indices,2]])
    mean.cor <- mean.cor / dataset.count
    result <- new("CorNeighbors", datasets=datasets, parameters=c(CC=CC, FoC=FoC), entrez_annotation=entrez.annotation, cor_neighbors=cor.neighbors, mean_cor=mean.cor)
    rm(cor.neighbors, mean.cor)
    invisible(gc())
    return(result)
}

