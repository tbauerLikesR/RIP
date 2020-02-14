createCorData <-
function(ExpressionSets=NULL, entrez.annotation=NULL, series="", path=".", overwrite.default=FALSE) {
    old.path <- getwd()
    on.exit(try(setwd(old.path)))
    setwd(path)
    ExpressionSets <- unique(ExpressionSets)
    if (! all(sapply(X=ExpressionSets, FUN=function(x) is(get(x), "ExpressionSet"))))
        stop("Object 'ExpressionSets' must be a list of available ExpressionSet names, each with identical features (probes/genes) in the rows!")
    if (! all(sapply(X=seq_along(ExpressionSets), FUN=function(x) identical(featureNames(get(ExpressionSets[[x]])), featureNames(get(ExpressionSets[[1]]))))))
        stop("Some ExpressionSets listed in 'ExpressionSets' have divergent features (probes/genes)!")
    if (is.null(entrez.annotation)) {
        entrez.annotation <- featureNames(get(ExpressionSets[[1]])) # entrez gene annotation is mandatory
    } else {
        if(length(entrez.annotation) != nrow(get(ExpressionSets[[1]]))) stop("Provided 'entrez.annotation' is not as long as the number of features in the ExpressionSets!")
    }
    entrez.annotation <- as.character(entrez.annotation)
    if (any(duplicated(entrez.annotation[!is.na(entrez.annotation)])))
        warning("Duplicated Entrezgene IDs found in annotation. This should be avoided by pre-processing!")
    entrez.indices <- entrez.annotation %in% RIP.universe$entrez$ID
    if (sum(entrez.indices) == 0)
        stop("No common Entrez gene IDs found in RIP-data and your provided annotation. Entrez annotation must be Entrez Gene IDs, either in 'entrez.annotation' or in the featureNames of the ExpressionSets.\nSee object 'RIP.universe$entrez$ID' for available genes.")
    if (sum(entrez.indices) > 1000)
        cat("\nPlease be aware that you will need some Gb of RAM and disc space for saving the correlation data if your ExpressionSets contain thousands of genes!\n\n")
    cat("Initializing.\n")
    cor.indices <- upper.tri(matrix(nrow=sum(entrez.indices), ncol=sum(entrez.indices)))
    filenames.CorData <- NULL
    filenames.CorData <- paste(ExpressionSets, "_CorData.RData", sep="")
    filenames.written <- sapply(X=seq_along(ExpressionSets), FUN=function(x, expr_data=ExpressionSets, ind=cor.indices) {
        cat("Processing dataset ", x, " of ", length(expr_data), "...\n", sep="")
        correlation <- current.filename <- NULL
        if (file_test(op="-f", x=filenames.CorData[x])) {
            flag.overwrite <- userQuery(msg=paste("File '", filenames.CorData[x], "' already exists, overwrite?", sep=""), default=ifelse(overwrite.default, "y", "n"))
            if (flag.overwrite == "n") {
                cat("\tdata is expected to already be in '", filenames.CorData[x], "'.\n", sep="")
                return(filenames.CorData[x])
            }
        }
        cat("\tadjusting dataset to available gene list: ", sum(entrez.indices), " genes.\n", sep="")
        correlation <- cor(t(exprs(get(expr_data[x])[entrez.indices,])), method="pearson")[ind]
        cordata <- new("CorData", dataset=expr_data[x], series=series, cor_coef=correlation, entrez_annotation=entrez.annotation[entrez.indices])
        cat("\tsaving data to '", filenames.CorData[x], "'.\n", sep="")
        save(cordata, file=filenames.CorData[x])
        rm(cordata, correlation)
        invisible(gc())
        return(filenames.CorData[x])
    })
    rm(cor.indices)
    invisible(gc())
    cat("All files can be found in '", path, "'.\nReturning file names of saved correlation data.\n", sep="")
    return(filenames.CorData)
}

