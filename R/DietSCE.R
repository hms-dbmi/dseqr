#' An lighter-weight SingleCellExperiment S4 class
#'
#' @slot metadata list.
#' @slot colData data.frame.
#' @slot reducedDims list.
#'
#' @return
#' @export
#'
setClass("DietSCE",
         representation(
             assays = "list",
             metadata = "list",
             colData = "data.frame",
             rowData = "data.frame",
             reducedDims = "list")
)

DietSCE <- function(assays, metadata, colData, rowData, reducedDims) {
    new("DietSCE",
        assays = assays,
        metadata = metadata,
        colData = colData,
        rowData = rowData,
        reducedDims = reducedDims)
}


setGeneric("reducedDimNames", function(x) standardGeneric("reducedDimNames"))
setMethod("reducedDimNames", "DietSCE", function(x) names(x@reducedDims))

setGeneric("reducedDim", function(x, type) standardGeneric("reducedDim"))
setMethod("reducedDim", "DietSCE", function(x, type) x@reducedDims[[type]])

setGeneric("reducedDim<-", function(x, type, value) standardGeneric("reducedDim<-"))
setMethod("reducedDim<-", "DietSCE", function(x, type, value) x@reducedDims[[type]] <- value)

setGeneric("colData", function(x) standardGeneric("colData"))
setMethod("colData", "DietSCE", function(x) x@colData)

setGeneric("rowData", function(x) standardGeneric("rowData"))
setMethod("rowData", "DietSCE", function(x) x@rowData)

setGeneric("rowData<-", function(x, value) standardGeneric("rowData<-"))
setMethod("rowData<-", "DietSCE", function(x, value) x@rowData <- value)

setGeneric("assay", function(x, type) standardGeneric("assay"))
setMethod("assay", "DietSCE", function(x, type) x@assays[[type]])

setGeneric("assay<-", function(x, type, value) standardGeneric("assay<-"))
setMethod("assay<-", "DietSCE", function(x, type, value) x@assays[[type]] <- value)

setGeneric("assayNames", function(x, type) standardGeneric("assayNames"))
setMethod("assayNames", "DietSCE", function(x, type) names(x@assays))

setGeneric("assayNames<-", function(x, type, value) standardGeneric("assayNames<-"))
setMethod("assayNames<-", "DietSCE", function(x, type, value) names(x@assays) <- value)

setGeneric("reducedDim<-", function(x, type, value) standardGeneric("reducedDim<-"))
setMethod("reducedDim<-", "DietSCE", function(x, type, value) x@reducedDims[[type]] <- value)

setMethod("dimnames", "DietSCE", function(x) {
    list(row.names(x@rowData), row.names(x@colData))
})

setMethod("row.names", "DietSCE", function(x) row.names(x@rowData))
setMethod("colnames", "DietSCE", function(x) row.names(x@colData))
setMethod("ncol", "DietSCE", function(x) nrow(x@colData))
setMethod("nrow", "DietSCE", function(x) nrow(x@rowData))

# $ gets colData
setMethod("$", "DietSCE", function(x, name) {
    eval(substitute(x@colData$NAME_ARG, list(NAME_ARG=name)))
})

# so that get hints
.DollarNames.DietSCE <- function(x, pattern = "") {
    grep(pattern, names(x@colData), value=TRUE)
}


# assign to colData
setReplaceMethod("$", "DietSCE", function(x, name, value) {
    x@colData[[name]] = value
    x
})

setMethod("[", "DietSCE", function(x, i, j, ..., drop = FALSE) {
    if (missing(drop))
        drop <- FALSE
    if (missing(i) && missing(j)) {
        if (!missing(...))
            stop("specify genes or samples to subset; use '",
                 substitute(x), "$", names(list(...))[[1]],
                 "' to access colData variables")
        return(x)
    }
    if (!missing(j))
        x@colData <- x@colData[j,, ..., drop = drop]

    if (!missing(i))
        x@rowData <- x@rowData[i,,..., drop=drop]
    x
})


setMethod("[[", "DietSCE", function(x, i, j, ...) x@colData[[i]])

setReplaceMethod("[[", "DietSCE",
                 function(x, i, j, ..., value) {
                     x@colData[[i, ...]] <- value
                     x
                 })

setGeneric("logcounts<-", function(x, value) standardGeneric("logcounts<-"))
setMethod("logcounts<-", "DietSCE", function(x, value) {
    x@assays$logcounts <- value
    return(x)
})

setGeneric("logcounts", function(x) standardGeneric("logcounts"))
setMethod("logcounts", "DietSCE", function(x) x@assays$logcounts)

setGeneric("counts<-", function(x, value) standardGeneric("counts<-"))
setMethod("counts<-", "DietSCE", function(x, value) {
    x@assays$counts <- value
    return(x)
})

setGeneric("counts", function(x) standardGeneric("counts"))
setMethod("counts", "DietSCE", function(x) x@assays$counts)

# blah2 <- new('DietSCE',
#              assays = list(),
#              metadata = metadata(blah),
#              colData = as.data.frame(colData(blah)),
#              rowData = as.data.frame(rowData(blah)),
#              reducedDims = as.list(reducedDims(blah)))
