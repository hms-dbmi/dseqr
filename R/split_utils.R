split_mtx <- function(mtx, samples) {

    mtxs <- list()
    for (sample in unique(samples)) {
        is.sample <- samples == sample
        mtxs[[sample]] <-  mtx[, is.sample]
    }

    return(mtxs)
}


save_split_mtxs <- function(mtxs, features, data_dir) {
    colnames(features)[1:2] <- c('enid', 'symbol')
    fnames <- c('barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz')

    unique.genes <- make.unique(features$symbol)

    out_dir <- file.path(data_dir, 'split')
    dir.create(out_dir)

    for (i in seq_along(mtxs)) {
        samplei <- names(mtxs)[i]
        mtxi <- mtxs[[i]]

        if (!identical(unique.genes, row.names(mtxi)))
            stop('Gene names not same')

        sample_dir <- file.path(out_dir, samplei)

        DropletUtils::write10xCounts(sample_dir,
                                     mtxi,
                                     gene.id = features$enid,
                                     gene.symbol = features$symbol,
                                     version = '3')
        file.rename(file.path(sample_dir, fnames),
                    file.path(out_dir, paste0(samplei, '_', fnames)))

        unlink(sample_dir, recursive = TRUE)
    }
}
