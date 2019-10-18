#' Get whitelist for high quality cells
#'
#' @param counts dgTMatrix of counts for all barcodes.
#' @param data_dir Directory that contains kallisto \code{'bus_output'} folder that counts were loaded from.
#'   whitelist is loaded from here if previously saved.
#'
#' @return Character vector of barcodes called as high quality cells.
#' @export
#' @keywords internal
get_scseq_whitelist <- function(counts, data_dir, overwrite = FALSE) {

  # check for previous whitelist
  whitelist_path <- file.path(data_dir, 'whitelist.txt')
  if (file.exists(whitelist_path) & !overwrite) {
    whitelist <- readLines(whitelist_path)
    return(whitelist)
  }

  # based on salmon alevin paper
  # get knee and keep at least 1000 below it
  knee <- get_knee(counts)
  ncount <- Matrix::colSums(counts)
  ncount.ord <- order(ncount, decreasing = TRUE)

  counts <- counts[, ncount.ord]
  ncount <- ncount[ncount.ord]
  nabove <- sum(ncount > knee)
  nbelow <- max(0.2*nabove, 1000)

  # if already filtered cellranger matrix file then keep everything
  if (min(ncount) > 10) {
    keep <- seq_len(ncol(counts))

  } else {
    keep <- 1:(nabove + nbelow)
  }

  # add qc metrics
  scseq <- Seurat::CreateSeuratObject(counts[, keep])
  sce <- srt_to_sce(scseq)
  sce <- add_scseq_qc_metrics(sce)

  # detect outliers based on PCA and mitochondrial content above knee
  sce_above <- scater::runPCA(sce[, 1:nabove], use_coldata = TRUE, detect_outliers = TRUE)
  mito.drop <- scater::isOutlier(sce_above$pct_counts_mito, nmads=3, type="higher")
  outl.drop <- sce_above$outlier

  # below knee is low quality
  df <- as.data.frame(sce@colData)


  # if filtered cellranger then just drop outliers, no modeling
  if (min(ncount) > 10) {
    df$quality <- 'high'
    df$quality[sce$total_counts <= knee] <- 'low'

    # set outliers above the knee as low quality
    df[colnames(sce_above)[mito.drop], 'quality'] <- 'low'
    df[colnames(sce_above)[outl.drop], 'quality'] <- 'low'

  } else {
    df$quality <- NA
    df$quality[sce$total_counts <= knee] <- 'low'

    # devide above evenly into high and ambiguous
    midpnt <- round(nabove/2)
    df$quality[1:midpnt] <- 'high'
    df$quality[(midpnt+1):nabove] <- 'ambig'

    # set outliers above the knee as low quality
    df[colnames(sce_above)[mito.drop], 'quality'] <- 'low'
    df[colnames(sce_above)[outl.drop], 'quality'] <- 'low'

    # clean df
    nunique <- function(x) { length(unique(na.omit(x))) != 1 }
    df <- df[, apply(df, 2, nunique)]
    df <- df[, !duplicated(t(df))]

    # run model/get preds
    df$quality <- factor(df$quality)
    train <- df[df$quality != 'ambig', ]
    test  <- df[df$quality == 'ambig', ]

    svm_model <- e1071::svm(quality ~ ., data = train)
    svm_preds <- predict(svm_model, newdata = test)

    # determine what to keep and save for future
    test$quality <- svm_preds
    df <- rbind(train, test)

  }

  whitelist <- row.names(df)[df$quality == 'high']
  kneelist  <- colnames(sce_above)

  writeLines(whitelist, whitelist_path)
  writeLines(kneelist, file.path(data_dir, 'kneelist.txt'))
  return(whitelist)
}


#' Get Roryk knee point and plot on barcode rank plot
#'
#' @param counts dgTMatrix of counts for all barcodes.
#'
#' @return count value corresponding to Roryk knee point
#' @export
#' @keywords internal
get_knee <- function(counts) {
  bcrank <- DropletUtils::barcodeRanks(counts)

  # Only showing unique points for plotting speed.
  uniq <- !duplicated(bcrank$rank)

  plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
       xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

  roryk_knee <- pick_roryk_cutoff(bcrank$total)

  abline(h=bcrank@metadata$inflection, col="darkgreen", lty=2)
  abline(h=bcrank@metadata$knee, col="dodgerblue", lty=2)
  abline(h=roryk_knee, col="red", lty=2)

  legend("bottomleft", legend=c("Inflection", "Knee", "Roryk Knee"),
         col=c("darkgreen", "dodgerblue", "red"), lty=2, cex=1.2)

  return(bcrank@metadata$inflection)
}

#' Pick Roryk knee point
#'
#' from https://github.com/COMBINE-lab/salmon/issues/362
#'
#' @param bcs Vector of total counts for all barcodes
#'
#' @return count value corresponding to Roryk knee point
#' @export
#' @keywords internal
pick_roryk_cutoff <- function(bcs){
  bcs_hist = hist(log10(bcs), plot=FALSE, n=50)
  mids = bcs_hist$mids
  vals = bcs_hist$count
  wdensity = vals * (10^mids) / sum(vals * (10^mids))
  baseline <- median(wdensity)

  # Find highest density in upper half of barcode distribution
  peak <- which(wdensity == max(wdensity[((length(wdensity)+1)/2):length(wdensity)]))

  # Cutoff is the point before the peak at which density falls below 2X baseline
  10^mids[max(which(wdensity[1:peak] < (2*baseline)))]
}
