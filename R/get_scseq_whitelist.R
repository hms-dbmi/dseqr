#' Get whitelist for high quality cells
#'
#' @param counts dgTMatrix of counts for all barcodes.
#' @inheritParams load_scseq
#'
#' @return Character vector of barcodes called as high quality cells.
#' @export
#' @keywords internal
get_scseq_whitelist <- function(counts, data_dir, overwrite = FALSE, knee_type = c('inflection', 'roryk', 'knee')) {

  # check for previous whitelist
  whitelist_path <- file.path(data_dir, 'whitelist.txt')
  if (file.exists(whitelist_path) & !overwrite) {
    whitelist <- readLines(whitelist_path)
    return(whitelist)
  }

  # based on salmon alevin paper
  # get knee and keep at least 1000 below it
  knee <- get_knee(counts, knee_type = knee_type)
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
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts[, keep]))
  sce <- add_scseq_qc_metrics(sce)

  # detect outliers based on PCA and mitochondrial content above knee
  df <- as.data.frame(sce@colData)
  stats <- df[1:nabove, c('sum', 'detected', 'subsets_mito_percent', 'subsets_ribo_percent')]
  stats$sum <- log10(stats$sum)
  stats$detected <- log10(stats$detected)

  outlying <- robustbase::adjOutlyingness(stats, only.outlyingness = TRUE)
  outl.drop <- scater::isOutlier(outlying, type = "higher")
  mito.drop <- scater::isOutlier(stats$subsets_mito_percent, type="higher")

  # below knee is low quality
  # if filtered cellranger then just drop outliers, no modeling
  if (min(ncount) > 10) {
    df$quality <- 'high'
    df$quality[sce$total <= knee] <- 'low'

    # set outliers above the knee as low quality
    df[row.names(stats)[mito.drop], 'quality'] <- 'low'
    df[row.names(stats)[outl.drop], 'quality'] <- 'low'

  } else {
    df$quality <- NA
    df$quality[sce$total <= knee] <- 'low'

    # devide above evenly into high and ambiguous
    midpnt <- round(nabove/2)
    df$quality[1:midpnt] <- 'high'
    df$quality[(midpnt+1):nabove] <- 'ambig'

    # set outliers above the knee as low quality
    df[row.names(stats)[mito.drop], 'quality'] <- 'low'
    df[row.names(stats)[outl.drop], 'quality'] <- 'low'

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
  kneelist  <- row.names(stats)

  writeLines(whitelist, whitelist_path)
  writeLines(kneelist, file.path(data_dir, 'kneelist.txt'))
  return(whitelist)
}


#' Get knee point and plot on barcode rank plot
#'
#' @inheritParams get_scseq_whitelist
#'
#' @return count value corresponding to knee point
#' @export
#' @keywords internal
get_knee <- function(counts, knee_type = c('inflection', 'roryk', 'knee')) {
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

  switch(knee_type[1],
         'inflection' = bcrank@metadata$inflection,
         'roryk' = roryk_knee,
         'knee' = bcrank@metadata$knee)
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
