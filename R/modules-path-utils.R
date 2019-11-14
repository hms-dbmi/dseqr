#' Get perturbation signature
#'
#' @param pert Name of perturbation signature.
#' @param pert_type One of \code{'cmap'}, \code{'l1000_genes'}, or \code{'l1000_drugs'}.
#' @export
#' @keywords internal
load_pert_signature <- function(pert, pert_type, pert_signature_dir) {
  type_dir <- file.path(pert_signature_dir, pert_type)
  if (!file.exists(type_dir)) dir.create(type_dir)


  fname <- paste0(pert, '.rds')
  sig_path <- file.path(type_dir, fs::path_sanitize(fname))
  if (!file.exists(sig_path)) dl_pert_signature(sig_path, pert_type)

  sig <- readRDS(sig_path)
  return(sig)
}

#' Download CMAP02/L1000 pert signature from S3
#'
#' @param sig_path Path to download file to.
#' @param pert_type One of \code{'cmap'}, \code{'l1000_drugs'}, or \code{'l1000_genes'}.
#'
#' @return NULL
#' @export
#' @examples
#' sig_path <- 'cmap_res_BRD-K45319408_PC3_5um_24h.rds'
#' dl_pert_result(res_path)
#'
dl_pert_signature <- function(sig_path, pert_type) {
  # name of the file being requested
  dl_url <- paste0('https://s3.us-east-2.amazonaws.com/drugseqr/drug_es_dir/', pert_type, '/', basename(sig_path))
  dl_url <- utils::URLencode(dl_url)
  download.file(dl_url, sig_path, )
}

#' Used by get_path_df to construct the return result
#'
#' @param top_table Filtered result of \code{\link[limma]{topTable}}
#' to limit number of plotted genes.
#' @export
#' @keywords internal
construct_path_df <- function(top_table) {

  data.frame(
    Gene = row.names(top_table),
    Dprime = top_table$dprime,
    sd = sqrt(top_table$vardprime),
    description = tx2gene$description[match(row.names(top_table), tx2gene$gene_name)],
    Link = paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", row.names(top_table), "'>", row.names(top_table), "</a>"),
    stringsAsFactors = FALSE
  ) %>%
    mutate(Gene = factor(Gene, levels = Gene))
}


#' Get data.frame for plotting gene expression values of a pathway
#'
#' @param anal Result of call to \code{\link{diff_expr_scseq}}.
#' @param path_id String with KEGG pathway id.
#' @param path_genes Character vector of custom genes to construct pathway data.frame for.
#' @param nmax Maximum number of genes to keep from CMAP02/L1000 common and CMAP02 only genes for Drug and genetic query genes. Default is 200
#'  so that all drug and genetic query genes are shown.
#'
#' @return \code{data.frame} with columns: \itemize{
#'  \item Gene gene names.
#'  \item Dprime standardized unbiased effect size values.
#'  \item sd standard deviations of \code{Dprime}.
#'  \item Link url to GeneCards page for gene.
#' }
#' @export
#' @keywords internal
get_path_df <- function(anal, path_id = NULL, pert_signature = NULL, nmax = 200) {

  # add dprimes and vardprime values
  anal <- add_es(anal)
  top_table <- anal$top_table
  top_table <- top_table[order(abs(top_table$dprime), decreasing = TRUE), ]

  # only show pathway if in kegg
  if (path_id %in% names(gslist.kegg)) {
    path_enids <- gslist.kegg[[path_id]]
    path_genes <- names(path_enids)
    top_table <- top_table[row.names(top_table) %in% path_genes, ]
  }

  path_df <- construct_path_df(top_table)
  path_df$color = 'black'

  # for drug/genetic queries: keep up to nmax in cmap/l1000
  is.cmap2 <- path_id == 'Query genes - CMAP02'
  is.l1000 <- path_id == 'Query genes - L1000'

  if (is.l1000 | is.cmap2) {
    if (is.cmap2) keep <- head(which(path_df$Gene %in% unlist(genes)), nmax)
    if (is.l1000) keep <- head(which(path_df$Gene %in% genes$common), nmax)

    path_df <- path_df[row.names(path_df) %in% keep, ]
  }

  # add pert signature data
  if (!is.null(pert_signature)) {
    pert_signature <- pert_signature[as.character(path_df$Gene)]
    pert_df <- path_df
    pert_df$Dprime <- pert_signature

    # get colors based on signs
    path_signs <- sign(path_df$Dprime)
    pert_signs  <- sign(pert_df$Dprime)
    pert_df$color <- ifelse(pert_signs == path_signs, 'red', 'blue')

    # TODO show pert variances
    pert_df$sd <- NA
    path_df <- rbind(path_df, pert_df)

  }

  return(path_df)
}

#' Generate data.frame of saved single cell RNA-Seq analyses
#'
#' Used together with \code{\link{load_bulk_anals}} to create data.frame
#' for selectizeInput choices.
#'
#' @param data_dir Directory to folder with single-cell analysis folders.
#' @param with_type if \code{TRUE} include column \code{'type'} with values \code{'Single Cell'}. Used for optgroupField in
#'  selectizeInput's. Default is \code{FALSE}
#'
#' @return data.frame with columns: \itemize{
#'  \item dataset_name \code{NA} included for consistency with \code{load_bulk_anals}.
#'  \item dataset_dir Directory name in \code{data_dir} with single cell analysis.
#'  \item anal_name Name of single cell analysis. Same as \code{dataset_dir}.
#'
#' }
#' @export
#' @keywords internal
load_scseq_anals <- function(data_dir, with_type = FALSE) {
  int_path <- file.path(data_dir, 'single-cell', 'integrated.rds')

  anals <- data.frame(matrix(ncol = 3, nrow = 0), stringsAsFactors = FALSE)
  colnames(anals) <- c("dataset_name", "dataset_dir", "anal_name")

  if (file.exists(int_path)) {
    integrated <- readRDS(int_path)
    for(anal in integrated) anals[nrow(anals) + 1, ] <- c(NA, anal, anal)

  }

  anals$label <- anals$anal_name
  anals$value <- seq_len(nrow(anals))

  if (with_type & nrow(anals)) anals$type <- 'Single Cell'

  return(anals)
}

#' Run PADOG pathway analysis on a single cell RNA-Seq dataset
#'
#' Analysis is performed comparing test to control cells for all selected clusters.
#'
#' @param scseq \code{Seurat} object subsetted to selected clusters
#' @param prev_anal Result of call to \code{diff_expr_scseq}
#' @param data_dir Directory to folders with single cell analyses
#' @param anal_name Name of folder in \code{data_dir} to save results to
#' @param clusters_name String with sorted comma seperated integers of selected clusters returned from \code{collapse_sorted}.
#'
#' @return result of \code{\link[PADOG]{padog}}
#' @export
#' @keywords internal
diff_path_scseq <- function(scseq, prev_anal, ambient, data_dir, anal_name, clusters_name, NI = 1000) {
  assay <- get_scseq_assay(scseq)
  Seurat::DefaultAssay(scseq) <- assay

  # load previous if exists
  fname <- paste0('diff_path_kegg_', clusters_name, '.rds')
  fpath <- file.path(data_dir, anal_name, fname)

  # for compatibility with previous versions
  fname_old <- paste0('diff_path_', clusters_name, '.rds')
  fpath_old <- file.path(data_dir, anal_name, fname_old)

  if(file.exists(fpath_old)) file.rename(fpath_old, fpath)
  if(file.exists(fpath)) return(readRDS(fpath))

  # subset to non-ambient geness
  scseq <-  scseq[!row.names(scseq) %in% ambient, ]

  # get groups
  group <- as.character(scseq$orig.ident)
  group <- ifelse(group == 'ctrl', 'c', 'd')

  # expression matrix
  # SCT corrected log counts
  esetm  <- scseq[[assay]]@data

  # already annotated with hgnc symbols
  gslist.kegg <- lapply(gslist.kegg, function(gs) {ret <- names(gs); names(ret) <- ret; return(ret)})

  # run padog
  padog_table <- PADOG::padog(esetm = esetm, group = group, parallel = TRUE, ncr = 4, gs.names = gs.names.kegg, gslist = gslist.kegg,
                              verbose = FALSE, rna_seq = FALSE, NI = NI)

  # save results
  saveRDS(padog_table, fpath)

  return(padog_table)
}

#' Run limma differential expression between test and control groups in single cell RNA-Seq dataset
#'
#' @param scseq \code{Seurat} object
#' @param clusters Character vector of clusters to include in analysis
#' @seealso \code{\link{get_ambient}}
#'
#' @return Named list with slots: \itemize{
#'  \item top_table result of limma::topTable.
#'  \item ebayes_sv results of \code{\link{fit_ebayes_scseq}}.
#'  \item pdata \code{data.frame} with column \code{'group'} indicating test and control samples.
#' }
#' @export
#' @keywords internal
diff_expr_scseq <- function(scseq, data_dir, anal_name, clusters_name) {

  # pseudo bulk analysis
  has_replicates <- length(unique(scseq$project)) > 2
  if (has_replicates) {
    pbulk <- diff_expr_pbulk(scseq = scseq,
                             data_dir = file.path(data_dir, anal_name),
                             clusters_name = clusters_name)
  }

  ebayes_sv <- fit_ebayes_scseq(scseq, ident.1 = 'test', ident.2 = 'ctrl')
  tt <- limma::topTable(ebayes_sv, coef = 1, number = Inf)

  # need ebayes_sv and pdata for add_es to get dprimes and vardprimes
  pdata <- scseq[['orig.ident']]
  pdata$group <- as.character(pdata$orig.ident)
  pdata$orig.ident <- NULL

  anal <- list(top_table = tt,
               ebayes_sv = ebayes_sv,
               pdata = pdata)

  # save to disk
  save_name <- paste("diff_expr_symbol_scseq", clusters_name, sep = "_")
  save_name <- paste0(save_name, ".rds")

  saveRDS(anal, file.path(data_dir, anal_name, save_name))
  return(anal)
}



#' Run pseudo bulk differential expression  analysis
#'
#' @param scseq \code{Seurat} object
#' @param data_dir Path to folder to save analysis in.
#' @param clusters_name String with comma seperated selected clusters.
diff_expr_pbulk <- function(scseq, data_dir, clusters_name) {

  summed <- scater::sumCountsAcrossCells(scseq[['RNA']]@counts, scseq$project)
  summed <- summed[edgeR::filterByExpr(summed), ]
  summed <- as.matrix(summed)

  # get pdata
  pdata <- scseq@meta.data %>%
    dplyr::select(project, orig.ident) %>%
    dplyr::distinct() %>%
    dplyr::rename(group = orig.ident)

  # get normalization factors
  row.names(pdata) <- pdata$project
  pdata <- pdata[colnames(summed), ]
  pdata$lib.size <- colSums(summed)
  pdata$norm.factors <- edgeR::calcNormFactors(summed)

  # construct eset
  fdata <- data.frame(SYMBOL = row.names(summed), row.names = row.names(summed))
  eset <- Biobase::ExpressionSet(summed,
                                 phenoData = Biobase::AnnotatedDataFrame(pdata),
                                 featureData = Biobase::AnnotatedDataFrame(fdata))

  anal <- diff_expr(eset, data_dir = data_dir, anal_name = paste0('pbulk_', clusters_name), prev_anal = list(pdata = pdata))
  return(anal)
}

#' Omit ambient genes from analysis or top table
#'
#' Must provide one of anal or top_table
#'
#' @param scseq \code{Seurat} object
#' @param anal Result of \code{diff_expr_scseq}
#' @param top_table Usually from \code{anal$top_table}
#'
#' @return One of \code{anal} or \code{top_table} with entries corresponding to ambient genes omited.
#' @export
#' @keywords internal
ambient.omit <- function(scseq, markers, cluster_markers) {
  ambient.genes <- get_ambient(scseq, markers = markers, cluster_markers = cluster_markers)

  is.ambient <- row.names(markers) %in% ambient.genes
  markers <- markers[!is.ambient, ]

  return(markers)
}


#' Fit limma eBayes on single cell RNA-seq dataset
#'
#' @param scseq \code{Seurat} object.
#' @param ident.1 String with control group name in \code{Idents(scseq)}.
#' @param ident.2 String with test group name in \code{Idents(scseq)}.
#'
#' @return result of call to \code{\link[limma]{eBayes}}
#' @export
#' @keywords internal
fit_ebayes_scseq <- function(scseq, ident.1, ident.2) {
  contrast = paste0(ident.1, '-', ident.2)

  # use SCT corrected logcounts
  dat <- scseq[['SCT']]@data
  group <- Seurat::Idents(scseq)
  design <- stats::model.matrix(~0 + group)
  colnames(design) <- levels(group)
  fit <- limma::lmFit(dat, design)
  cont.matrix <- limma::makeContrasts(contrasts = contrast, levels = design)
  fit <- limma::contrasts.fit(fit, cont.matrix)
  return(limma::eBayes(fit))
}

#' Get ambient genes to exclude for diff_expr_scseq
#'
#' Excludes upregulated genes if they are ambient in the test group.
#' Excludes downregulated genes if they are ambient in the control group.
#'
#' @param scseq \code{Suerat} object.
#' @param top_table \code{data.frame} with results from \code{\link[limma]{topTable}}.
#'
#' @return Character vector of HGNC symbols representing ambient genes to exclude
#' @export
#' @keywords internal
get_ambient <- function(scseq, markers, cluster_markers) {

  fts <- scseq[['SCT']]@meta.features
  rns <- row.names(markers)
  fts <- fts[rns, ]
  test.ambient <- rns[fts$test_ambient]
  ctrl.ambient <- rns[fts$ctrl_ambient]

  # exclude test/ctrl ambient if positive/negative effect size
  # opposite would decrease extent of gene expression difference but not direction
  pos.test <- markers[test.ambient, 't'] > 0
  neg.ctrl <- markers[ctrl.ambient, 't'] < 0

  test.exclude <- test.ambient[pos.test]
  ctrl.exclude <- ctrl.ambient[neg.ctrl]

  ambient <- c(test.exclude, ctrl.exclude)

  # dont include cluster markers in ambient
  ambient <- setdiff(ambient, row.names(cluster_markers))

  return(ambient)
}

#' Get direction of pathways
#'
#' @param top_table
#'
#' @return a \code{data.frame} with columns: \itemize{
#'  \item is.up Boolean indicating if direction is mostly up
#'  \item label Character vector with values \code{'Mostly Up'} or \code{'Mostly Down'}.
#'  These are used for selectize input option group field.
#' }
#' @export
get_path_directions <- function(top_table) {

  is.up <- sapply(gslist.kegg, function(gs) {
    in.gs <- row.names(top_table) %in% names(gs)
    mean(top_table[in.gs, 't']) > 0
  })

  is.up <- is.up[!is.na(is.up)]

  data.frame(label = ifelse(is.up, 'Mostly Up', 'Mostly Down'),
             is.up = is.up,
             row.names = names(is.up),
             stringsAsFactors = FALSE)
}

#' Toggle disabled state for multiple ids
#'
#' @param ids Character vector of ids to disable
#' @export
#' @keywords internal
toggleAll <- function(ids){
  for(id in ids) shinyjs::toggleState(id)
}
