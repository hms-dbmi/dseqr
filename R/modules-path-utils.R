

#' Used by get_all_df and get_path_df to construct the return result
#'
#' @param top_table Filtered result of \code{\link[limma]{topTable}}
#' @param nmax Maximum number of rows to keep. Default is all. Used be \code{construct_all_df}
#' to limit number of plotted genes.
#' @export
#' @keywords internal
construct_path_df <- function(top_table, nmax = nrow(top_table)) {

  # show up to nmax genes
  nkeep <- min(nmax, nrow(top_table))

  path_df <- data.frame(
    Gene = row.names(top_table),
    Dprime = top_table$dprime,
    sd = sqrt(top_table$vardprime),
    Link = paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", row.names(top_table), "'>", row.names(top_table), "</a>"), stringsAsFactors = FALSE
  )

  path_df <- path_df %>%
    arrange(desc(abs(Dprime))) %>%
    mutate(Gene = factor(Gene, levels = Gene)) %>%
    head(nkeep)


  return(path_df)
}


#' Get data.frame for plotting gene expression values of top genes
#'
#' @param anal Result of call to \code{\link{diff_expr_scseq}}
#' @param show_up Boolean. if \code{TRUE}, will return only upregulated genes.
#' If \code{FALSE} will return only downregulated genes.
#'
#' @return \code{data.frame} with columns: \itemize{
#'  \item Gene gene names.
#'  \item Dprime standardized unbiased effect size values.
#'  \item sd standard deviations of \code{Dprime}.
#'  \item Link url to GeneCards page for gene.
#' }
#' @export
get_all_df <- function(anal, show_up) {

  # add dprimes and vardprime values
  anal <- add_es(anal)
  top_table <- anal$top_table

  is.up <- top_table$t > 0
  filter <- if (show_up) is.up else !is.up

  top_table <- top_table[filter, ]

  construct_path_df(top_table, nmax = 200)
}


#' Get data.frame for plotting gene expression values of a pathway
#'
#' @param path_id String with KEGG pathway id.
#' @param anal Result of call to \code{\link{diff_expr_scseq}}
#'
#' @return \code{data.frame} with columns: \itemize{
#'  \item Gene gene names.
#'  \item Dprime standardized unbiased effect size values.
#'  \item sd standard deviations of \code{Dprime}.
#'  \item Link url to GeneCards page for gene.
#' }
#' @export
#' @keywords internal
get_path_df <- function(path_id, anal) {

  # add dprimes and vardprime values
  anal <- add_es(anal)
  top_table <- anal$top_table

  path_enids <- gslist[[path_id]]
  path_genes <- names(path_enids)

  # subset top table to genes in the pathway
  top_table <- top_table[row.names(top_table) %in% path_genes, ]

  construct_path_df(top_table)
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
#' Analysis is performed comparing test to control cells for all included clusters.
#'
#' @param scseq \code{Seurat} object
#' @param prev_anal Result of call to \code{diff_expr_scseq}
#' @param data_dir Directory to folders with single cell analyses
#' @param anal_name Name of folder in \code{data_dir} to save results to
#' @param clusters Character vector of clusters in \code{scseq$seurat_clusters} to include
#'
#' @return result of \code{\link[PADOG]{padog}}
#' @export
#' @keywords internal
diff_path_scseq <- function(scseq, prev_anal, data_dir, anal_name, clusters) {

  # load previous if exists
  clusters_name <- paste(sort(clusters), collapse = ',')
  fname <- paste0('diff_path_', clusters_name, '.rds')
  fpath <- file.path(data_dir, anal_name, fname)

  if(file.exists(fpath)) return(readRDS(fpath))

  Seurat::DefaultAssay(scseq) <- 'SCT'
  Seurat::Idents(scseq) <- scseq$orig.ident

  # subset to analysed gened (excludes ambient)
  genes <- row.names(prev_anal$top_table)
  in.clusters <- scseq$seurat_clusters %in% clusters

  scseq <-  scseq[genes, in.clusters]

  # get groups
  group <- as.character(scseq$orig.ident)
  group <- ifelse(group == 'ctrl', 'c', 'd')

  # expression matrix
  esetm  <- scseq[['SCT']]@data

  # already annotated with hgnc symbols
  gslist <- lapply(gslist, function(gs) {ret <- names(gs); names(ret) <- ret; return(ret)})

  # run padog
  padog_table <- PADOG::padog(esetm = esetm, group = group, parallel = TRUE, ncr = 4, gs.names = gs.names, gslist = gslist,
                              verbose = FALSE, rna_seq = FALSE)

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
diff_expr_scseq <- function(scseq, clusters, pseudo_bulk = FALSE) {
  Seurat::Idents(scseq) <- scseq$orig.ident

  # exclude non-selected clusters
  scseq <-  scseq[, scseq$seurat_clusters %in% clusters]

  if (pseudo_bulk) {
    # pseudo bulk counts
    summed <- scater::sumCountsAcrossCells(scseq[['SCT']]@counts, scseq$project)
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

    anal <- diff_expr(eset, getwd(), anal_name = 'blah', prev_anal = list(pdata = pdata))
    return(anal)

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
ambient.omit <- function(scseq, anal = NULL, top_table = anal$top_table) {
  ambient.genes <- get_ambient(scseq, top_table)
  is.ambient <- row.names(top_table) %in% ambient.genes
  top_table <- top_table[!is.ambient, ]

  if (!is.null(anal)) {
    anal$top_table <- top_table
    # need to make work with add_es
    anal$ebayes_sv$df.residual <- anal$ebayes_sv$df.residual[!is.ambient]
    return(anal)

  } else {
    return(top_table)
  }
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
get_ambient <- function(scseq, top_table) {

  fts <- scseq[['SCT']]@meta.features
  test.ambient <- row.names(fts)[fts$test_ambient]
  ctrl.ambient <- row.names(fts)[fts$ctrl_ambient]

  # exclude test/ctrl ambient if positive/negative effect size
  # opposite would decrease extent of gene expression difference but not direction
  pos.test <- top_table[test.ambient, 't'] > 0
  neg.ctrl <- top_table[ctrl.ambient, 't'] < 0

  test.exclude <- test.ambient[pos.test]
  ctrl.exclude <- ctrl.ambient[neg.ctrl]

  ambient <- c(test.exclude, ctrl.exclude)

  return(unique(ambient))
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

  is.up <- sapply(gslist, function(gs) {
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

