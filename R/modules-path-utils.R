#' Get perturbation signature
#'
#' @param pert Name of perturbation signature.
#' @param pert_type One of \code{'cmap'}, \code{'l1000_genes'}, or \code{'l1000_drugs'}.
#' @export
#' @keywords internal
load_pert_signature <- function(pert, pert_type, pert_signature_dir) {
  sig <- NULL
  type_dir <- file.path(pert_signature_dir, pert_type)
  if (!file.exists(type_dir)) dir.create(type_dir)


  fname <- paste0(pert, '.rds')
  sig_path <- file.path(type_dir, fs::path_sanitize(fname))
  if (!file.exists(sig_path)) dl_pert_signature(sig_path, pert_type)
  if (file.exists(sig_path)) sig <- readRDS(sig_path)
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

  # don't error if pert_type updates but sig_path corresponds to previous pert_type
  try(download.file(dl_url, sig_path))
}

#' Used by get_path_df to construct the return result
#'
#' @param top_table Filtered result of \code{\link[limma]{toptable}}
#' to limit number of plotted genes.
#' @export
#' @keywords internal
construct_path_df <- function(top_table) {

  data.frame(
    Gene = row.names(top_table),
    Dprime = top_table$dprime,
    sd = sqrt(top_table$vardprime),
    description = tx2gene$description[match(row.names(top_table), tx2gene$gene_name)],
    Link = paste0("<a class='xaxis-tooltip' href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", row.names(top_table), "'>", row.names(top_table), "</a>"),
    stringsAsFactors = FALSE
  ) %>%
    mutate(Gene = factor(Gene, levels = Gene))
}


#' Get data.frame for plotting gene expression values of a pathway
#'
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
get_path_df <- function(top_table, path_id = NULL, pert_signature = NULL, nmax = 200, ambient = NULL) {

  # single cell ambient genes excluded from query
  is.ambient <- row.names(top_table) %in% ambient
  top_table <- top_table[!is.ambient, ]

  # only show pathway if in GO
  if (path_id %in% names(gslist.go)) {
    path_enids <- gslist.go[[path_id]]
    path_genes <- names(path_enids)
    top_table <- top_table[row.names(top_table) %in% path_genes, ]
    top_table <- head(top_table, nmax)
  }

  path_df <- construct_path_df(top_table)
  path_df$point_color <- 'black'

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

    # show pert as arrow going to sum of disease + pert
    path_df$dprime_sum <- path_df$Dprime + pert_signature

    # get arrow colors based on if drug pushes gene in right direction
    path_signs <- sign(path_df$Dprime)
    pert_signs  <- sign(pert_signature)
    path_df$arrow_color <- ifelse(pert_signs == path_signs, '#f64f5a', '#1863e6')

  }

  return(path_df)
}

#' Generate data.frame of integrated single cell RNA-Seq datasets
#'
#' Used together with \code{\link{load_bulk_datasets}} to create data.frame
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
load_scseq_datasets <- function(data_dir) {
  int_path <- file.path(data_dir, 'single-cell', 'integrated.rds')

  datasets <- data.frame(matrix(ncol = 5, nrow = 0), stringsAsFactors = FALSE)
  colnames(datasets) <- c("dataset_name", "dataset_dir", "label", "value", "type")

  if (file.exists(int_path)) {
    integrated <- readRDS(int_path)
    for(ds in integrated) datasets[nrow(datasets) + 1, ] <- c(ds, file.path('single-cell', ds), ds, NA, 'Single Cell')
  }

  return(datasets)
}


#' Run limma fit for clusters in single cell RNA-Seq dataset
#'
#' @param obj \code{SingleCellExperiment} object
#' @param dataset_dir Path to folder to save results and pseudobulk eset to
#' @param is_summed Boolean indicating if \code{obj} is a pseudobulk \code{SingleCellExperiment} object.
#'
#' @return Named list with slots: \itemize{
#'  \item fit result of \link[limma]{lmFit}.
#'  \item mod \code{model.matrix} used for \code{fit}.
#' }
#' @export
#' @keywords internal
run_limma_scseq <- function(obj) {

  if ('ExpressionSet' %in% class(obj)) {
    # run as pseudobulk
    lm_fit <- run_limma(obj, prev_anal = list(pdata = Biobase::pData(obj)))
    lm_fit$has_replicates <- TRUE

  } else {
    lm_fit <- fit_lm_scseq(obj)
    lm_fit$has_replicates <- FALSE
  }

  return(lm_fit)
}

#' Loads markers for single cell cluster
#'
#' @param selected_clusters Character vector of integers
#' @param dataset_dir Directory to folder with single-cell dataset

#' @return data.frame with markers
#' @export
#' @keywords internal
get_cluster_markers <- function(selected_clusters, dataset_dir) {

  # get cluster markers
  clusters_name <- collapse_sorted(selected_clusters)
  fname <- paste0("markers_", clusters_name, '.rds')
  fpath <- file.path(dataset_dir, fname)

  if (file.exists(fpath)) {
    cluster_markers <- readRDS(fpath)

  } else {
    # get markers for multi-cluster selections
    browser()

  }

  return(cluster_markers)
}



#' Run pseudo bulk limma fit
#'
#' @param summed Pseudobulk \code{SingleCellExperiment}
#' @return Normalized \code{ExpressionSet}
#'
#'@export
#'@keywords internal
construct_pbulk_eset <- function(summed, species = 'Homo sapiens', release = '94') {

  # discard outlier samples
  y <- edgeR::DGEList(counts(summed), samples = summed@colData)
  outl <- scater::isOutlier(y$samples$lib.size, log=TRUE, type="lower")
  y <- y[, !outl]

  # filter genes
  keep <- edgeR::filterByExpr(y, group = summed$cluster[!outl])
  y <- y[keep, ]

  # normalize for composition
  y <- edgeR::calcNormFactors(y)

  # construct eset
  annot <- get_ensdb_package(species, release)
  fdata <- setup_fdata(species, release)
  eset <- construct_eset(y, fdata, annot)

  # use e.g. test_1 as sample label (test sample cluster 1)
  pdata <- Biobase::pData(eset)
  pdata$group <- paste(pdata$orig.ident, pdata$cluster, sep = '_')
  Biobase::pData(eset) <- pdata

  # add vst transformed values
  eset <- add_vsd(eset, pbulk = TRUE)
  return(eset)
}

#' Move genes to bottom of markers data.frame
#'
#' Used to supress ambient outliers and positive genes for single-cell sample comparisons
#'
#' @param markers data.frame with genes as row names.
#' @param supress Genes to move to bottom of \code{markers}.
#'
#' @return \code{markers} with genes in \code{supress} at bottom.
#' @export
#' @keywords internal
supress.genes <- function(markers, supress) {
  other <- setdiff(row.names(markers), supress)
  markers[c(other, supress), ]
}


#' Fit limma eBayes on single cell RNA-seq dataset
#'
#' @param scseq \code{SingleCellExperiment} object.
#'
#' @return List with model matrix and result of call to \code{\link[limma]{lmFit}}
#' @export
#' @keywords internal
fit_lm_scseq <- function(scseq) {

  # use multiBatchNorm logcounts
  dat <- SingleCellExperiment::logcounts(scseq)
  group <- paste(scseq$orig.ident, as.numeric(scseq$cluster), sep = '_')
  mod <- stats::model.matrix(~0 + group)
  colnames(mod) <- gsub('^group', '', colnames(mod))
  fit <- limma::lmFit(dat, mod)

  # add enids for goana/kegga pathway analyses
  data.table::setkey(hs, SYMBOL_9606)
  fit$genes <- hs[row.names(scseq), list(ENTREZID)]

  return(list(fit = fit, mod = mod))
}

#' Get ambient genes to exclude
#'
#' Excludes upregulated genes if they are ambient in the test group.
#' Excludes downregulated genes if they are ambient in the control group.
#'
#' @param scseq \code{SingleCellExperiment} object.
#'
#' @return Character vector of HGNC symbols representing ambient genes to exclude
#' @export
#' @keywords internal
decide_ambient <- function(ambient, top_table, cluster_markers) {

  # exclude test ambient if up
  # exclude ctrl ambient if down
  # opposite would decrease extent of gene expression difference but not direction
  pos.test <- top_table[ambient$test, 'logFC'] > 0
  neg.ctrl <- top_table[ambient$ctrl, 'logFC'] < 0

  test.exclude <- ambient$test[pos.test]
  ctrl.exclude <- ambient$ctrl[neg.ctrl]

  ambient <- c(test.exclude, ctrl.exclude)

  # exclude top cluster markers from ambient
  # asumes in-cell >> ambient for markers
  keep <- cluster_markers$FDR < 0.05
  ambient <- setdiff(ambient, row.names(cluster_markers)[keep])

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

#' Get heatmap of genes in a GO pathways
#'
#' @param eset ExpressionSet
#' @param top_table limma topTable
#' @param contrast_groups groups to include in heatmap
#' @param path_id GO pathway id string
#'
#' @return \code{pheatmap} object.
#' @export
get_pathway_heatmap <- function(eset, top_table, contrast_groups, path_id) {

  # subset to pathway genes
  adj <- Biobase::assayDataElement(eset, 'adjusted')
  path_genes <- names(gslist.go[[path_id]])

  # subset to max 1000 top genes
  sig.genes <- head(row.names(top_table[path_genes, ]), 1000)
  adj <- adj[row.names(adj) %in% sig.genes, ]

  # get group colors and levels
  pdata <- Biobase::pData(eset)
  group_levels <- get_group_levels(pdata)
  group_colors <- get_group_colors(group_levels)

  # subset eset to selected groups
  in.sel <- pdata$group %in% contrast_groups
  adj <- adj[, in.sel]
  if ((nrow(adj) %% 2) == 1) adj <- rbind(adj, rnorm(ncol(adj)))

  # cluster rows on correlation
  dist.rows <- as.dist(1 - cor(t(adj)))
  annot <- data.frame(Group = pdata$group[in.sel], row.names = colnames(adj))
  annotation_colors <- list(Group = group_colors)
  names(annotation_colors$Group) <- group_levels
  annotation_colors$Group <- annotation_colors$Group[contrast_groups]

  pheatmap::pheatmap(adj,
                     color = gplots::greenred(10),
                     clustering_distance_rows = dist.rows,
                     show_colnames = TRUE,
                     show_rownames = FALSE,
                     border_color = NA,
                     scale = 'row',
                     treeheight_row = 0,
                     annotation_colors = annotation_colors,
                     annotation_names_col = FALSE,
                     annotation_col = annot,na_col = 'black',
                     cellwidth=18, fontsize = 12)


}
