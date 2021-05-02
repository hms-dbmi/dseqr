.datatable.aware = TRUE

#' Get perturbation signature
#'
#' @param pert Name of perturbation signature.
#' @param pert_type One of \code{'cmap'}, \code{'l1000_genes'}, or \code{'l1000_drugs'}.
#' @param pvals If \code{TRUE} returns adjusted pvalues for signature. If \code{FALSE} (default) returns signature.
#'
#' @keywords internal
load_pert_signature <- function(pert, pert_type, pert_signature_dir, pvals = FALSE) {
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
#' @keywords internal
#'
#' @return NULL
#'
#' @examples
#' data_dir <- tempdir()
#' sig_path <- file.path(data_dir, '0317956-0000_PC3_1e-06M_6h.rds')
#' dseqr:::dl_pert_signature(sig_path, pert_type = 'cmap')
#'
dl_pert_signature <- function(sig_path, pert_type) {
  # name of the file being requested
  dl_url <- paste0('https://s3.us-east-2.amazonaws.com/dseqr/drug_es_dir/', pert_type, '/', basename(sig_path))
  dl_url <- utils::URLencode(dl_url)
  dl_url <- gsub('+', '%2B', dl_url, fixed = TRUE)

  # don't error if pert_type updates but sig_path corresponds to previous pert_type
  try(utils::download.file(dl_url, sig_path))
}


#' Used by get_path_df to construct the return result
#'
#' @param top_table Filtered result of \code{\link[limma]{toptable}}
#' to limit number of plotted genes.
#'
#' @keywords internal
construct_path_df <- function(top_table) {
  Gene <- NULL

  tx2gene <- dseqr.data::load_tx2gene()

  df <- data.frame(
    Gene = row.names(top_table),
    Dprime = signif(top_table$dprime, digits = 3),
    stringsAsFactors = FALSE)

  # split up assignment because custom query signatures may not have all
  df$logfc <- safe.fun(signif, top_table$logFC, digits = 3)
  df$sd <- safe.fun(signif, safe.fun(sqrt, top_table$vardprime), digits = 3)
  df$pval <- safe.fun(format.pval, top_table$P.Value, eps = 0.005, digits = 2)
  df$fdr <- safe.fun(format.pval, top_table$adj.P.Val, eps = 0.005, digits = 2)
  df$description <- tx2gene$description[match(row.names(top_table), tx2gene$gene_name)]
  df$Link <- paste0("<a class='xaxis-tooltip' ",
                   "href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                   row.names(top_table), "'>", row.names(top_table), "</a>")

  dplyr::mutate(df, Gene = factor(Gene, levels = Gene))
}

safe.fun <- function(fun, x, ...) {
  if (!length(x)) return(NULL)
  fun(x, ...)

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
#'
#' @keywords internal
get_path_df <- function(top_table, path_id = NULL, pert_signature = NULL, nmax = 200, ambient = NULL) {

  # single cell ambient genes excluded from query
  is.ambient <- row.names(top_table) %in% ambient
  top_table <- top_table[!is.ambient,, drop=FALSE]

  path_df <- construct_path_df(top_table)
  path_df$color <- 'black'

  # for drug/genetic queries: keep up to nmax in cmap/l1000
  is.cmap2 <- path_id == 'CMAP02'
  is.l1000 <- path_id == 'L1000'

  if (is.l1000 | is.cmap2) {
    if (is.cmap2) keep <- utils::head(which(path_df$Gene %in% unlist(genes)), nmax)
    if (is.l1000) keep <- utils::head(which(path_df$Gene %in% genes$common), nmax)

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
#'
#' @keywords internal
load_scseq_datasets <- function(data_dir) {
  sc_dir <- file.path(data_dir, 'single-cell')
  int_path <- file.path(sc_dir, 'integrated.qs')

  datasets <- data.frame(matrix(ncol = 5, nrow = 0), stringsAsFactors = FALSE)
  colnames(datasets) <- c("dataset_name", "dataset_dir", "label", "value", "type")

  if (file.exists(int_path)) {
    integrated <- qread.safe(int_path)
    has.scseq <- check_has_scseq(integrated, sc_dir)
    integrated <- integrated[has.scseq]

    for(ds in integrated) datasets[nrow(datasets) + 1, ] <- c(ds, file.path('single-cell', ds), ds, NA, 'Single Cell')
  }

  return(datasets)
}


#' Run limma fit for clusters in single cell RNA-Seq dataset
#'
#' @param obj \code{SingleCellExperiment} or list of \code{ExpressionSet} objects for pseudobulk.
#' @param annot Annotation to use. Default \code{'gene_name'} works for both mouse and human.
#'
#' @return Named list with slots: \itemize{
#'  \item fit result of \link[limma]{lmFit}.
#'  \item mod \code{model.matrix} used for \code{fit}.
#' }
#'
#' @keywords internal
run_limma_scseq <- function(obj, annot = 'gene_name') {

  if (is.list(obj)) {
    # run as pseudobulk per eset
    lm_fit <- lapply(obj, function(x) {
      x <- crossmeta::run_limma_setup(x, list(pdata = Biobase::pData(x)))
      crossmeta::run_limma(x,  annot = annot, filter = FALSE, quick = TRUE)
    })

  } else {
    lm_fit <- fit_lm_scseq(obj)
  }

  return(lm_fit)
}

#' Loads markers for single cell cluster
#'
#' @param selected_clusters Character vector of integers
#' @param dataset_dir Directory to folder with single-cell dataset

#' @return data.frame with markers
#'
#' @keywords internal
get_cluster_markers <- function(selected_clusters, dataset_dir) {

  # get cluster markers
  clusters_name <- collapse_sorted(selected_clusters)
  fname <- paste0("markers_", clusters_name, '.qs')
  fpath <- file.path(dataset_dir, fname)

  if (file.exists(fpath)) {
    cluster_markers <- qs::qread(fpath)

  }

  return(cluster_markers)
}



#' Run pseudo bulk limma fit
#'
#' @param summed Pseudobulk \code{SingleCellExperiment}
#' @return Normalized \code{ExpressionSet}
#'
#'
#'@keywords internal
construct_pbulk_esets <- function(summed, species = 'Homo sapiens') {

  release <- switch(species,
                    'Homo sapiens' = '94',
                    'Mus musculus' = '98')

  y <- edgeR::DGEList(SingleCellExperiment::counts(summed), samples = summed@colData)
  clusters <- y$samples$cluster

  annot <- dseqr.data::get_ensdb_package(species, release)
  fdata <- rkal::setup_fdata(species, release)

  esets <- list()
  for (clust in levels(clusters)) {
    yi <- y[, clusters == clust]

    # skip if no replicates
    group <- yi$samples$orig.ident
    neach <- table(group)
    if (any(neach == 0) || sum(neach) < 3) next

    # filter genes
    keep <- edgeR::filterByExpr(yi, group = group)
    yi <- yi[keep, ]

    # skip if less than 2 genes (required by voom)
    if (nrow(yi) < 2) next

    # normalize for composition
    yi <- edgeR::calcNormFactors(yi)

    # construct eset
    eset <- rkal::construct_eset(yi, fdata, annot)

    # use test and ctrl as group
    eset$group <- eset$orig.ident

    # add vst transformed values
    eset <- add_vsd(eset, pbulk = TRUE)
    esets[[clust]] <- eset
  }

  return(esets)
}


#' Fit limma eBayes on single cell RNA-seq dataset
#'
#' @param scseq \code{SingleCellExperiment} object.
#'
#' @return List with model matrix and result of call to \code{\link[limma]{lmFit}}
#'
#' @keywords internal
fit_lm_scseq <- function(scseq) {
  SYMBOL_9606 <- ENTREZID <- NULL

  hs <- readRDS(system.file('extdata', 'hs.rds', package = 'dseqr.data'))
  data.table::setkey(hs, SYMBOL_9606)

  # one fit per cluster
  fits <- list()
  clusters <- scseq$cluster
  species <- scseq@metadata$species

  for (clust in levels(clusters)) {

    # use multiBatchNorm logcounts
    yi <- scseq[, clusters == clust]
    dat <- SingleCellExperiment::logcounts(yi)
    dat <- dat[Matrix::rowSums(dat) > 0, ]

    # require at least 5 cells per group
    group <- yi$orig.ident
    neach <- table(group)
    if (any(neach < 5)) next

    mod <- stats::model.matrix(~0 + group)
    colnames(mod) <- gsub('^group', '', colnames(mod))
    fit <- limma::lmFit(dat, mod)

    # add enids for go/kegg pathway analyses
    rn <- row.names(dat)
    if (species == 'Homo sapiens') {
      fit$genes <- hs[rn, list(ENTREZID)]

    } else if (species == 'Mus musculus') {
      tx2gene_mouse <- dseqr.data::load_tx2gene(species)
      tx <- tx2gene_mouse
      tx$ENTREZID <- as.character(tx$entrezid)
      ind <- match(rn, tx2gene_mouse$gene_name)
      fit$genes <- tx[ind, 'ENTREZID']
    }

    fits[[clust]] <- list(fit = fit, mod = mod)
  }
  return(fits)
}



#' Get gslist for pathway analysis
#'
#' @param species Species identifier
#' @param type either 'go' or 'kegg'
#' @param gs_dir Directory to save gslist to
#'
#' @return gslist
#' @keywords internal
get_gslist <- function(species = 'Hs', type = 'go', gs_dir = '/srv/dseqr/gs_dir') {
  PathwayID <- GeneID <- SYMBOL <- NULL

  if (!dir.exists(gs_dir)) dir.create(gs_dir)

  fname <- paste('gslist', type, species, 'qs', sep = '.')
  gslist_path <- file.path(gs_dir, fname)

  if (file.exists(gslist_path)) {
    gslist <- qs::qread(gslist_path)

  } else if (type == 'go') {
    orgPkg <- paste0("org.",species,".eg.db")
    require(orgPkg, character.only = TRUE, quietly = TRUE)

    #	Get access to package of GO terms
    suppressPackageStartupMessages(OK <- requireNamespace("GO.db",quietly=TRUE))
    if(!OK) stop("GO.db package required but is not installed (or can't be loaded)")

    #	Get access to required annotation functions
    suppressPackageStartupMessages(OK <- requireNamespace("AnnotationDbi",quietly=TRUE))
    if(!OK) stop("AnnotationDbi package required but is not installed (or can't be loaded)")

    #	Load appropriate organism package
    suppressPackageStartupMessages(OK <- requireNamespace(orgPkg,quietly=TRUE))
    if(!OK) stop(orgPkg," package required but is not installed (or can't be loaded)")

    #	Get GO to Entrez Gene mappings
    obj <- paste0("org.",species,".egGO2ALLEGS")
    egGO2ALLEGS <- tryCatch(utils::getFromNamespace(obj,orgPkg), error=function(e) FALSE)
    if(is.logical(egGO2ALLEGS)) stop("Can't find gene ontology mappings in package ",orgPkg)

    # Entrez gene to symbol
    # TODO for Hs: get from toupper(hs[EG.GO$gene_id, SYMBOL_9606]) so that consistent with original annotation
    EG.GO <- AnnotationDbi::toTable(egGO2ALLEGS)
    EG.GO$SYMBOL <- AnnotationDbi::mapIds(get(orgPkg), EG.GO$gene_id, column = 'SYMBOL', keytype = 'ENTREZID')
    gslist <- split(EG.GO, EG.GO$go_id)
    gslist <- lapply(gslist, function(df) {tmp <- df$gene_id; names(tmp) <- df$SYMBOL; tmp})
    qs::qsave(gslist, gslist_path)

  } else if (type == 'kegg') {

    orgPkg <- paste0("org.",species,".eg.db")
    require(orgPkg, character.only = TRUE, quietly = TRUE)

    kegg_species <- get_kegg_species(species)
    gkl <- limma::getGeneKEGGLinks(kegg_species, convert = TRUE)
    # TODO for Hs: get from toupper(hs[EG.GO$gene_id, SYMBOL_9606]) so that consistent with original annotation
    gkl$SYMBOL <- AnnotationDbi::mapIds(get(orgPkg), gkl$GeneID, column = 'SYMBOL', keytype = 'ENTREZID')
    gkl <- gkl %>%
      dplyr::group_by(PathwayID) %>%
      dplyr::summarise(gslist = list(GeneID),
                       symbols = list(SYMBOL))

    gslist <- gkl$gslist
    symbols <- gkl$symbols
    gslist <- lapply(seq_along(gslist), function(i) {tmp <- gslist[[i]]; names(tmp) <- symbols[[i]]; tmp})

    names(gslist) <- gkl$PathwayID
    qs::qsave(gslist, gslist_path)

  }

  return(gslist)
}

#' Convert from species identifier to KEGG species identifier
#'
#' @param species species identifier
#'
#' @return KEGG species identifier
#' @keywords internal
get_kegg_species <- function(species) {
  species <- match.arg(species, c("Ag", "At", "Bt", "Ce", "Dm", "Dr", "EcK12", "EcSakai", "Gg", "Hs", "Mm", "Mmu", "Pf", "Pt", "Rn", "Ss", "Xl"))
  #	Convert from Bioconductor to KEGG species codes
  species.KEGG <- switch(species, "Ag"="aga", "At"="ath", "Bt"="bta", "Ce"="cel", "Cf"="cfa", "Dm"="dme", "Dr"="dre", "EcK12"="eco", "EcSakai"="ecs", "Gg"="gga", "Hs"="hsa", "Mm"="mmu", "Mmu"="mcc", "Pf"="pfa", "Pt"="ptr", "Rn"="rno", "Ss"="ssc", "Xl"="xla")

  return(species.KEGG)
}


#' Get names of gene set
#'
#' @param gslist result of \code{\link{get_gslist}}
#' @param type either 'go' or 'kegg'
#' @param species species identifier
#' @param gs_dir Directory to save results to
#'
#' @return Description of \code{gslist} gene sets
#' @keywords internal
#'
get_gs.names <- function(gslist, type = 'go', species = 'Hs', gs_dir = '/srv/dseqr/gs_dir') {
  if (!dir.exists(gs_dir)) dir.create(gs_dir)

  fname <- paste('gs.names', type, species, 'qs', sep = '.')
  gs.names_path <- file.path(gs_dir, fname)

  if (file.exists(gs.names_path)) {
    gs.names <- qs::qread(gs.names_path)

  } else if (type == 'go') {
    GOID <- names(gslist)
    TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db,keys=GOID,columns="TERM"))
    gs.names <- TERM$TERM
    names(gs.names) <- TERM$GOID
    qs::qsave(gs.names, gs.names_path)

  } else if (type == 'kegg') {
    kegg_species <- get_kegg_species(species)
    gs.names <- limma::getKEGGPathwayNames(kegg_species, remove=TRUE)
    row.names(gs.names) <- gs.names$PathwayID
    gs.names <- gs.names[names(gslist), 'Description']
    names(gs.names) <- names(gslist)
    qs::qsave(gs.names, gs.names_path)

  }

  return(gs.names)
}

