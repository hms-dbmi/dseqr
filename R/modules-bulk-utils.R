#' Validate pdata before differential expression analysis
#'
#' @param pdata data.frame with column \code{'Group'}
#'
#' @return NULL if valid, otherwise a character vector indicating what's wrong
#'
#' @keywords internal
validate_pdata <- function(pdata) {
  group <- pdata$Group
  group <- group[!is.na(group)]

  if (length(unique(group)) != 2) {
    msg <- 'Analysis requires test and control groups'

  } else if (length(group) < 3) {
    msg <- 'At least three samples are required for analysis'

  } else {
    msg <- NULL
  }
  return(msg)
}

#' Generate boxplotly for vsd normalized gene plots.
#'
#' @param df \code{data.frame} with columns: \itemize{
#'  \item x Factor for x-labels.
#'  \item y Numeric column used for y-values.
#'  \item text Character column used for hoverinfo.
#'  \item name Character column used for legend names of \code{x} values.
#'  \item color Factor column used to generate ordered colors for boxplots for each \code{'x'} value.
#' }
#' @param boxgap Used for plotly layout.
#' @param boxgroupgap Used for plotly layout.
#' @param plot_fname Name to save plot as.
#' @param ytitle Y axis title.
#' @param xtitle X axis title.
#'
#' @return plotly
#'
#' @keywords internal
boxPlotly <- function(df, boxgap, boxgroupgap, plot_fname, ytitle, xtitle) {

  # legend styling
  l <- list(
    font = list(
      family = "sans-serif",
      size = 12,
      color = "#000"),
    bgcolor = "#f8f8f8",
    bordercolor = "#e7e7e7",
    borderwidth = 1)


  df %>%
    plotly::plot_ly() %>%
    plotly::add_trace(x = ~x,
                      y = ~y,
                      color = ~color,
                      colors = get_palette(levels(df$color)),
                      text = ~text,
                      name = ~name,
                      type = 'violin',
                      points = 'all',
                      jitter = 1,
                      pointpos = 0,
                      fillcolor = 'transparent',
                      hoverinfo = 'text',
                      hoveron = 'points',
                      marker = list(color = "rgba(0, 0, 0, 0.6)")) %>%
    plotly::layout(violinmode = 'group', violingroupgap = boxgroupgap, violingap = boxgap,
                   xaxis = list(fixedrange=TRUE, title = xtitle),
                   yaxis = list(fixedrange=TRUE, title = ytitle),
                   legend = l) %>%
    plotly::config(displaylogo = FALSE,
                   displayModeBar = 'hover',
                   modeBarButtonsToRemove = c('lasso2d',
                                              'select2d',
                                              'toggleSpikelines',
                                              'hoverClosestCartesian',
                                              'hoverCompareCartesian'),
                   toImageButtonOptions = list(format = "png", filename = plot_fname))

}

#' Generate boxplotly for cell type deconvolution plots.
#'
#' @param df \code{data.frame} with columns: \itemize{
#'  \item x Factor for x-labels.
#'  \item y Numeric column used for y-values.
#'  \item text Character column used for hoverinfo.
#'  \item name Character column used for legend names of \code{x} values.
#'  \item color Factor column used to generate ordered colors for boxplots for each \code{'x'} value.
#' }
#' @param boxgap Used for plotly layout.
#' @param boxgroupgap Used for plotly layout.
#' @param plot_fname Name to save plot as.
#' @param ytitle Y axis title.
#' @param xtitle X axis title.
#'
#' @return plotly
#'
#' @keywords internal
boxPlotlyCells <- function(df, boxgap, boxgroupgap, plot_fname, ytitle, xtitle) {

  # legend styling
  l <- list(
    font = list(
      family = "sans-serif",
      size = 12,
      color = "#000"),
    bgcolor = "#f8f8f8",
    bordercolor = "#e7e7e7",
    borderwidth = 1)

  # total width of plot
  nx <- length(unique(df$x))
  ng <- length(unique(df$name))
  plot_width <- max(923, (ng+2)*25*nx)


  df %>%
    plotly::plot_ly() %>%
    plotly::add_trace(x = ~x,
                      y = ~y,
                      color = ~color,
                      colors = get_palette(levels(df$color)),
                      text = ~text,
                      name = ~name,
                      type = 'box',
                      boxpoints = 'all',
                      jitter = 1,
                      pointpos = 0,
                      fillcolor = 'transparent',
                      hoverinfo = 'text',
                      hoveron = 'points',
                      marker = list(color = "rgba(0, 0, 0, 0.6)", size = 3)) %>%
    plotly::layout(boxmode = 'group', boxgroupgap = boxgroupgap, boxgap = boxgap,
                   width = plot_width,
                   height = 620,
                   xaxis = list(fixedrange=TRUE, title = xtitle),
                   yaxis = list(fixedrange=TRUE, title = ytitle),
                   legend = l) %>%
    plotly::config(displaylogo = FALSE,
                   displayModeBar = 'hover',
                   modeBarButtonsToRemove = c('lasso2d',
                                              'select2d',
                                              'toggleSpikelines',
                                              'hoverClosestCartesian',
                                              'hoverCompareCartesian'),
                   toImageButtonOptions = list(format = "png", filename = plot_fname))

}

#' Get arguments for gene boxplotly
#'
#' @param eset ExpressionSet object with \code{'adjusted'} assayDataElement
#' @param explore_genes Character vector of genes to plot
#' @param dataset_name Name of bulk dataset.
#'
#' @return List with items \code{'df'}, \code{'boxgap'}, \code{'boxgroupgap'}, and \code{'plot_fname'}.
#' @keywords internal
#'
get_boxplotly_gene_args <- function(eset, explore_genes, dataset_name) {
  dat <- Biobase::assayDataElement(eset, 'vsd')
  pdata <- Biobase::pData(eset)

  dfs <- list()
  for (gene in explore_genes) {
    dfs[[gene]] <- data.frame(text = row.names(pdata),
                              x = gene,
                              y = dat[gene, row.names(pdata)],
                              name = as.character(pdata$`Group name`),
                              color = pdata$Group,
                              stringsAsFactors = FALSE)
  }

  df <- do.call(rbind, dfs)

  group_levels <- as.character(sort(unique(df$color)))

  # adjust gaps within/between group based on number of boxs (chosen by trial and error)
  nbox <- length(group_levels) * length(explore_genes)
  boxgap <- ifelse(nbox > 5, 0.4, 0.6)
  boxgroupgap <- ifelse(nbox > 6, 0.3, 0.6)

  df$color <- factor(df$color, levels = group_levels)
  df$x <- factor(df$x, levels = explore_genes)

  # name for saving plot
  fname <- paste(explore_genes, collapse = '_')
  fname <- paste('bulk', dataset_name, fname, Sys.Date(), sep='_')

  return(list(df = df,
              boxgap = boxgap,
              boxgroupgap = boxgroupgap,
              plot_fname = fname))
}


#' Get arguments for cells boxplotly
#'
#' @param pdata data.frame of phenotype data.
#' @param dtangle_est data.frame of proportion estimates from dtangle.
#' @param dataset_name Name of bulk dataset.
#'
#' @return List with items \code{'df'}, \code{'boxgap'}, \code{'boxgroupgap'}, and \code{'plot_fname'}.
#' @keywords internal
#'
get_boxplotly_cell_args <- function(pdata, dtangle_est, dataset_name) {

  common <- intersect(row.names(dtangle_est), row.names(pdata))
  dtangle_est <- as.data.frame(dtangle_est)[common, ]
  pdata <- pdata[common, c('Title', 'Group name', 'Group')]
  colnames(pdata) <- c('text', 'name', 'color')

  df <- utils::stack(dtangle_est)
  colnames(df) <- c('y', 'x')

  df <- cbind(pdata, df)
  group_levels <- as.character(sort(unique(df$color)))
  clus_levels <- colnames(dtangle_est)

  # adjust gaps within/between group based on number of boxs (chosen by trial and error)
  nbox <- length(group_levels) * length(clus_levels)
  boxgap <- ifelse(nbox > 5, 0.4, 0.6)
  boxgroupgap <- ifelse(nbox > 6, 0.3, 0.6)

  df$color <- factor(df$color, levels = group_levels)
  df$x <- factor(df$x, levels = clus_levels)

  # name for saving plot
  fname <- paste(clus_levels, collapse = '_')
  fname <- gsub(' ', '', fname)
  fname <- paste(dataset_name, fname, Sys.Date(), sep='_')

  return(list(df = df,
              boxgap = boxgap,
              boxgroupgap = boxgroupgap,
              plot_fname = fname))
}



#' Load previous bulk datasets dataframe
#'
#' @param data_dir Path to folder container \code{'bulk'} and \code{'single-cell'} directories
#' @return data.frame with columns "dataset_name" and "dataset_dir"
#'
#' @keywords internal
load_bulk_datasets <-function(data_dir) {

  datasets <- data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = FALSE)
  colnames(datasets) <- c("dataset_name", "dataset_dir")

  dataset_names <- list.dirs(file.path(data_dir, 'bulk'), full.names = FALSE, recursive = FALSE)
  has.eset <- file.exists(file.path(data_dir, 'bulk', dataset_names, 'eset.qs'))
  dataset_names <- dataset_names[has.eset]

  datasets <- data.frame(dataset_name = dataset_names,
                         dataset_dir = file.path('bulk', dataset_names), stringsAsFactors = FALSE)

  datasets$value <-  datasets$label <- datasets$dataset_name
  if (nrow(datasets)) datasets$type <- datasets$group <- 'Bulk Data'

  return(datasets)
}


#' Delete stale dataset files after changing sample groups or running sva
#'
#' @param data_dir Path to folder with dataset files
#' @param patterns patterns to remove. Default is all.
#'
#'
#' @keywords internal
remove_dataset_files <- function(data_dir, patterns = c('^adjusted_\\d+svs.qs$',
                                                        '^iqr_keep_\\d+svs.qs$',
                                                        '^vsd.qs$',
                                                        '^svobj.qs$',
                                                        '^numsv.qs$',
                                                        '^lm_fit_\\d+svs.qs$',
                                                        '^mds_\\d+svs.qs$',
                                                        'cmap_res_.+_\\d+svs.qs$',
                                                        'go_.+_\\d+svs.qs$',
                                                        'goana_.+_\\d+svs.qs$',
                                                        'kegg_.+_\\d+svs.qs$',
                                                        'kegga_.+_\\d+svs.qs$',
                                                        'l1000_drugs_res_.+_\\d+svs.qs$',
                                                        'l1000_genes_res_.+_\\d+svs.qs$',
                                                        'diff_expr_symbol_.+_\\d+svs.qs$'), exclude = NULL) {

  #TODO: preface everything with hash based on group names so that don't have to delete

  for (pattern in patterns) {
    fpaths <- list.files(data_dir, pattern)
    if (!is.null(exclude)) {
      keep <- grepl(exclude, fpaths)
      fpaths <- fpaths[!keep]
    }
    unlink(file.path(data_dir, fpaths))
  }
}


#' Format uploaded annotation
#'
#' @keywords internal
format_up_annot <- function(up, ref) {
  row.names(up) <- up$Title
  up[up == ''] <- NA

  # Group in order of Group name
  # allows changing color of groups by changing order or samples
  group <- up$`Group name`
  levels <- unique(group[!is.na(group)])
  group <- as.numeric(factor(group, levels =levels))
  up <- tibble::add_column(up, Group = group, .before = 1)

  up$Pair <- factor(up$Pair)

  # in case order of sample was changed
  up <- up[row.names(ref), ]

  # restore rna seq specific things
  up$lib.size <- ref$lib.size
  up$norm.factors <- ref$norm.factors

  return(up)

}

#' Format downloaded annotation
#'
#' @keywords internal
format_dl_annot <- function(annot) {

  add_pair <- function(df) {
    pair <- df$Pair
    if (is.null(pair)) pair <- df$pair
    if (is.null(pair)) pair <- NA

    df$pair <- df$Pair <- NULL
    tibble::add_column(df, Pair = pair, .after = 'Title')
  }

  annot <- add_pair(annot)
  annot <- annot[, !colnames(annot) %in% c('Group', 'lib.size', 'norm.factors')]
  return(annot)

}

#' Validate uploaded bulk annotation
#'
#' @keywords internal
validate_up_annot <- function(up, ref) {
  msg <- NULL


  req_cols <- c('Title', 'Group name', 'Pair')
  miss_cols <- req_cols[!req_cols %in% colnames(ref)]

  group <- up$`Group name`
  group <- group[!is.na(group)]
  ngroups <- length(unique(group))

  if (length(miss_cols)) {
    msg <- paste('Missing columns:', paste(miss_cols, collapse = ', '))

  } else if (!all(up$Title %in% ref$Title)) {
    msg <- 'Do not change Title column'

  } else if (ngroups < 2) {
    msg <- 'Need at least 2 groups'

  } else if (length(group) < 3) {
    msg <- 'Need at least 3 grouped samples'

  } else if (!is_invertible(up)) {
    msg <- 'Group name and Pair combination not solvable'
  }

  return(msg)
}

#' Check uploaded bulk pdata to make sure the study design is invertible
#'
#' @keywords internal
is_invertible <- function(pdata) {
  pdata <- pdata[!is.na(pdata$`Group name`), ]

  pair <- pdata$Pair
  if (length(unique(pair)) > 1) pdata$pair <- pair

  pdata$group <- pdata$`Group name`

  mod <- crossmeta::get_sva_mods(pdata)$mod

  methods::is(try(solve.default(t(mod) %*% mod),silent=T), 'matrix')
}


#' Get and save pathway results for ebfit object
#'
#' Used to avoid code reuse for single-cell and bulk
#'
#' @param de Result of \code{\link[crossmeta]{fit_ebayes}} or
#'   \code{link[crossmeta]{get_top_table}}
#' @param goana_path Path to save goana Gene Ontology result
#'
#' @return List with GO and KEGG results
#'
#' @keywords internal
get_path_res <- function(de, goana_path, gs_dir, species = 'Hs', genego = NULL, gonames = NULL, coef = ncol(de), nmin = 50, cutoff = 0.05) {

  is.marray <- is(de, 'MArrayLM')
  if (!is.marray) {
    universe <- de$ENTREZID
    is.up <- de$logFC > 0
    is.sig <- de$adj.P.Val < cutoff

  } else {
    # order by pvals
    de <- de[order(de$p.value[, coef]), ]

    universe <- de$genes$ENTREZID
    is.up <- de$coefficients[, coef] > 0
    is.sig <- p.adjust(de$p.value[,coef], method = "BH") < 0.05
  }

  names(universe) <- row.names(de)

  # keep at least nmin up and down genes
  nsig.up <- sum(is.sig & is.up)
  nsig.dn <- sum(is.sig & !is.up)

  up <- head(universe[is.up], max(nmin, nsig.up))
  dn <- head(universe[!is.up], max(nmin, nsig.dn))

  de <- list(Up = up, Down = dn)

  universe <- universe[!duplicated(universe)]
  go_res <- run_go(de, universe, species, gs_dir, genego, gonames)
  return(go_res)
}


run_go <- function(de, universe, species, gs_dir, genego = NULL, gonames = NULL, min.genes = 4, max.genes = 250, min.diff = 3) {

  nsets <- length(de)

  # Ensure de gene IDs are unique and of character mode
  for (s in 1:nsets) de[[s]] <- unique(as.character(de[[s]]))

  # Restrict DE genes to universe
  for (s in 1:nsets) de[[s]] <- de[[s]][de[[s]] %in% universe]

  #	Check universe isn't empty
  NGenes <- length(universe)
  if (NGenes < 1L) stop("No annotated genes found in universe")

  # load gene to GO map and GO names
  if (is.null(genego)) genego <- get_genego(species, gs_dir)

  # Restrict pathways to universe
  genego <- genego[genego$gene_id %in% universe, ]
  if (is.null(gonames)) gonames <- get_gonames(genego, gs_dir)

  #	Overlap pathways with DE genes
  #	Create incidence matrix (X) of gene.pathway by DE sets
  X <- matrix(1, nrow(genego),nsets+1)
  colnames(X) <- c("N", names(de))
  for (s in 1:nsets) X[, s+1] <- (genego[[1]] %in% de[[s]])

  #	Count genes and DE genes for each pathway
  S <- rowsum(X, group=genego[[2]], reorder=FALSE)

  # get gene identities
  gnames <- names(universe)
  names(gnames) <- universe

  gids <- tibble::tibble(cbind(genego, X))
  gids <- gids %>%
    dplyr::filter(Up|Down) %>%
    dplyr::mutate(gene_name = gnames[gene_id]) %>%
    dplyr::group_by(go_id) %>%
    dplyr::summarise(Up.Genes = name_list(gene_name[Up > 0], go_id[1]),
                     Down.Genes = name_list(gene_name[Down > 0], go_id[1]))

  # Overlap tests
  PValue <- matrix(0, nrow=nrow(S), ncol=nsets)
  colnames(PValue) <- rep("P.Value", 2)
  nde <- lengths(de, use.names=FALSE)

  #	Fisher's exact test
  for (j in seq_len(nsets))
    PValue[,j] <- phyper(S[,1L+j]-0.5, nde[j], NGenes-nde[j], S[,"N"], lower.tail=FALSE)

  # Assemble output
  GOID <- rownames(S)
  TERM <- gonames[GOID]$TERM
  go_up <- data.frame(Term=TERM, Ont='BP', S, P.Value=PValue[, 1], stringsAsFactors=FALSE)
  go_dn <- data.frame(Term=TERM, Ont='BP', S, P.Value=PValue[, 2], stringsAsFactors=FALSE)

  go_up <- go_up[order(go_up$P.Value), ]
  go_dn <- go_dn[order(go_dn$P.Value), ]

  # limit to terms with min.genes
  go_up <- go_up[go_up$Up >= min.genes, ]
  go_dn <- go_dn[go_dn$Down >= min.genes, ]

  # limit to terms with lt max.genes
  go_up <- go_up[go_up$N < max.genes, ]
  go_dn <- go_dn[go_dn$N < max.genes, ]

  # limit to terms with gt min.diff per 10
  ndiff_up <- go_up$Up - go_up$Down
  ndiff_dn <- go_dn$Down - go_dn$Up
  nmin_up <- pmax(1, go_up$Up/10)*min.diff
  nmin_dn <- pmax(1, go_dn$Down/10)*min.diff

  go_up <- go_up[ndiff_up > nmin_up, ]
  go_dn <- go_dn[ndiff_dn > nmin_dn, ]

  # adjust
  go_up$FDR <- p.adjust(go_up$P.Value, 'BH')
  go_dn$FDR <- p.adjust(go_dn$P.Value, 'BH')

  # filter by cutoff/min.go
  nsig.up <- sum(go_up$FDR < 0.05)
  nsig.dn <- sum(go_dn$FDR < 0.05)
  go_up <- head(go_up, nsig.up)
  go_dn <- head(go_dn, nsig.dn)

  # filter terms both up and down
  both <- intersect(row.names(go_up), row.names(go_dn))
  go_up <- go_up[!row.names(go_up) %in% both, ]
  go_dn <- go_dn[!row.names(go_dn) %in% both, ]


  # for adding gene names
  up_genes <- gids %>%
    dplyr::filter(go_id %in% row.names(go_up)) %>%
    dplyr::mutate(idx = match(go_id, row.names(go_up))) %>%
    dplyr::arrange(idx) %>%
    dplyr::pull(Up.Genes)

  dn_genes <- gids %>%
    dplyr::filter(go_id %in% row.names(go_dn)) %>%
    dplyr::mutate(idx = match(go_id, row.names(go_dn))) %>%
    dplyr::arrange(idx) %>%
    dplyr::pull(Down.Genes)


  # simplify using similarity
  # also formats FDRs and adds gene names
  go_up <- simplify(go_up, up_genes) %>% dplyr::rename(Genes = genes)
  go_dn <- simplify(go_dn, dn_genes) %>% dplyr::rename(Genes = genes)

  return(list(up=go_up, dn=go_dn))
}


name_list <- function(genes, name) {
  res <- list(genes)
  names(res) <- name
  return(res)
}

jaccard <- function (regmat) {

  nart <- ncol(regmat)
  jdist <- rep(0, nart * nart)
  dim(jdist) <- c(nart, nart)
  reg.col.sum <- apply(regmat, 2, sum)
  reg.aggrement <- t(regmat) %*% regmat
  reg.aggrement/(reg.col.sum - t(t(reg.aggrement) - reg.col.sum))
}


simplify <- function(go_res, genes, cutoff = 0.7, by = 'FDR') {
  go_res$genes <- rep(NA, nrow(go_res))
  if (nrow(go_res) == 0) return(go_res)

  genes <- genes[row.names(go_res)]
  go_res$ID <- row.names(go_res)

  x <- table(stack(genes))
  sim <- jaccard(x)

  go1 <- go2 <- similarity <- NULL

  sim.df <- as.data.frame(sim)
  sim.df$go1 <- row.names(sim.df)
  sim.df <- tidyr::gather(sim.df, go2, similarity, -go1)

  sim.df <- sim.df[!is.na(sim.df$similarity),]

  ## feature 'by' is attached to 'go1'
  sim.df <- merge(sim.df, go_res[, c("ID", by)], by.x="go1", by.y="ID")
  sim.df$go2 <- as.character(sim.df$go2)

  ID <- go_res$ID
  go_res$group <- 0
  go_res$genes <- NA

  used <- character()
  for (i in seq_along(ID)) {
    ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff)
    if (length(ii) < 2) next

    sim_subset <- sim.df[ii,]

    jj <- which.min(sim_subset[, by])
    usedi <- sim_subset$go1
    usedi <- usedi[!usedi %in% used]

    if (length(usedi)) {
      go_res[usedi, 'group'] <- max(go_res$group)+1

      # put genes common to group first
      gni <- genes[usedi]
      common <- Reduce(intersect, gni)
      gni <- lapply(gni, function(x) c(common, setdiff(x, common)))
      go_res[usedi, 'genes'] <- unlist(lapply(gni, function(gs) paste(gs, collapse=', ')))

      used <- c(used, usedi)
    }
  }
  no.grp <- go_res$group == 0
  max.grp <- max(go_res$group)
  go_res$group[no.grp] <- seq(max.grp+1, max.grp+sum(no.grp))
  go_res$genes[no.grp] <- unlist(lapply(genes[no.grp], function(gs) paste(gs, collapse=', ')))

  go_res <- go_res %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(gmin = min(FDR)) %>%
    dplyr::arrange(FDR) %>%
    dplyr::mutate(Term = group_terms(Term)) %>%
    dplyr::mutate(FDR = format.pval(FDR, eps = 0.0001, digits = 1)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(gmin, group) %>%
    dplyr::select(-gmin, -group) %>%
    as.data.frame()

  row.names(go_res) <- go_res$ID
  go_res$ID <- NULL

  return(go_res)
}

group_terms <- function(terms) {
  nterms <- length(terms)
  if (nterms == 1) return(terms)
  if (nterms == 2) return(c(terms[1], paste0('      ° -- ', terms[2])))

  c(terms[1],
    paste0('      ¦ -- ', terms[2:(nterms-1)]),
    paste0('      ° -- ', terms[nterms]))
}


#' Get group levels for bulk data plots
#'
#' @param pdata Data.frame of phenotype data
#'
#'
#' @keywords internal
get_group_levels <- function(pdata) {
  group <- pdata$`Group name`
  group_order <- order(unique(pdata$Group))
  unique(group)[group_order]
}


#' Check if newly uploaded pdata is the same as previously uploaded
#'
#' @keywords internal
check_bulk_changed <- function(prev, pdata) {

  # don't re-run if Group name or Pairs the same
  changed <-
    !identical(prev$`Group name`, pdata$`Group name`) |
    !identical(prev$Pair, pdata$Pair)


  return(changed)
}

add_vsd <- function(eset, vsd_path = NULL, rna_seq = TRUE) {

  # for cases where manually added (e.g. nanostring dataset)
  els <- Biobase::assayDataElementNames(eset)
  if ('vsd' %in% els) return(eset)

  if (!rna_seq) {
    # for microarray use exprs
    vsd <- Biobase::assayDataElement(eset, 'exprs')

  } else if (!is.null(vsd_path) && file.exists(vsd_path)) {
    vsd <- qs::qread(vsd_path)

  } else {
    vsd <- crossmeta::get_vsd(eset)
    vsd <- SummarizedExperiment::assay(vsd)
    qs::qsave(vsd, vsd_path)
  }

  Biobase::assayDataElement(eset, 'vsd') <- vsd
  return(eset)
}

add_adjusted <- function(eset, adj_path, svobj = list(sv = NULL), numsv = 0) {

  if (file.exists(adj_path)) {
    Biobase::assayDataElement(eset, 'adjusted') <- qs::qread(adj_path)

  } else {
    eset <- crossmeta::add_adjusted(eset, svobj, numsv)
    adj <- Biobase::assayDataElement(eset, 'adjusted')
    qs::qsave(adj, adj_path)
  }

  return(eset)
}

iqr_replicates <- function(eset, keep_path, annot = "SYMBOL", rm.dup = FALSE) {

  if (file.exists(keep_path)) {
    keep <- qs::qread(keep_path)
    eset <- eset[keep, ]
    Biobase::featureNames(eset) <- Biobase::fData(eset)[, annot]

  } else {
    Biobase::fData(eset)$iqr_rows <- 1:nrow(eset)
    eset <- crossmeta::iqr_replicates(eset, annot, rm.dup)
    keep <- Biobase::fData(eset)$iqr_rows
    qs::qsave(keep, keep_path)
  }

  return(eset)
}
