
get_path_df <- function(path_id, anal) {

  # add dprimes and vardprime values
  anal <- add_es(anal)
  top_table <- anal$top_table

  path_enids <- gslist[[path_id]]
  path_genes <- names(path_enids)

  # subset top table to genes in the pathway
  top_table <- top_table[row.names(top_table) %in% path_genes, ]

  path_df <- data.frame(
    Gene = row.names(top_table),
    Dprime = top_table$dprime,
    sd = sqrt(top_table$vardprime),
    Link = paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", row.names(top_table), "'>", row.names(top_table), "</a>"), stringsAsFactors = FALSE
  )

  path_df <- path_df %>%
    arrange(desc(abs(Dprime))) %>%
    mutate(Gene = factor(Gene, levels = Gene))

  return(path_df)
}

pathPage <- function(input, output, session, new_anal, data_dir) {
  form <- callModule(pathForm, 'form',
                     new_anal = new_anal,
                     data_dir)



  # the gene plot
  pl <- reactive({

    diffs <- form$diffs()
    path_id <- form$pathway()
    anal <- diffs$anal

    req(path_id, anal)

    path_df <- get_path_df(path_id, anal)
    # 30 pixels width per gene in pathway
    plot_width <- nrow(path_df)*25 + 125


    plotly::plot_ly(data = path_df,
                    y = ~Dprime,
                    x = ~Gene,
                    text = ~Gene,
                    customdata = ~sd,
                    type = 'scatter',
                    mode = 'markers',
                    width = plot_width,
                    height = 550,
                    marker = list(size = 5, color = '#000000'),
                    error_y = ~list(array = sd, color = '#000000', thickness = 0.5, width = 0),
                    hovertemplate = paste0(
                      '<span style="color: crimson; font-weight: bold;">Gene</span>: %{text}<br>',
                      '<span style="color: crimson; font-weight: bold;">Dprime</span>: %{y:.2f}<br>',
                      '<span style="color: crimson; font-weight: bold;">SD</span>: %{customdata:.2f}',
                      '<extra></extra>')
    ) %>%
      plotly::config(displayModeBar = FALSE) %>%
      plotly::layout(yaxis = list(fixedrange = TRUE),
                     xaxis = list(fixedrange = TRUE,
                                  range = c(-2, nrow(path_df) + 1),
                                  tickmode = 'array',
                                  tickvals = 0:nrow(path_df),
                                  ticktext = ~Link,
                                  tickangle = -45),
                     autosize = FALSE)

  })

  # inside observe to allow dynamic width
  output$path_plot <- plotly::renderPlotly({
    pl()
  })


}

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

  if (with_type) anals$type <- 'Single Cell'

  return(anals)
}

scPathClusters <- function(input, output, session, data_dir, anal, is_sc) {

  contrast_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))

  # differential and pathway analyses
  diffs <- reactiveVal()

  sc_dir <- reactive({
    req(is_sc())
    file.path(data_dir, 'single-cell')
  })

  scseq <- reactive({
    scseq_path <- scseq_part_path(sc_dir(), anal()$anal_name, 'scseq')
    readRDS(scseq_path)
  })

  # TODO: update if annotation change from Single Cell tab
  annot <- reactive({
    annot_path <- scseq_part_path(sc_dir(), anal()$anal_name, 'annot')
    readRDS(annot_path)
  })


  cluster_choices <- reactive({
    # value is original cluster number so that saved pathway analysis name
    # isn't affected by updates to cluster annotation
    scseq <- scseq()
    value <- levels(Seurat::Idents(scseq))
    get_cluster_choices(clusters = annot(), scseq = scseq, value = value)
  })




  # update UI for contrast/cluster choices
  observeEvent(cluster_choices(), {
    updateSelectizeInput(session, 'selected_clusters',
                         choices = cluster_choices(),
                         options = contrast_options, server = TRUE)
  })

  observeEvent(input$run_comparison, {

    # run differential expression analysis ----
    # set idents to ctrl and test
    scseq <- scseq()
    clusters <- input$selected_clusters
    anal <- diff_expr_sc(scseq, clusters = clusters)


    # run pathway analysis
    panal <- diff_path_sc(scseq,
                          prev_anal = anal,
                          data_dir = sc_dir(),
                          anal_name = anal()$anal_name,
                          clusters = clusters)

    diffs(list(anal = anal, path = panal))
  })

  return(list(
    diffs = diffs
  ))


}

diff_path_sc <- function(scseq, prev_anal, data_dir, anal_name, clusters) {

  # load previous if exists
  clusters_name <- paste(clusters, collapse = ',')
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


diff_expr_sc <- function(scseq, clusters, exclude_ambient = TRUE) {
  Seurat::Idents(scseq) <- scseq$orig.ident

  # exclude non-selected clusters
  scseq <-  scseq[, scseq$seurat_clusters %in% clusters]
  ebayes_sv <- fit_ebayes_scseq(scseq, ident.1 = 'test', ident.2 = 'ctrl')
  tt <- limma::topTable(ebayes_sv, coef = 1, number = Inf)


  if (exclude_ambient) {
    ambient <- get_ambient(scseq, tt)

    scseq <- scseq[!row.names(scseq) %in% ambient, ]
    ebayes_sv <- fit_ebayes_scseq(scseq, ident.1 = 'test', ident.2 = 'ctrl')
    tt <- limma::topTable(ebayes_sv, coef = 1, number = Inf)
  }

  # need ebayes_sv and pdata for add_es to get dprimes and vardprimes
  pdata <- scseq[['orig.ident']]
  pdata$group <- as.character(pdata$orig.ident)
  pdata$orig.ident <- NULL

  anal <- list(top_table = tt,
               ebayes_sv = ebayes_sv,
               pdata = pdata)

  return(anal)
}

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

get_ambient <- function(scseq, markers) {

  fts <- scseq[['SCT']]@meta.features
  test.ambient <- row.names(fts)[fts$test_ambient]
  ctrl.ambient <- row.names(fts)[fts$ctrl_ambient]

  # exclude test/ctrl ambient if positive/negative effect size
  # opposite would decrease extent of gene expression difference but not direction
  pos.test <- markers[test.ambient, 't'] > 0
  neg.ctrl <- markers[ctrl.ambient, 't'] < 0

  test.exclude <- test.ambient[pos.test]
  ctrl.exclude <- ctrl.ambient[neg.ctrl]

  ambient <- c(test.exclude, ctrl.exclude)

  return(unique(ambient))
}

pathForm <- function(input, output, session, new_anal, data_dir) {

  # reload analysis choices if new analysis
  anals <- reactive({
    new_anal()
    bulk_anals <- load_bulk_anals(data_dir, with_type = TRUE)
    scseq_anals <- load_scseq_anals(data_dir, with_type = TRUE)

    anals <- rbind(bulk_anals, scseq_anals)
    anals$value <- seq_len(nrow(anals))
    return(anals)
  })


  # update analysis choices
  observe({
    anals <- anals()
    req(anals)
    updateSelectizeInput(session, 'anal', choices = rbind(rep(NA, 5), anals), server = TRUE)
  })

  # get directory/name info about analysis
  anal <- reactive({
    row_num <- input$anal
    anals <- anals()
    req(row_num, anals)

    anals[row_num, ]
  })

  # show/hide single cell stuff
  is_sc <- reactive({
    anal <- anal()
    anal$type == 'Single Cell'
  })

  observe({
    shinyjs::toggle('sc_clusters_container', condition = is_sc())
  })

  # get single-cell data
  sc_inputs <- callModule(scPathClusters, 'sc_clusters',
                          data_dir = data_dir,
                          anal = anal,
                          is_sc = is_sc)



  # file paths to pathway and analysis results
  fpaths <- reactive({
    anal <- anal()
    is_sc <- is_sc()

    if (is_sc()) {
      return(NULL)
    }

    dataset_dir <-  file.path(data_dir, 'bulk', anal$dataset_dir)
    anal_name <- anal$anal_name

    list(
      diff_anal = file.path(dataset_dir, paste0('diff_expr_symbol_', anal_name, '.rds')),
      diff_path = file.path(dataset_dir, paste0('diff_path_', anal_name, '.rds'))
    )
  })

  # load pathway and analysis results
  diffs <- reactive({
    fpaths <- fpaths()
    if (is_sc()) return (sc_inputs$diffs())


    list(
      path = readRDS(fpaths$diff_path),
      anal = readRDS(fpaths$diff_anal)
    )
  })

  show_up <- reactive({
    input$show_up %% 2 == 0
  })

  observe({
    name <- ifelse(show_up(),  'chevron-up', 'chevron-down')
    updateActionButton(session, 'show_up', icon = icon(name, 'fa-fw'))
  })

  path_directions <- reactive({
    diffs <- diffs()
    res <- diffs$path$res

    if (is.null(res)) return(NULL)

    path_directions <- get_path_directions(diffs$anal$top_table)
    path_directions[res$ID, ]
  })

  path_choices <- reactive({
    diffs <- diffs()
    res <- diffs$path$res
    directions <- path_directions()

    if (is.null(res)) return(NULL)
    is.up <- directions$is.up

    path_choices <- data.frame(
      name = res$Name,
      value = res$ID,
      label = res$Name,
      direction_label = directions$label,
      is.up = is.up,
      fdr = format.pval(res$Ppadog, eps = 0.001, digits = 2),
      stringsAsFactors = FALSE)

    filter <- if (show_up()) is.up else !is.up

    return(path_choices[filter, ])
  })

  # update pathway dropdown
  observe({
    updateSelectizeInput(session, 'pathway',
                         choices = path_choices(),
                         options = list(render= I('{option: pathOptions}')),
                         server = TRUE)
  })


  # open KEGG when click button
  observeEvent(input$kegg, {
    path_id <- input$pathway
    req(path_id)
    kegg_link <- paste0('https://www.genome.jp/kegg-bin/show_pathway?map', path_id)
    runjs(paste0("window.open('", kegg_link, "')"))
  })

  return(list(
    diffs = diffs,
    pathway = reactive(input$pathway)
  ))
}

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

get_path_choices <- function(res) {
  data.frame(
    name = res$Name,
    value = res$ID,
    fdr = res$Ppadog
  )
}

server <- function(input, output, session) {

  # get arguments from calling function
  # defaults for server
  # base directory contains data_dir folder
  data_dir <- getShinyOption('data_dir', '/srv/shiny-server/drugseqr/data_dir')

  # for testing don't seem to be able to pass arguments as options
  if (isTRUE(getOption('shiny.testmode'))) {

    # reset data for testing
    data_dir <- 'tests/data/test'
    static_dir <- 'tests/data/static'
    unlink(data_dir, recursive = TRUE)
    dir.create(data_dir)
    file.copy(list.files(static_dir, full.names = TRUE), data_dir, recursive = TRUE)
  }

  sc_dir <- file.path(data_dir, 'single-cell')
  bulk_dir <- file.path(data_dir, 'bulk')

  dir.create(sc_dir, showWarnings = FALSE)
  dir.create(bulk_dir, showWarnings = FALSE)


  # single cell analysis and options
  scPage <- callModule(scPage, 'sc',
                       sc_dir = sc_dir)

  dsPage <- callModule(dsPage, 'datasets',
                       data_dir = data_dir)


  drugsPage <- callModule(drugsPage, 'drug',
                          new_anal = dsPage$new_anal,
                          data_dir = data_dir)

  pathPage <- callModule(pathPage, 'pathways',
                         new_anal = dsPage$new_anal,
                         data_dir = data_dir)



}
