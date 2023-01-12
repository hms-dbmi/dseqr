#' Logic for Bulk Data Tab
#'
#' to be called with \link[shiny]{callModule}
#'
#' @param input,output,session standard shiny module boilerplate
#' @param project_dir path to folder with project files
#' @param sc_dir sub folder of \code{project_dir} where single-cell data is stored
#' @param bulk_dir sub folder of \code{project_dir} where bulk data is stored
#' @param add_bulk reactive that triggers modal to upload a bulk dataset
#' @param remove_bulk reactive that triggers modal for deleting bulk datasets
#' @inheritParams run_dseqr
#'
#' @return list with reactive \code{new_dataset} that is triggered with a new
#'   bulk dataset is added.
#'
#' @export
bulkPage <- function(input, output, session, project_dir, sc_dir, bulk_dir, tx2gene_dir, indices_dir, gs_dir, add_bulk, remove_bulk) {

  msg_quant <- reactiveVal()

  eset <- reactive(qread.safe(file.path(bulkForm$dataset_dir(), 'eset.qs')))


  explore_eset <- exploreEset(eset = eset,
                              dataset_dir = bulkForm$dataset_dir,
                              explore_pdata = dsExploreTable$pdata,
                              numsv = bulkForm$numsv_r,
                              svobj = bulkForm$svobj_r)


  bulkForm <- callModule(bulkForm, 'form',
                         project_dir = project_dir,
                         sc_dir = sc_dir,
                         bulk_dir = bulk_dir,
                         tx2gene_dir = tx2gene_dir,
                         gs_dir = gs_dir,
                         msg_quant = msg_quant,
                         explore_eset = explore_eset,
                         indices_dir = indices_dir,
                         add_bulk = add_bulk,
                         remove_bulk = remove_bulk)

  # mds plotly with different orientations for mobile/desktop

  bulkMDS <- callModule(bulkMDS, 'bulk_mds',
                        explore_eset = explore_eset,
                        dataset_name = bulkForm$dataset_name,
                        numsv = bulkForm$numsv_r,
                        bulk_dir = bulk_dir)

  callModule(bulkMDSplotly, 'mds_plotly_unadjusted',
             explore_eset = explore_eset,
             dataset_name = bulkForm$dataset_name,
             numsv = bulkForm$numsv_r,
             mds = bulkMDS$mds,
             group_colors = bulkMDS$group_colors,
             adjusted = FALSE)

  callModule(bulkMDSplotly, 'mds_plotly_adjusted',
             explore_eset = explore_eset,
             dataset_name = bulkForm$dataset_name,
             numsv = bulkForm$numsv_r,
             mds = bulkMDS$mds,
             group_colors = bulkMDS$group_colors,
             adjusted = TRUE)


  # toggle tables
  observe({
    toggle('anal_table_container', condition = bulkForm$show_anal())
  })

  # toggle plots
  sel_genes <- reactive(length(bulkForm$explore_genes() > 0))

  observe({
    toggle('mds_plotly_container', condition = !sel_genes() & !bulkForm$show_deconv())
    toggle('gene_plotly_container', condition = sel_genes() & !bulkForm$show_deconv())
    toggle('cells_plotly_container', condition = bulkForm$show_deconv())
  })


  up_annot <- callModule(bulkAnnot, 'anal',
                         dataset_name = bulkForm$dataset_name,
                         annot = dsExploreTable$pdata)


  dsExploreTable <- callModule(bulkExploreTable, 'explore',
                               eset = eset,
                               up_annot = up_annot,
                               project_dir = project_dir,
                               dataset_dir = bulkForm$dataset_dir,
                               dataset_name = bulkForm$dataset_name,
                               svobj_r = bulkForm$svobj_r,
                               numsv_r = bulkForm$numsv_r)


  callModule(bulkGenePlotly, 'gene_plotly',
             eset = explore_eset,
             explore_genes = bulkForm$explore_genes,
             dataset_name = bulkForm$dataset_name)

  callModule(bulkCellsPlotly, 'cells_plotly',
             est_prop = bulkForm$est_prop,
             pdata = dsExploreTable$pdata,
             dataset_name = bulkForm$dataset_name,
             contrast_groups = bulkForm$contrast_groups)




  return(list(
    new_dataset = bulkForm$new_dataset
  ))
}


#' Logic for Bulk Data MDS data
#'
#' @keywords internal
#' @noRd
bulkMDS <- function(input, output, session, explore_eset, dataset_name, numsv, bulk_dir) {

  group <- reactive({
    eset <- explore_eset()
    if (is.null(eset)) return(NULL)
    pdata <- Biobase::pData(eset)

    # setup group factor and colors
    group <- pdata$`Group name`
    group_order <- order(unique(pdata$Group))
    group_levels <- unique(group)[group_order]
    factor(group, levels = group_levels)
  })


  group_colors <- reactive({
    get_palette(levels(group()))
  })

  mds_path <- reactive({
    file.path(bulk_dir(), dataset_name(), paste0('mds_', numsv(), 'svs.qs'))
  })


  mds <- reactive({
    eset <- explore_eset()
    req(eset)

    mds_path <- isolate(mds_path())
    if (file.exists(mds_path)) {
      mds <- qs::qread(mds_path)

    } else {
      progress <- Progress$new(session, min=0, max = 2)
      on.exit(progress$close())
      progress$set(message = "Scaling for MDS plots", value = 1)

      vsd <- Biobase::assayDataElement(eset, 'vsd')
      adj <- Biobase::assayDataElement(eset, 'adjusted')
      mds <- get_mds(vsd, adj, group())
      progress$set(value = 2)

      qs::qsave(mds, mds_path)
    }

    return(mds)
  })

  return(list(
    mds = mds,
    group_colors = group_colors
  ))
}


#' Logic for Bulk Data MDS plotly
#'
#' @keywords internal
#' @noRd
bulkMDSplotly <- function(input, output, session, explore_eset, dataset_name, numsv, mds, group_colors, adjusted) {

  # MDS plot
  plot <- reactive({
    mds <- mds()
    plotlyMDS(mds$scaling, mds$scaling_adj, group_colors = group_colors(), adjusted = adjusted)
  })

  base_fname <- reactive(paste0(dataset_name(), '_', numsv(), 'SV'))

  content <- function(file) {
    #go to a temp dir to avoid permission issues
    owd <- setwd(tempdir())
    on.exit(setwd(owd))

    base_fname <- base_fname()
    eset <- explore_eset()
    vsd <- Biobase::assayDataElement(eset, 'vsd')
    pdata <- Biobase::pData(eset)

    vsd_fname <- paste0(base_fname, '.csv')
    pdata_fname <- paste0(base_fname, '_pdata.csv')
    utils::write.csv(vsd, vsd_fname)
    utils::write.csv(pdata, pdata_fname)

    #create the zip file
    utils::zip(file, c(vsd_fname, pdata_fname))
  }

  filename <- function() {paste0(base_fname(), '_', Sys.Date(), '.zip')}


  # not currently downloadable so use default fname_fun and data_fun
  callModule(shinydlplot::downloadablePlotly, 'plotly',
             plot = plot,
             content = content,
             filename = filename,
             title = 'Download adjusted data')
}


#' Logic for Bulk Data gene plotly
#'
#' @keywords internal
#' @noRd
bulkGenePlotly <- function(input, output, session, eset, explore_genes, dataset_name) {

  boxplotly_args <- reactive({
    # need eset and at least one gene
    explore_genes <- explore_genes()
    req(explore_genes)
    eset <- eset()

    # prevent warning when switching between dataset (eset updated before genes)
    req(all(explore_genes %in% row.names(eset)))
    get_boxplotly_gene_args(eset, explore_genes, dataset_name())
  })

  plot <- reactive({

    args <- boxplotly_args()

    boxPlotly(df = args$df,
              boxgap = args$boxgap,
              boxgroupgap = args$boxgroupgap,
              plot_fname = args$plot_fname,
              ytitle = 'Normalized Expression',
              xtitle = 'Gene')
  })

  filename <- function() {
    paste(dataset_name(), '_', paste(explore_genes(), collapse = '_'), '_', Sys.Date(), ".csv", sep = "")
  }


  content <- function(file) {
    args <- boxplotly_args()
    df <- args$df
    df$color <- NULL
    colnames(df) <- c('Sample', 'Gene', 'Normalized Expression', 'Group')
    utils::write.csv(df, file, row.names = FALSE)
  }

  callModule(shinydlplot::downloadablePlotly,
             'plotly',
             plot = plot,
             filename = filename,
             content = content)



}


#' Logic for Bulk Data cell type deconvolution plotly
#'
#' @keywords internal
#' @noRd
bulkCellsPlotly <- function(input, output, session, est_prop, pdata, dataset_name, contrast_groups) {

  boxplotly_args <- reactive({
    # need at least two groups

    pdata <- pdata()
    pdata <- pdata[!is.na(pdata$Group), ]
    req(length(unique(pdata$Group)) > 1)

    est_prop <- est_prop()
    req(est_prop)

    contrast <- contrast_groups()

    get_boxplotly_cell_args(pdata, est_prop, dataset_name(), contrast)
  })

  plot <- reactive({

    args <- boxplotly_args()
    contrast <- contrast_groups()
    is.comparison <- length(contrast) == 2
    ytitle <- 'Estimated Proportion'
    if (is.comparison) ytitle <- paste(ytitle, '(FDR)')

    boxPlotlyCells(df = args$df,
                   boxgap = args$boxgap,
                   boxgroupgap = args$boxgroupgap,
                   pvals = args$pvals,
                   plot_fname = args$plot_fname,
                   ytitle = ytitle,
                   xtitle = 'Cluster')
  })

  filename <- function() {
    clusters <- colnames(est_prop())
    clusters <- gsub(' ', '', clusters)
    paste(dataset_name(), '_', paste(clusters, collapse = '_'), '_', Sys.Date(), ".csv", sep = "")
  }


  content <- function(file) {
    args <- boxplotly_args()
    df <- args$df
    df$color <- df$pair <- NULL
    colnames(df) <- c('Sample', 'Group', 'Proportion', 'Cluster')
    utils::write.csv(df, file, row.names = FALSE)
  }

  callModule(shinydlplot::downloadablePlotly, 'plotly', plot = plot, filename = filename, content = content)


}


#' Logic for Bulk Data form
#'
#' @keywords internal
#' @noRd
bulkForm <- function(input, output, session, project_dir, sc_dir, bulk_dir, tx2gene_dir, msg_quant, explore_eset, indices_dir, gs_dir, add_bulk, remove_bulk) {

  dataset <- callModule(bulkDataset, 'selected_dataset',
                        project_dir = project_dir,
                        sc_dir = sc_dir,
                        bulk_dir = bulk_dir,
                        tx2gene_dir = tx2gene_dir,
                        indices_dir = indices_dir,
                        explore_eset = explore_eset,
                        add_bulk = add_bulk,
                        remove_bulk = remove_bulk)


  show_anal <- reactive(dataset$dataset_name() != '')

  observe(toggle('anal_dataset_panel', condition = show_anal()))


  anal <- callModule(bulkFormAnal, 'anal_form',
                     project_dir = project_dir,
                     dataset_dir = dataset$dataset_dir,
                     gs_dir = gs_dir,
                     dataset_name = dataset$dataset_name,
                     explore_eset = explore_eset,
                     numsv_r = dataset$numsv_r,
                     svobj_r = dataset$svobj_r)



  return(list(
    fastq_dir = dataset$fastq_dir,
    anal_labels = anal$labels,
    explore_genes = anal$explore_genes,
    dataset_name = dataset$dataset_name,
    dataset_dir = dataset$dataset_dir,
    show_anal = show_anal,
    is.explore = anal$is.explore,
    est_prop = dataset$est_prop,
    show_deconv = dataset$show_deconv,
    numsv_r = dataset$numsv_r,
    svobj_r = dataset$svobj_r,
    contrast_groups = anal$contrast_groups
  ))

}


#' Logic for selected dataset part of bulkFrom
#'
#' @keywords internal
#' @noRd
bulkDataset <- function(input, output, session, sc_dir, bulk_dir, tx2gene_dir, project_dir, indices_dir, explore_eset, add_bulk, remove_bulk) {

  new_dataset <- reactiveVal()
  options <- list(optgroupField = 'type',
                  render = I('{option: bulkDatasetOptions, item: bulkDatasetItem}'))



  # only show nsv/deconv toggle if dataset with groups
  have_groups <- reactive(isTruthy(explore_eset()))

  observe({
    shinypanel::toggleSelectizeButtons('dataset_name',
                                       button_ids = c('show_nsv', 'show_deconv'),
                                       condition = have_groups())
  })


  datasets <- reactive({
    new_dataset()
    load_bulk_datasets(project_dir())
  })


  observe({
    datasets <- datasets()
    req(datasets)
    datasets <- datasets_to_list(datasets)
    updateSelectizeInput(session, 'dataset_name', selected = isolate(input$dataset_name), choices = c("", datasets), options = options)
  })


  # Delete Bulk Dataset
  # ---

  # show remove bulk datasets modal
  observeEvent(remove_bulk(), {
    ds <- datasets()
    ds <- tibble::as_tibble(ds)

    choices <- ds$dataset_name
    showModal(deleteModal(session, choices, type = 'Bulk'))
  })

  allow_delete <- reactive(isTruthy(input$remove_datasets) & input$confirm_delete == 'delete')

  observe({
    shinyjs::toggleState('delete_dataset', condition = allow_delete())
    shinyjs::toggleClass('delete_dataset', class = 'btn-danger', condition = allow_delete())
  })

  observe({
    shinyjs::toggle('confirm_delete_container', condition = isTruthy(input$remove_datasets))
  })

  observeEvent(input$delete_dataset, {
    remove_datasets <- input$remove_datasets
    unlink(file.path(bulk_dir(), remove_datasets), recursive = TRUE)
    updateTextInput(session, 'confirm_delete', value = '')
    removeModal()
    new_dataset(paste0(remove_datasets, '_deleted'))
  })


  # Add Bulk Dataset
  # ---

  uploads_table <- reactiveVal()
  pairs <- reactiveVal()
  reps <- reactiveVal()

  error_msg_bulk_file <- reactiveVal()
  error_msg_labels <- reactiveVal()


  # initial uploads DT
  # updates in order to remove/add Pair column
  empty_table <- reactiveVal()

  observe({
    prev <- empty_table()
    df <- isolate(uploads_table_html())
    is_eset <- is_eset()

    if (is.null(df)) {
      df <- data.frame(
        ' ' = character(0),
        Pair = character(0),
        Replicate = character(0),
        File = character(0),
        Size = character(0),
        check.names = FALSE)
    }

    prev_paired <- !is.null(prev$Pair)
    curr_paired <- detected_paired()
    prev_eset <- !'Replicate' %in% colnames(prev)
    if (!curr_paired) df$Pair <- NULL
    if (is_eset) df$Pair <- df$Replicate <- NULL

    if (is.null(prev) || prev_paired != curr_paired || prev_eset != is_eset) empty_table(df)
  })


  # render upload DT
  output$uploads_table <- DT::renderDataTable({

    df <- empty_table()
    if (is.null(df)) return(NULL)
    targets <- 1
    if ('Pair' %in% colnames(df)) targets <- c(1, 2)
    if (!'Replicate' %in% colnames(df)) targets <- NULL
    DT::datatable(df,
                  class = 'cell-border dt-fake-height',
                  rownames = FALSE,
                  escape = FALSE, # to allow HTML in table
                  selection = 'multiple',
                  options = list(
                    columnDefs = list(list(className = 'dt-nopad', targets = targets)),
                    scrollX = TRUE,
                    dom = 't',
                    paging = FALSE
                  )) %>%
      DT::formatStyle('Size', `text-align` = 'right') %>%
      DT::formatStyle(c('File', 'Size'), color = 'gray')
  })

  # update rendered uploads table using proxy
  proxy <- DT::dataTableProxy('uploads_table')

  observe({
    table <- uploads_table_html()
    if (!detected_paired()) table$Pair <- NULL
    if (is_eset()) table$Pair <- table$Replicate <- NULL
    DT::replaceData(proxy, table, rownames = FALSE)
  })


  # open add dataset modal
  observeEvent(add_bulk(), {



    showModal(uploadBulkModal(
      session,
      have_uploads(),
      is_eset(),
      input$import_dataset_name,
      detected_paired()
    ))
  })


  # set message if tried to upload anything but specified
  observeEvent(input$up_raw_errors, {
    msg <- 'Only specified file types can be uploaded.'
    error_msg_bulk_file(msg)
  })

  # validate uploaded files
  observeEvent(uploads_table(), {
    up_df <- uploads_table()

    msg <- validate_bulk_uploads(up_df, is_eset())
    error_msg_bulk_file(msg)
  })

  # show any errors with uploads
  observe({
    msg <- error_msg_bulk_file()
    html('error_msg_bulk_file', html = msg)
    shinyjs::toggleClass('validate-up-fastq', 'has-error', condition = isTruthy(msg))
  })

  is_eset <- reactiveVal(FALSE)

  # append to uploads table
  # also add placeholder for pairs and replicates
  observeEvent(input$up_raw, {
    prev <- uploads_table()
    new <- input$up_raw
    new <- new[file.exists(new$datapath), ]

    is.eset <- grep('[.]qs$|[.]rds$', new$name)

    if (length(is.eset)) {
      # take first eset
      is_eset(TRUE)
      reps(NULL)
      pairs(NULL)
      uploads_table(new[is.eset[1], ])
    } else {
      is_eset(FALSE)
      # extend reps/pairs
      placeholder <- rep(NA, nrow(new))
      reps(c(reps(), placeholder))
      pairs(c(pairs(), placeholder))

      uploads_table(rbind.data.frame(prev, new))

    }

  })

  # handle delete row button
  observeEvent(input$delete_row, {
    selected_row <- as.numeric(strsplit(input$delete_row, "_")[[1]][3])
    reps <- reps()
    pairs <- pairs()
    uploads <- uploads_table()

    unlink(uploads$datapath[selected_row])

    # clear any associated replicate/apir
    rep <- reps[selected_row]
    pair <- pairs[selected_row]
    reps[reps == rep] <- NA
    pairs[pairs == pair] <- NA

    # remove selected row
    reps <- reps[-selected_row]
    pairs <- pairs[-selected_row]
    uploads <- uploads[-selected_row, ]

    if (!nrow(uploads))
      reps <- pairs <- uploads <- NULL

    reps(reps)
    pairs(pairs)
    uploads_table(uploads)
  })

  # update/initialize uploads table with html
  uploads_table_html <- reactiveVal()

  observe({
    df <- uploads_table()
    is_eset <- is_eset()
    if (is.null(df)) {
      uploads_table_html(NULL)
      return()
    }

    df <- df[, c('name', 'size')]
    df$size <- sapply(df$size, utils:::format.object_size, units = 'auto')
    colnames(df) <- c('File', 'Size')

    if (!is_eset)
      df <- dplyr::mutate(df, ' ' = NA, Pair = NA, Replicate = NA,  .before = 1)
    else
      df <- dplyr::mutate(df, ' ' = NA,  .before = 1)

    df$` ` <- getDeleteRowButtons(session, nrow(df))

    uploads_table_html(df)
  })


  # auto detected if paired
  detected_paired <- reactive({
    up_df <- uploads_table()
    if (is.null(up_df)) return(FALSE)
    if (is_eset()) return(FALSE)

    fastqs <- up_df$datapath
    fastqs <- fastqs[file.exists(fastqs)]

    # auto-detect if paired
    fastq_id1s <- rkal::get_fastq_id1s(fastqs)
    if (!length(fastq_id1s)) return(FALSE)

    pairs <- rkal::get_fastq_pairs(fastq_id1s)
    uniq.pairs <- unique(pairs)
    paired <- setequal(c("1", "2"), uniq.pairs) || uniq.pairs != '1'
    return(paired)
  })



  # dashed lines for 'Pair' html
  background <- 'url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAPklEQVQoU43Myw0AIAgEUbdAq7VADCQaPyww55dBKyQiHZkzBIwQLqQzCk9E4Ytc6KEPMnTBCG2YIYMVpHAC84EnVbOkv3wAAAAASUVORK5CYII=) repeat'

  # update Pair/Replicate html
  observe({

    # things that trigger update
    pdata <- uploads_table_html()
    req(pdata)

    reps <- reps()
    pairs <- pairs()

    # clear previous html
    pdata$Replicate <- pdata$Pair <- NA

    # update pdata Replicate column
    rep_nums <- sort(unique(setdiff(reps, NA)))
    rep_nums <- as.numeric(rep_nums)
    rep_colors <- get_palette(reps)

    for (rep_num in rep_nums) {
      color <- rep_colors[rep_num]
      rows <- which(reps == rep_num)
      pdata[rows, 'Replicate'] <- paste('<div style="background-color:', color, ';"></div>')
    }

    # update pdata Pair column
    if (detected_paired()) {
      pair_nums <- sort(unique(setdiff(pairs, NA)))
      pair_nums <- as.numeric(pair_nums)
      pair_colors <- get_palette(pairs)

      for (pair_num in pair_nums) {
        color <- pair_colors[pair_num]
        rows <- which(pairs == pair_num)
        pdata[rows, 'Pair'] <- paste('<div style="background:', color, background, ';"></div>')
      }
    }

    uploads_table_html(pdata)
  })


  # click 'Paired'
  shiny::observeEvent(input$pair, {

    reps <- reps()
    pairs <- pairs()

    # get rows
    rows  <- input$uploads_table_rows_selected

    # check for incomplete/wrong input
    msg <- rkal::validate_pairs(pairs, rows, reps)
    error_msg_labels(msg)

    if (is.null(msg)) {
      # add rows as a pair
      pairs[rows] <- get_next_number(pairs)
      pairs(pairs)
    }
  })

  get_next_number <- function(x) {
    numbers <- unique(setdiff(x, NA))
    if (!length(numbers)) return(1)

    setdiff(seq_len(max(numbers)+1), numbers)[1]
  }

  # click 'Replicate'
  shiny::observeEvent(input$rep, {

    reps <- reps()
    pairs <- pairs()

    # get rows
    rows  <- input$uploads_table_rows_selected
    msg <- rkal::validate_reps(pairs, rows, reps)
    error_msg_labels(msg)

    if (is.null(msg)) {
      # add rows as replicates
      reps[rows] <- get_next_number(reps)
      reps(reps)
    }
  })

  # click 'Reset'
  shiny::observeEvent(input$reset, {

    rows  <- input$uploads_table_rows_selected

    new_reps <- reps()
    new_pairs <- pairs()
    new_reps[rows] <- NA
    new_pairs[rows] <- NA

    reps(new_reps)
    pairs(new_pairs)
  })


  # show any errors with Pair/Replicate labels
  observe({
    msg <- error_msg_labels()
    html('error_msg_labels', html = msg)
    shinyjs::toggleClass('fastq_labels_container', 'has-error', condition = isTruthy(msg))
  })

  # hide/show inputs based on user selections
  have_uploads <- reactive(!is.null(uploads_table()))

  observe({
    shinyjs::toggleCssClass('fastq_labels-1', class = 'hidden', condition = !detected_paired())
    shinyjs::toggleCssClass('rep', class = 'radius-left', condition = !detected_paired())
  })

  observe({
    shinyjs::toggleCssClass('uploads_table_container', 'invisible-height', condition = !have_uploads())
  })

  observe({
    shinyjs::toggle('fastq_labels_container', condition = have_uploads() & !is_eset())
  })

  observe({
    shinyjs::toggle('import_name_container', condition = have_uploads())
  })

  observe({
    shinyjs::toggleState('import_bulk_dataset', condition = have_uploads())
  })


  # Run Import
  # ---


  # show missing name error
  error_msg_name <- reactiveVal()

  observe({
    msg <- error_msg_name()
    html('error_msg_name', html = msg)
    toggleClass('validate-up-name', 'has-error', condition = isTruthy(msg))
  })


  # validate that can import
  observeEvent(input$import_bulk_dataset, {

    reps <- reps()
    pairs <- pairs()
    up_df <- uploads_table()
    paired <- detected_paired()
    import_name <- input$import_dataset_name
    is_eset <- is_eset()

    # need a name

    msg_labels <- NULL
    if (!is_eset) msg_labels <- validate_bulk_labels(up_df, reps, pairs, paired)
    msg_bulk_uploads <- validate_bulk_uploads(up_df, is_eset)

    msg_name <- validate_bulk_name(import_name)

    error_msg_name(msg_name)
    error_msg_labels(msg_labels)
    error_msg_bulk_file(msg_bulk_uploads)

    if (is.null(msg_name) && is.null(msg_labels) && is.null(msg_bulk_uploads)) {

      removeModal()
      Sys.sleep(1)

      if (!is_eset) {
        showModal(confirmImportBulkModal(session))

      } else {
        progress <- Progress$new(session, min=0, max = 2)
        on.exit(progress$close())
        progress$set(message = "Importing ExpressionSet", value = 1)
        import_bulk_eset(up_df, import_name, bulk_dir())
        new_dataset(import_name)
        uploads_table(NULL)
        progress$set(value = 2)
      }
    }
  })

  import_bulk_eset <- function(up_df, import_name, bulk_dir) {
    data_dir <- file.path(bulk_dir, import_name)
    unlink(data_dir, recursive = TRUE)
    dir.create(data_dir)

    new_fpath <- file.path(data_dir, 'eset.qs')
    rds_upload <- grepl('[.]rds$', up_df$name)

    if (rds_upload) {
      eset <- readRDS(up_df$datapath)
      qs::qsave(eset, new_fpath)

    } else {
      file.move(up_df$datapath, new_fpath)
    }
  }



  # pdata used for quantification
  pdata <- reactive({
    up_df <- uploads_table()
    paired <- detected_paired()

    pdata <- tibble::tibble(
      'Pair' = pairs(),
      'Replicate' = reps(),
      'File Name' = up_df$name
    )

    return(pdata)
  })


  # run quantification upon confirmation
  quants <- reactiveValues()
  pquants <- reactiveValues()

  observeEvent(input$confirm, {

    removeModal()

    # check for index
    index_path <- rkal::get_kallisto_index(indices_dir)
    if(!length(index_path)) {
      shinyjs::alert('No kallisto index. See github README.')
      return(NULL)
    }

    # setup
    pdata <- pdata()
    paired <- detected_paired()
    dataset_name <- input$import_dataset_name
    fastq_dir <- file.path(bulk_dir(), dataset_name)

    # move files
    up_df <- uploads_table()
    unlink(fastq_dir, recursive = TRUE)
    dir.create(fastq_dir, showWarnings = FALSE)

    for (i in seq_len(nrow(up_df))) {
      dpath <- up_df$datapath[i]
      fpath <- file.path(fastq_dir, up_df$name[i])
      file.move(from = dpath, to = fpath)
    }

    quants[[dataset_name]] <- callr::r_bg(
      func = run_kallisto_bulk_bg,
      package = 'dseqr',
      args = list(
        indices_dir = indices_dir,
        data_dir = fastq_dir,
        quant_meta = pdata,
        paired = paired
      )
    )

    # Create a Progress object
    progress <- Progress$new(session, min=0, max = nrow(pdata)+2)
    progress$set(message = "Running pseudoalignment", value = 1)
    pquants[[dataset_name]] <- progress

    # clear upload inputs
    uploads_table(NULL)
    pairs(NULL)
    reps(NULL)
  })

  observe({
    invalidateLater(5000, session)
    handle_bulk_progress(quants, pquants, new_dataset)
  })



  # Selected Existing Dataset
  # ---

  # directory to existing dataset for exploration
  dataset_dir <- reactive({
    dataset_name <- input$dataset_name
    if (is.null(dataset_name)) return(NULL)

    datasets <- datasets()
    req(dataset_name, datasets)
    req(dataset_name %in% datasets$dataset_name)
    dir <- datasets[datasets$dataset_name == dataset_name, 'dataset_dir']
    file.path(project_dir(), dir)
  })

  numsv_path <- reactive(file.path(dataset_dir(), 'numsv.qs'))
  svobj_path <- reactive(file.path(dataset_dir(), 'svobj.qs'))

  # initialize svobj
  svobj_r <- reactiveVal()

  observe({
    svobj <- NULL
    svobj_path <- svobj_path()

    if (file.exists(svobj_path)) {
      svobj <- qs::qread(svobj_path)
    }

    svobj_r(svobj)
  })

  # initialize numsv
  numsv_r <- reactiveVal()
  maxsv_r <- reactiveVal()

  # initialize selected and max number of svs
  observe({
    new_dataset()
    numsv <- maxsv <- 0
    svobj <- svobj_r()
    if (!is.null(svobj$n.sv)) maxsv <- svobj$n.sv

    numsv_path <- numsv_path()
    if (file.exists(numsv_path)) numsv <- qs::qread(numsv_path)
    else qs::qsave(0, numsv_path)

    maxsv_r(maxsv)
    numsv_r(numsv)
  })

  # update number of surrogate variables slider
  observe({
    updateSliderInput(session, 'selected_nsv', value = numsv_r(), min = 0, max = maxsv_r())
  })

  observeEvent(input$selected_nsv, {
    numsv_r(input$selected_nsv)
    qs::qsave(input$selected_nsv, numsv_path())
  }, ignoreInit = TRUE)


  # update button icon with selected number of surrogate variables
  observe({
    updateActionButton(session, 'show_nsv', label = htmltools::doRenderTags(tags$span(input$selected_nsv, class='fa fa-fw')))
  })

  # toggle cell-type deconvolution
  show_deconv <- reactive(input$show_deconv %% 2 != 0)

  observe({
    shinyjs::toggleClass(id = "show_deconv", 'btn-primary', condition = show_deconv())
  })

  deconvForm <- callModule(
    deconvForm, 'deconv',
    show_deconv = show_deconv,
    new_dataset = new_dataset,
    sc_dir = sc_dir,
    bulk_dir = bulk_dir,
    tx2gene_dir = tx2gene_dir,
    dataset_dir = dataset_dir)

  return(list(
    dataset_name = reactive(input$dataset_name),
    dataset_dir = dataset_dir,
    est_prop = deconvForm$est_prop,
    show_deconv = show_deconv,
    selected_nsv = reactive(input$selected_nsv),
    numsv_r = numsv_r,
    svobj_r = svobj_r
  ))

}



#' Logic for differential expression analysis part of bulkForm
#'
#' @keywords internal
#' @noRd
bulkFormAnal <- function(input, output, session, project_dir, dataset_name, dataset_dir, gs_dir, explore_eset, numsv_r, svobj_r) {


  # run surrogate variable analysis if required
  observeEvent(explore_eset(), {

    eset <- explore_eset()
    pdata <- Biobase::pData(eset)
    if (is.null(eset)) return(NULL)

    prev_path <- file.path(dataset_dir(), 'pdata_explore_prev.qs')
    pdata_path <- file.path(dataset_dir(), 'pdata_explore.qs')

    if (file.exists(prev_path)) {
      prev <- qs::qread(prev_path)
      saved <- qs::qread(pdata_path)
      changed <- check_bulk_changed(prev, saved)
      req(changed)
    }

    group <- pdata$group
    req(length(unique(group)) > 1)

    # Create a Progress object
    progress <- Progress$new(session, min=0, max = 2)
    on.exit(progress$close())
    progress$set(message = "Running SVA", value = 1)

    mods <- crossmeta::get_sva_mods(eset)
    svobj <- crossmeta::run_sva(mods, eset)

    # add row names so that can check sync during dataset switch
    if (!is.null(svobj$sv)) row.names(svobj$sv) <- colnames(eset)

    # save current pdata_explore  so that can tell if changed
    file.copy(pdata_path, prev_path, overwrite = TRUE)


    # update saved svobj
    qs::qsave(svobj, file.path(dataset_dir(), 'svobj.qs'))

    svobj_r(svobj)
    numsv_r(NULL)

    progress$set(value = 2)

  })

  # Gene choices
  # ---

  have_groups <- reactive(!is.null(explore_eset()))
  observe(toggle('explore_genes_container', condition = have_groups()))


  # order by top table if have
  gene_choices <- reactive({
    tt <- bulkAnal$top_table()

    if (!is.null(tt)) {
      features <- row.names(tt)

      choices <- data.table::data.table(
        Feature = get_html_features(features),
        logFC = round(tt$logFC, 2),
        FDR = round(tt$adj.P.Val, 4),
        fdr =  signif(tt$adj.P.Val, digits = 3),
        feature = features
      )

    } else {
      eset <- explore_eset()
      features <- row.names(eset)

      choices <- data.table::data.table(
        Feature = get_html_features(features),
        feature = features
      )
    }

    return(choices)
  })


  output$gene_table <- DT::renderDataTable({

    gene_table <- gene_choices()
    if (is.null(gene_table)) return(NULL)

    # non-html feature column is hidden and used for search
    # different ncol if contrast
    cols <- colnames(gene_table)
    vis_targ <- which(cols %in% c('fdr', 'feature'))-1
    search_targs <- 0

    # prevent sort/filter when qc_first
    sort_targs <- 0
    filter <- list(position='top', clear = TRUE, vertical = TRUE, opacity = 0.85)

    fdr_col_idx <- which(cols == 'FDR')-1
    fdr_val_idx <- which(cols == 'fdr')-1

    fdr_title_js <- DT::JS(
      sprintf(paste0(
        "function (row, data, rowIndex) {",
        "  const cells = $('td', row);",
        "  $(cells[%s]).attr('title', data[%s]);",
        "}"), fdr_col_idx, fdr_val_idx
      )
    )

    qc_first <- all(colnames(gene_table) %in% c('Feature', 'feature'))
    if (qc_first) {
      sort_targs <- '_all'
      filter = list(position='none')
    }

    regex <- input$gene_search
    if (grepl(', ', regex)) regex <- format_comma_regex(regex)

    search_cols <- isolate(input$gene_table_search_columns)
    search_cols <- lapply(search_cols, function(str) if (str == '') return(NULL) else list(search = str))

    dt <- DT::datatable(
      gene_table,
      class = 'cell-border',
      rownames = FALSE,
      escape = FALSE, # to allow HTML in table
      extensions = c('Scroller'),
      filter = filter,
      options = list(
        keys = TRUE,
        select = list(style = "multiple", items = "row"),
        deferRender = TRUE,
        scroller = TRUE,
        scrollCollapse = TRUE,
        dom = '<"hidden"f>t',
        bInfo = 0,
        scrollY=250,
        search = list(regex = TRUE, search = regex),
        searchCols = search_cols,
        language = list(search = 'Select feature to plot:'),
        columnDefs = list(
          list(visible = FALSE, targets = vis_targ),
          list(searchable = FALSE, targets = search_targs),
          list(sortable = FALSE, targets = sort_targs)
        ),
        rowCallback = fdr_title_js
      )
    )

    return(dt)

  }, server = TRUE)

  # Comparison Groups and Differential Analyses
  # ---

  pdata <- reactive({
    req(explore_eset())
    Biobase::pData(explore_eset())
  })

  explore_genes <- reactive({
    rows <- input$gene_table_rows_selected

    choices <- gene_choices()
    return(choices$feature[rows])
  })


  bulkAnal <- callModule(bulkAnal, 'ds',
                         pdata = pdata,
                         dataset_name = dataset_name,
                         eset = explore_eset,
                         svobj = svobj_r,
                         numsv = numsv_r,
                         dataset_dir = dataset_dir,
                         gs_dir = gs_dir)


  return(list(
    explore_genes = explore_genes,
    contrast_groups = bulkAnal$contrast_groups
  ))
}

#' Logic for deconvolution form
#'
#' @keywords internal
#' @noRd
deconvForm <- function(input, output, session, show_deconv, new_dataset, sc_dir, bulk_dir, tx2gene_dir, dataset_dir) {
  include_options <- list(render = I('{option: contrastOptions, item: contrastItem}'))
  input_ids <- c('exclude_clusters', 'deconv_dataset', 'submit_deconv', 'deconv_method')

  est_prop <- reactiveVal()
  observeEvent(dataset_dir(), {new_deconv(NULL); est_prop(NULL)})
  observeEvent(input$deconv_method, {new_deconv(NULL); est_prop(NULL)})

  # show deconvolution form toggle
  observe({
    toggle(id = "deconv_form", anim = TRUE, condition = show_deconv())
  })

  prev_dataset <- reactive(qread.safe(file.path(sc_dir(), 'prev_dataset.qs')))

  # available single cell datasets for deconvolution
  ref_datasets <- reactive({
    datasets <- get_sc_dataset_choices(sc_dir(), prev_dataset())
    return(datasets)
  })


  # update reference dataset choices
  options <- list(
    render = I('{option: scDatasetOptions, item: scDatasetItem, optgroup_header: scDatasetOptGroup}'),
    searchField = c('optgroup', 'label'))

  observe({
    datasets <- ref_datasets()
    sel <- isolate(input$deconv_dataset)

    datasets <- add_optgroup_type(datasets)
    datasets <- datasets_to_list(datasets)
    updateSelectizeInput(session, 'deconv_dataset', selected = sel, choices = datasets, options = options)
  })


  deconv_dataset <- reactive({
    sel_idx <- input$deconv_dataset
    ds <- ref_datasets()
    ds$name[ds$value == sel_idx]
  })

  # get path to resolution subdir being used
  deconv_subdir <- reactive({
    anal_name <- deconv_dataset()
    if (is.null(anal_name)) return(NULL)
    req(anal_name)

    dataset_dir <- file.path(sc_dir(), anal_name)
    resoln_path <- file.path(dataset_dir, 'resoln.qs')
    resoln <- qs::qread(resoln_path)
    file.path(dataset_dir, get_resoln_dir(resoln))
  })

  annot <- reactive(qs::qread(file.path(deconv_subdir(), 'annot.qs')))

  # update exclude cluster choices
  include_choices <- reactive({
    clusters <- annot()
    get_cluster_choices(clusters, with_all = TRUE, resoln_dir = deconv_subdir())
  })

  observe({
    choices <- include_choices()
    updateSelectizeInput(session, 'exclude_clusters', choices = choices, options = include_options, server = TRUE)
  })

  new_deconv <- reactiveVal()
  is_disabled <- reactiveVal(FALSE)
  deconvs <- reactiveValues()
  pdeconvs <- reactiveValues()

  deconv_path <- reactive({
    deconv_hash <- calc_deconv_hash(deconv_dataset(), sc_dir(), input$exclude_clusters)
    deconv_fname <- paste0('deconv_', input$deconv_method, '_', deconv_hash, '.qs')
    file.path(dataset_dir(), deconv_fname)
  })

  is_integrated <- reactive({
    dataset_name <- deconv_dataset()
    req(dataset_name)
    integrated <- get_integrated_datasets(sc_dir())
    return(dataset_name %in% integrated)
  })

  observe({
    is.int <- is_integrated()
    choices <- c('MuSiC', 'DWLS')[c(is.int, TRUE)]
    updateSelectizeInput(session, 'deconv_method', choices = choices)
  })

  observeEvent(input$submit_deconv, {

    bulk_dataset_dir <- dataset_dir()
    scseq_dataset_name <- deconv_dataset()
    scseq_dataset_dir <- file.path(sc_dir(), scseq_dataset_name)
    exclude_clusters <- input$exclude_clusters
    method <- input$deconv_method

    if (!isTruthyAll(bulk_dataset_dir, scseq_dataset_name, method)) return(NULL)

    # check if already have
    deconv_path <- deconv_path()

    if (file.exists(deconv_path)) {
      new_deconv(deconv_path)

    } else {

      # disable inputs
      disableAll(input_ids)
      is_disabled(TRUE)

      # get method function
      method_fun <- c('MuSiC' = deconv_bulk_music, 'DWLS' = deconv_bulk_dwls)[[method]]

      deconvs[[deconv_path]] <- callr::r_bg(
        method_fun,
        package = 'dseqr',
        args = list(
          bulk_dataset_dir = bulk_dataset_dir,
          scseq_dataset_dir = scseq_dataset_dir,
          exclude_clusters = exclude_clusters,
          tx2gene_dir = tx2gene_dir,
          deconv_path = deconv_path
        ))

      progress <- Progress$new(max=3)
      progress$set(message = "Deconvoluting", value = 1)
      pdeconvs[[deconv_path]] <- progress

    }
  })

  observeEvent(new_deconv(), {

    deconv_path <- deconv_path()
    req(deconv_path)

    est <- qs::qread(deconv_path)
    est_cols <- as.numeric(colnames(est))

    annot <- annot()
    colnames(est) <- annot[est_cols]

    est_prop(est)
  })


  observe({
    invalidateLater(5000, session)
    handle_sc_progress(deconvs, pdeconvs, new_deconv)
  })

  # enable when complete (only one transfer at a time)
  observe({
    invalidateLater(1000, session)
    doing <- reactiveValuesToList(deconvs)
    doing <- names(doing)[!sapply(doing, is.null)]

    if (!length(doing) && is_disabled()) {
      enableAll(input_ids)
      is_disabled(FALSE)
    }
  })


  return(list(
    est_prop = est_prop
  ))
}

format_markers_for_dwls <- function(markers) {
  for (i in seq_along(markers)) {
    df <- markers[[i]]
    row.names(df) <- df$feature
    df <- dplyr::rename(df, avg_log2FC = .data$logFC, p_val_adj = .data$pval)
    markers[[i]] <- df
  }

  return(markers)
}

deconv_bulk_dwls <- function(bulk_dataset_dir, scseq_dataset_dir, exclude_clusters, tx2gene_dir, deconv_path) {

  # load markers
  resoln_dir <- load_resoln(scseq_dataset_dir)
  markers_path <- file.path(scseq_dataset_dir, resoln_dir, 'markers.qs')
  have.markers <- file.exists(markers_path)

  # need logs to get markers, counts for deconvolution
  scseq <- load_scseq_qs(scseq_dataset_dir, with_logs = !have.markers, with_counts = TRUE)

  species <- scseq@metadata$species
  is.human <- species == 'Homo sapiens'
  if (!is.human) scseq <- convert_species(scseq, tx2gene_dir, species)

  # subset to selected clusters
  if (length(exclude_clusters)) {
    all_clusters <- levels(scseq$cluster)
    keep_clusters <- setdiff(all_clusters, exclude_clusters)
    scseq <- scseq[, scseq$cluster %in% keep_clusters]
    scseq$cluster <- droplevels(scseq$cluster)
  }

  if (have.markers & is.human) {
    markers <- qs::qread(markers_path)
  } else {
    markers <- get_presto_markers(scseq)
  }

  markers <- markers[!names(markers) %in% exclude_clusters]
  markers <- format_markers_for_dwls(markers)

  # create signature for DWLS
  counts <- SingleCellExperiment::counts(scseq)
  sig <- DWLS::buildSignatureMatrix(counts, scseq$cluster, markers)

  eset <- qs::qread(file.path(bulk_dataset_dir, 'eset.qs'))
  bulk <- Biobase::exprs(eset)

  # trim to common genes
  tr <- DWLS::trimData(sig, bulk)

  # deconvolution
  res <- apply(tr$bulk, 2, function(x) DWLS::solveDampenedWLS(tr$sig, x))
  qs::qsave(t(res), deconv_path)
  return(TRUE)
}

calc_deconv_hash <- function(dataset_name, sc_dir, exclude_clusters, method = 'MuSiC') {
  dataset_dir <- file.path(sc_dir, dataset_name)
  resoln_name <- load_resoln(dataset_dir)
  resoln_dir <- file.path(dataset_dir, resoln_name)
  clusters_path <- file.path(resoln_dir, 'clusters.qs')
  clusters <- qs::qread(clusters_path)

  deconv_hash <- digest::digest(list(clusters, sort(exclude_clusters), dataset_name, method))
  return(deconv_hash)

}

deconv_bulk_music <- function(bulk_dataset_dir, scseq_dataset_dir, exclude_clusters, tx2gene_dir, deconv_path) {
  # MuSiC doesn't seem to properly import SingleCellExperiment::counts
  require('SingleCellExperiment')

  # load scseq (reference data)
  scseq <- load_scseq_qs(scseq_dataset_dir, with_logs = FALSE, with_counts = TRUE)

  species <- scseq@metadata$species
  if (species != 'Homo sapiens')
    scseq <- convert_species(scseq, tx2gene_dir, species)


  # subset to selected clusters
  if (length(exclude_clusters)) {
    all_clusters <- levels(scseq$cluster)
    keep_clusters <- setdiff(all_clusters, exclude_clusters)
    scseq <- scseq[, scseq$cluster %in% keep_clusters]
    scseq$cluster <- droplevels(scseq$cluster)
  }

  # load bulk ExpressionSet (to deconvolute)
  eset <- qs::qread(file.path(bulk_dataset_dir, 'eset.qs'))

  # MuSiC deconvolution
  bulk.mtx <- Biobase::exprs(eset)
  est.prop <- MuSiC::music_prop(bulk.mtx, sc.sce = scseq, clusters = 'cluster', samples = 'batch')
  est.prop <- est.prop$Est.prop.weighted

  qs::qsave(est.prop, deconv_path)

  return(TRUE)
}

deconv_bulk_dtangle <- function(bulk_dataset_dir, scseq_dataset_dir, exclude_clusters, tx2gene_dir, deconv_path) {

  # load scseq (reference data)
  scseq <- load_scseq_qs(scseq_dataset_dir, with_logs = TRUE)
  scseq <- downsample_clusters(scseq)

  species <- scseq@metadata$species
  if (species != 'Homo sapiens')
    scseq <- convert_species(scseq, tx2gene_dir, species)


  # subset to selected clusters
  if (length(exclude_clusters)) {
    all_clusters <- levels(scseq$cluster)
    keep_clusters <- setdiff(all_clusters, exclude_clusters)
    scseq <- scseq[, scseq$cluster %in% keep_clusters]
    scseq$cluster <- droplevels(scseq$cluster)
  }

  # load bulk ExpressionSet (to deconvolute)
  eset <- qs::qread(file.path(bulk_dataset_dir, 'eset.qs'))
  pdata <- qs::qread(file.path(bulk_dataset_dir, 'pdata_explore.qs'))
  vsd_path <- file.path(bulk_dataset_dir, 'vsd.qs')
  eset <- normalize_eset(eset, pdata, vsd_path)

  svobj <- qs::qread(file.path(bulk_dataset_dir, 'svobj.qs'))
  numsv <- qs::qread(file.path(bulk_dataset_dir, 'numsv.qs'))
  adj_path <- file.path(bulk_dataset_dir, paste0('adjusted_', numsv, 'svs.qs'))
  keep_path <- file.path(bulk_dataset_dir, paste0('iqr_keep_', numsv, 'svs.qs'))

  # adjust for pairs/surrogate variables
  eset <- add_adjusted(eset, adj_path, svobj, numsv)

  # use SYMBOL as annotation
  # keep unique symbol based on row IQRs
  eset <- iqr_replicates(eset, keep_path)

  # get normalized values
  adj <- Biobase::assayDataElement(eset, 'adjusted')

  # common genes only
  commongenes <- intersect(rownames(adj), rownames(scseq))
  adj <- adj[commongenes, ]
  scseq <- scseq[commongenes, ]

  y <- cbind(as.matrix(SingleCellExperiment::logcounts(scseq)), adj)
  y <- limma::normalizeBetweenArrays(y)
  y <- t(y)

  include_clusters <- levels(scseq$cluster)
  pure_samples <- list()
  for (i in seq_along(include_clusters))
    pure_samples[[include_clusters[i]]] <- which(scseq$cluster == include_clusters[i])

  # markers for each included cluster
  marker_list = dtangle::find_markers(y,
                                      pure_samples = pure_samples,
                                      data_type = "rna-seq",
                                      marker_method='ratio')

  # use markers in top 10th quantile with a minimum of 3
  q = 0.1
  quantiles = lapply(marker_list$V,function(x) stats::quantile(x,1-q))
  K = length(pure_samples)
  n_markers = sapply(seq_len(K),function(i){
    min(
      length(marker_list$L[[i]]),
      max(3, which(marker_list$V[[i]] > quantiles[[i]]))
    )
  })

  # run deconvolution and get proportion estimates
  marks <- marker_list$L
  dc <- dtangle::dtangle(y,
                         pure_samples = pure_samples,
                         n_markers = n_markers,
                         data_type = 'rna-seq',
                         markers = marks)

  dc <- dc$estimates[colnames(eset), ]
  qs::qsave(dc, deconv_path)

  return(TRUE)
}



#' Logic for differential expression analysis table
#'
#' @keywords internal
#' @noRd
bulkExploreTable <- function(input, output, session, eset, up_annot, project_dir, dataset_dir, dataset_name, svobj_r, numsv_r) {

  # things user will update and return
  pdata_r <- reactiveVal()
  pdata_path <- reactive(file.path(dataset_dir(), 'pdata_explore.qs'))

  eset_pdata <- reactive({
    eset <- eset()
    req(eset)

    pdata <- Biobase::pData(eset) %>%
      tibble::add_column(Group = NA, 'Group name' = NA, Title = colnames(eset),  Pair = NA, .before = 1)

    return(pdata)
  })

  # update when upload new annotation
  observeEvent(up_annot(), {
    up <- up_annot()
    pdata_path <- pdata_path()
    req(up)


    if (file.exists(pdata_path)) {
      prev <- qs::qread(pdata_path)
      changed <- check_bulk_changed(prev, up)

      if (changed) {
        remove_dataset_files(dataset_dir())
        svobj_r(NULL)
        numsv_r(NULL)
      }
    }

    pdata_r(up)
    qs::qsave(up, pdata_path)
  })


  # reset when new eset
  observe({
    pdata_path <- pdata_path()
    pdata <- eset_pdata()
    req(pdata, pdata_path)

    # load pdata from previous if available
    if (file.exists(pdata_path)) {
      pdata <- qs::qread(pdata_path)

      # TODO remove
      # need for legacy purposes
      pair <- pdata$Pair
      if (is.null(pair)) pair <- pdata$pair
      if (is.null(pair)) pair <- NA
      pdata$Pair <- NULL
      pdata <- tibble::add_column(pdata, Pair = pair, .after = 'Title')
    }

    pdata_r(pdata)
  })



  # format pdata for table
  html_pdata <- reactive({

    # things that trigger update
    pdata <- pdata_r()
    group <- pdata$Group
    name  <- pdata$`Group name`
    req(pdata)


    # update pdata Group column
    not.na <- !is.na(group)
    group_nums  <- unique(group[not.na])
    group_names <- unique(name[not.na])
    ind <- order(group_nums)

    group_nums  <- group_nums[ind]
    group_names <- group_names[ind]

    group_colors <- get_palette(group_nums)
    for (i in seq_along(group_nums)) {
      group_num <- group_nums[i]
      group_name <- group_names[i]

      # plotly color bug when two groups
      color <- group_colors[group_num]

      rows <- which(group == group_num)
      pdata[rows, 'Group'] <- paste('<div style="background-color:', color, ';"></div>')
      pdata[rows, 'Group name'] <- group_name
    }

    return(pdata)
  })

  # redraw table when new pdata
  output$pdata <- DT::renderDataTable({

    pdata <- html_pdata()
    hide_target <- which(colnames(pdata) %in% c('lib.size', 'norm.factors', 'pair')) - 1

    DT::datatable(
      pdata,
      class = 'cell-border dt-fake-height',
      rownames = FALSE,
      escape = FALSE, # to allow HTML in table
      options = list(
        columnDefs = list(
          list(className = 'dt-nopad', targets = 0),
          list(targets = hide_target, visible=FALSE)),
        scrollX = TRUE,
        paging = FALSE,
        bInfo = 0
      )
    )
  })


  return(list(
    pdata = pdata_r
  ))

}

#' Logic for downloading and uploading bulk annotation
#'
#' @keywords internal
#' @noRd
bulkAnnot <- function(input, output, session, dataset_name, annot) {


  observeEvent(input$click_up, {
    error_msg(NULL)
  })

  observeEvent(input$click_dl, {
    error_msg(NULL)
  })

  fname <- reactive(paste0(dataset_name(), '_annot.csv'))

  output$dl_annot <- downloadHandler(
    filename = fname,
    content = function(con) {

      utils::write.csv(format_dl_annot(annot()), con, row.names = FALSE)
    }
  )

  # uploaded annotation
  up_annot <- reactiveVal()
  error_msg <- reactiveVal()

  # reset when change dataset
  observe({
    dataset_name()
    up_annot(NULL)
  })

  observe({
    msg <- error_msg()
    html('error_msg', html = msg)
    shinyjs::toggleClass('validate-up', 'has-error', condition = isTruthy(msg))
  })

  observeEvent(input$up_annot, {
    ref <- annot()
    req(ref)

    infile <- input$up_annot
    if (!isTruthy(infile)){
      res <- msg <- NULL

    } else {
      res <- utils::read.csv(infile$datapath, check.names = FALSE, stringsAsFactors = FALSE)
      msg <- validate_up_annot(res, annot())

      res <- if (is.null(msg)) format_up_annot(res, ref) else NULL
    }

    error_msg(msg)
    up_annot(res)
  })

  return(up_annot)

}


#' Logic for bulk group analyses for Bulk, Drugs, and Pathways tabs
#'
#' @keywords internal
#' @noRd
bulkAnal <- function(input, output, session, pdata, dataset_name, eset, numsv, svobj, dataset_dir, gs_dir, is_bulk = function()TRUE) {
  contrast_options <- list(render = I('{option: bulkContrastOptions, item: bulkContrastItem}'))
  input_ids <- c('click_dl', 'contrast_groups')

  # hide group selector if no groups
  have_groups <- reactive(!is.null(eset()))
  observe(toggle('run_anal_container', condition = have_groups()))


  # group levels used for selecting test and control groups
  group_levels <- reactive({
    req(is_bulk())
    get_group_levels(pdata())
  })

  group_colors <- reactive(get_palette(group_levels()))

  group_choices <- reactive({

    data.frame(
      name = group_levels(),
      value = group_levels(),
      color = group_colors(), stringsAsFactors = FALSE
    )
  })

  # re-select previous group choices
  contrast_groups <- reactiveVal()
  prev_path <- reactive(file.path(dataset_dir(), 'prev_groups.qs'))


  observe(contrast_groups(qread.safe(prev_path())))

  observeEvent(input$contrast_groups, {
    groups <- input$contrast_groups
    if (is.null(groups)) return(NULL)

    attr(groups, 'dataset_name') <- dataset_name()
    prev <- contrast_groups()

    if (!identical(prev, groups)) {
      qs::qsave(groups, prev_path())
      contrast_groups(groups)
    }

  })

  observe({
    prev <- qread.safe(prev_path())
    updateSelectizeInput(session,
                         'contrast_groups',
                         choices = group_choices(),
                         selected = prev,
                         server = TRUE,
                         options = contrast_options)
  })

  valid_contrast <- reactive({
    groups <- contrast_groups()
    dataset <- attr(groups, 'dataset_name')
    length(groups) == 2 && !is.null(dataset) && dataset == dataset_name()
  })

  anal_name <- reactive({
    req(valid_contrast())
    groups <- contrast_groups()
    return(paste0(groups[1], '_vs_', groups[2]))
  })

  # path to lmfit and drug query results
  numsv_str <- reactive(paste0(numsv(), 'svs'))

  lmfit_path <- reactive({
    req(is_bulk())
    lmfit_file <- paste0('lm_fit_', numsv_str(), '.qs')
    file.path(dataset_dir(), lmfit_file)
  })

  drug_paths <- reactive({
    suffix <- paste(anal_name(), numsv_str(), sep = '_')
    get_drug_paths(dataset_dir(), suffix)
  })

  goana_path <- reactive({
    fname <- paste0('goana_', anal_name(), '_',
                    numsv_str(), '_', input$min_abs_logfc, '_', input$max_fdr, '.qs')
    file.path(dataset_dir(), fname)
  })

  # do we have drug query results?
  saved_drugs <- reactive(file.exists(drug_paths()$cmap))

  # load lm_fit if saved or run limma
  lm_fit <- reactive({

    if (file.exists(lmfit_path())) {
      lm_fit <- qs::qread(lmfit_path())

    } else {

      # visual that running
      disableAll(input_ids)

      progress <- Progress$new(session, min=0, max = 3)
      on.exit(progress$close())

      # check for previous lm_fit
      dataset_dir <- dataset_dir()
      numsv <- numsv()

      # get what need
      eset <- eset()
      svobj <- svobj()
      req(eset, dataset_dir)

      prev_anal <- list(pdata = Biobase::pData(eset))

      # run differential expression
      progress$set(message = "Fitting limma model", value = 1)
      eset <- crossmeta::run_limma_setup(eset, prev_anal)
      lm_fit <- crossmeta::run_limma(eset,
                                     svobj = svobj,
                                     numsv = numsv,
                                     filter = FALSE)

      save_lmfit(lm_fit, dataset_dir, numsv = numsv)

      # visual that done
      progress$inc(1)
      enableAll(input_ids)
    }

    return(lm_fit)
  })

  drug_queries <- reactive({

    if (!valid_contrast()) {
      res <- NULL

    } else if (saved_drugs()) {
      paths <- drug_paths()
      res <- lapply(paths, qs::qread)

    } else {
      tt <- top_table()
      if (!isTruthy(tt)) return(NULL)

      disableAll(input_ids)
      progress <- Progress$new(session, min = 0, max = 3)
      progress$set(message = "Querying drugs", value = 1)
      on.exit(progress$close())

      es <- load_drug_es()
      progress$inc(1)
      res <- run_drug_queries(tt, drug_paths(), es)
      progress$inc(1)
      enableAll(input_ids)
    }
    return(res)
  })

  # differential expression top table
  top_table <- reactive({
    if (!valid_contrast()) return(NULL)
    lm_fit <- lm_fit()
    groups <- contrast_groups()

    # loses sync when groups selected and change dataset
    if (!all(groups %in% colnames(lm_fit$mod))) return(NULL)

    crossmeta::get_top_table(lm_fit, groups)
  })

  # go/kegg pathway result
  path_res <- reactive({
    req(valid_contrast())
    goana_path <- goana_path()

    if (file.exists(goana_path)) {
      res <- qs::qread(goana_path)

    } else {
      max_fdr <- input$max_fdr
      min_abs_logfc <- input$min_abs_logfc

      lm_fit <- lm_fit()

      # visual that running
      disableAll(input_ids)
      progress <- Progress$new(session, min=0, max = 2)
      on.exit(progress$close())
      progress$set(message = "Running pathway analysis", value = 1)
      groups <- contrast_groups()

      # loses sync when groups selected and change dataset
      req(all(groups %in% colnames(lm_fit$mod)))

      groups <- make.names(groups)

      contrast <- paste0(groups[1], '-', groups[2])
      ebfit <- crossmeta::fit_ebayes(lm_fit, contrast)
      res <- get_path_res(ebfit, goana_path, gs_dir, max_fdr = max_fdr, min_abs_logfc = min_abs_logfc)
      qs::qsave(res, goana_path)

      progress$inc(1)
      enableAll(input_ids)
    }

    return(res)
  })



  # enable download
  observe({
    shinyjs::toggleState('download', condition = valid_contrast())
  })

  fname_str <- reactive({
    numsv_str <- paste0(numsv(), 'SV')
    paste('bulk', dataset_name(), anal_name(), numsv_str, sep='_')
  })

  filter_str <- reactive({
    fdr_str <- paste0('FDR', input$max_fdr)
    logfc_str <- paste0('logFC', input$min_abs_logfc)

    paste(fdr_str, logfc_str, sep='_')
  })


  data_fun <- function(file) {
    #go to a temp dir to avoid permission issues
    owd <- setwd(tempdir())
    on.exit(setwd(owd))

    tt_fname <- 'top_table.csv'
    tt_fname_all <- 'top_table_all.csv'
    goup_fname <- 'go_up.csv'
    godn_fname <- 'go_dn.csv'

    path_res <- path_res()
    utils::write.csv(filtered_tt(), tt_fname)
    utils::write.csv(top_table(), tt_fname_all)
    utils::write.csv(path_res$up, goup_fname)
    utils::write.csv(path_res$dn, godn_fname)

    #create the zip file
    utils::zip(file, c(tt_fname, tt_fname_all, goup_fname, godn_fname))
  }

  output$download <- downloadHandler(
    filename = function() {
      paste0(fname_str(), '_', filter_str(), '_', Sys.Date(), '.zip')
    },
    content = data_fun
  )

  prev_max_fdr <- reactiveVal(0.05)
  prev_min_abs_logfc <- reactiveVal(0)

  observeEvent(input$click_dl, {
    showModal(downloadResultsModal(session, prev_max_fdr(), prev_min_abs_logfc()))
  })

  observe({
    max_fdr <- input$max_fdr
    req(is.numeric(max_fdr))
    prev_max_fdr(input$max_fdr)
  })

  observe({
    min_abs_logfc <- input$min_abs_logfc
    req(is.numeric(min_abs_logfc))
    prev_min_abs_logfc(input$min_abs_logfc)
  })

  callModule(volcanoPlotOutput, 'volcano_plot',
             top_table = top_table,
             max_fdr = reactive(input$max_fdr),
             min_abs_logfc = reactive(input$min_abs_logfc))

  filtered_tt <- reactive({
    tt <- top_table()
    min_abs_logfc <- input$min_abs_logfc
    max_fdr <- input$max_fdr

    tt <- tt[tt$adj.P.Val < max_fdr, ]
    tt <- tt[abs(tt$logFC) > min_abs_logfc, ]
    return(tt)
  })

  output$ngenes_up <- renderText({
    tt <- filtered_tt()
    tt <- tt[tt$logFC > 0,]
    return(format(nrow(tt), big.mark =','))
  })

  output$ngenes_dn <- renderText({
    tt <- filtered_tt()
    tt <- tt[tt$logFC < 0,]
    return(format(nrow(tt), big.mark =','))
  })

  # download can timeout so get objects before clicking
  observeEvent(input$confirm_dl_anal, {
    removeModal()
    pres <- path_res()
    shinyjs::click("download")
  })

  return(list(
    name = anal_name,
    contrast_groups = contrast_groups,
    lm_fit = lm_fit,
    top_table = top_table,
    drug_queries = drug_queries,
    path_res = path_res
  ))
}



volcanoPlotOutput <- function(input, output, session, top_table, max_fdr, min_abs_logfc) {


  output$plot <- shiny::renderPlot({
    tt <- top_table()
    min_abs_logfc <- min_abs_logfc()
    max_fdr <- max_fdr()

    tt$gene_name <- row.names(tt)

    tt$color <- rgb(0.76, 0.76, 0.76, .3)
    dn <- tt$logFC < -min_abs_logfc & tt$adj.P.Val < max_fdr
    up <- tt$logFC >  min_abs_logfc & tt$adj.P.Val < max_fdr
    tt$color[dn] <- rgb(0, 0, 1, .3)
    tt$color[up] <- rgb(1, 0, 0, .3)

    tt$gene_name[!dn & !up] <- NA

    max_abs_logfc <- max(abs(tt$logFC))
    xlims <- c(-max_abs_logfc-0.25, max_abs_logfc+0.25)

    mar <- par()$mar
    mar[3] <- 1
    par(mar=mar)

    plot(tt$logFC, -log10(tt$adj.P.Val), pch=19, col=tt$color, xlim=xlims,
         ylab=bquote(~-log[10] ~ FDR),  xlab= "logFC", bty="l")

    abline(h = -log10(max_fdr), lty=2)
    abline(v = -min_abs_logfc, lty=2)
    abline(v = min_abs_logfc, lty=2)

  }, width = 380, height = 300)



}

#' Logic to setup explore_eset for Bulk Data plots
#'
#' @keywords internal
#' @noRd
exploreEset <- function(eset, dataset_dir, explore_pdata, numsv, svobj) {

  vsd_path <- reactive(file.path(dataset_dir(), 'vsd.qs'))
  adj_path <- reactive(file.path(dataset_dir(), paste0('adjusted_', numsv(), 'svs.qs')))
  keep_path <- reactive(file.path(dataset_dir(), paste0('iqr_keep_', numsv(), 'svs.qs')))


  norm_eset <- reactive({
    eset <- eset()
    pdata <- explore_pdata()
    vsd_path <- vsd_path()
    normalize_eset(eset, pdata, vsd_path)
  })

  # explore_eset used for all plots
  explore_eset <- reactive({
    eset <- norm_eset()
    if (is.null(eset)) return(NULL)

    numsv <- numsv()

    # adjust for pairs/surrogate variables
    svobj <- svobj()

    # can lose sync when switching datasets
    if (!is.null(svobj$sv)) {
      req(numsv <= svobj$n.sv)
      req(identical(row.names(svobj$sv), colnames(eset)))
    }

    eset <- add_adjusted(eset, adj_path(), svobj, numsv)

    # use SYMBOL as annotation
    # keep unique symbol based on row IQRs
    eset <- iqr_replicates(eset, keep_path())

    return(eset)
  })
  return(explore_eset)
}

normalize_eset <- function(eset, pdata, vsd_path) {

  # pdata and eset lose sync when switch datasets
  in.sync <- identical(colnames(eset), row.names(pdata))
  if (!in.sync) return(NULL)

  # need that pdata_explore has more than two groups
  keep <- row.names(pdata)[!is.na(pdata$Group)]
  pdata <- pdata[keep, ]
  if (!length(unique(pdata$Group)) > 1) return(NULL)
  if (!length(keep) > 2) return(NULL)

  # determine if this is rna seq data
  rna_seq <- 'norm.factors' %in% colnames(Biobase::pData(eset))

  # subset eset and add group/pair
  eset <- eset[, keep]
  pdata$group <- pdata$`Group name`

  pair <- pdata$Pair
  if (any(!is.na(pair))) pdata$pair <- pair

  Biobase::pData(eset) <- pdata

  # filter rows by expression
  if (rna_seq) eset <- rkal::filter_genes(eset)

  # cpm normalize
  eset <- add_vsd(eset, vsd_path, rna_seq)
  return(eset)
}

