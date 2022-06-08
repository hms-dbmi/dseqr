# to debug:
# 1) add browser() inside test_that()
# 2) run devtools::test_active_file()

suppressWarnings(suppressPackageStartupMessages({
  library(shinytest2)
  library(SingleCellExperiment)
  library(dseqr)
}))

# create single-cell dataset for upload
mock_10x_files <- function(dir_name) {

  unlink(dir_name, recursive = TRUE)
  set.seed(0)
  ngenes <- 3000
  sce <- scDblFinder::mockDoubletSCE(ncells = c(50, 75, 25, 50, 75), ngenes = ngenes)
  counts <- as(counts(sce), 'dgCMatrix')
  t2g <- dseqr.data::load_tx2gene()
  DropletUtils::write10xCounts(
    path = dir_name,
    x = counts,
    gene.id = t2g$gene_id[seq_len(ngenes)],
    gene.symbol = t2g$gene_name[seq_len(ngenes)],
  )

  for (file in list.files(dir_name, full.names = TRUE)) R.utils::gzip(file)
}

# 20 mins
timeout <- 1000*60*20

test_that("{shinytest2} recording: Single-Cell Tab", {
  sample_dir <- 'mock_10x'
  data_dir <- 'test_data_dir'
  unlink(sample_dir, recursive = TRUE)
  unlink(data_dir, recursive = TRUE)

  suppressWarnings(
    app <- AppDriver$new(
      name = "import_sc",
      options = list(shiny.maxRequestSize = 30*1024*1024^2),
      width = 1619,
      height = 909,
      seed = 0,
      load_timeout = timeout)
  )

  list_files <- function()
    file.path(data_dir, list.files(data_dir, recursive = TRUE, all.files = TRUE, include.dirs = TRUE))

  # created expected files/folders
  init_files <- list_files()
  expect_snapshot(init_files)

  mock_10x_files(sample_dir)
  on.exit(unlink(sample_dir, recursive = TRUE), add = TRUE)
  on.exit(unlink(data_dir, recursive = TRUE), add = TRUE)
  on.exit(app$stop(), add = TRUE)


  app$upload_file(`sc-form-dataset-up_raw` = file.path(sample_dir, "barcodes.tsv.gz"))
  app$upload_file(`sc-form-dataset-up_raw` = file.path(sample_dir, "genes.tsv.gz"))
  app$upload_file(`sc-form-dataset-up_raw` = file.path(sample_dir, "matrix.mtx.gz"))
  app$wait_for_idle()
  app$set_inputs(`sc-form-dataset-up_table_rows_selected` = 1:3, allow_no_input_binding_ = TRUE)
  app$set_inputs(`sc-form-dataset-sample_name` = "mock_10x")

  # clears sample name after add
  app$click("sc-form-dataset-add_sample")
  expect_equal(app$get_value(input = 'sc-form-dataset-sample_name'), '')

  # auto detects as human
  app$click("sc-form-dataset-import_datasets")
  app$wait_for_idle()
  expect_equal(app$get_value(input = 'sc-form-dataset-import_species'), 'Homo sapiens')

  # auto selected 'all and none'
  metrics <- app$wait_for_value(input = 'sc-form-dataset-qc_metrics')
  expect_equal(metrics, 'all and none')

  # switch to 'all'
  app$set_inputs(`sc-form-dataset-qc_metrics` = 'all')
  app$wait_for_idle()

  # creates mock_10x sample
  app$click("sc-form-dataset-confirm_import_datasets")
  app$wait_for_value(export = 'sc-form-dataset-dataset_names', timeout = timeout, ignore = list(character(0)))
  expect_setequal(app$get_value(export = 'sc-form-dataset-dataset_names'), 'mock_10x')

  # didn't auto-select dataset
  expect_equal(app$get_value(input = 'sc-form-dataset-selected_dataset'), '')

  # added files/folder for dataset
  dataset_files <- setdiff(list_files(), init_files)
  expect_snapshot(dataset_files)

  # clusters show up after selecting dataset
  app$set_inputs(`sc-form-dataset-selected_dataset` = "1", wait_ = FALSE)
  app$wait_for_value(export = 'sc-form-sample_clusters-annot', ignore = list(NULL))
  expect_equal(
    app$get_value(export = 'sc-form-sample_clusters-annot'),
    as.character(1:5))

  # saved prev_dataset
  prev_file <- setdiff(list_files(), c(init_files, dataset_files))
  expect_equal(qs::qread(prev_file), sample_dir)

  # didn't auto-select cluster
  app$wait_for_idle()
  expect_equal(app$get_value(input = 'sc-form-cluster-selected_cluster'), '')
  app$set_inputs(`sc-form-cluster-selected_cluster` = "1")

  # created markers file
  app$wait_for_value(export = 'sc-form-cluster-have_selected_markers', timeout = timeout, ignore = list(FALSE, NULL))
  all_prev_files <- c(init_files, dataset_files, prev_file)
  markers_file <- setdiff(list_files(), all_prev_files)
  expect_equal(markers_file, "test_data_dir/test_user/default/single-cell/mock_10x/snn1/markers.qs")
  expect_equal(unname(tools::md5sum(markers_file)), "fadb32f7c1c082fdefc057158e30e5eb")

  # clusters and markers plot have same coordinates
  app$set_inputs(`sc-form-gene_clusters-gene_table_rows_selected` = 1, allow_no_input_binding_ = TRUE)
  cluster_coords <- jsonlite::fromJSON(app$wait_for_value(output = 'sc-cluster_plot-cluster_plot'))$x$coords

  marker_plot <- app$wait_for_value(output = 'sc-marker_plot_cluster-marker_plot', timeout = timeout)
  markers_coords <- jsonlite::fromJSON(marker_plot)$x$coords

  if (is.null(markers_coords)) {
    marker_plot <- app$wait_for_value(output = 'sc-marker_plot_cluster-marker_plot', timeout = timeout, ignore = list(marker_plot))
    markers_coords <- jsonlite::fromJSON(marker_plot)$x$coords
  }

  cluster_x_ord <- order(cluster_coords$x)
  markers_x_ord <- order(markers_coords$x)

  # same x coords exist, y is same for each x
  expect_equal(cluster_coords$x[cluster_x_ord], markers_coords$x[markers_x_ord])
  expect_equal(cluster_coords$y[cluster_x_ord], markers_coords$y[markers_x_ord])

  # can change cluster name
  app$click("sc-form-cluster-show_rename")
  app$set_inputs(`sc-form-cluster-new_cluster_name` = "CD14 Mono")
  app$click("sc-form-cluster-rename_cluster")
  app$wait_for_idle()
  expect_equal(
    app$get_value(export = 'sc-form-sample_clusters-annot'),
    c("CD14 Mono", as.character(2:5)))

  # can compare one cluster vs one other cluster
  app$click("sc-form-cluster-show_contrasts")
  contrast_choices <- app$get_value(export = "sc-form-cluster-choices")
  expect_snapshot(contrast_choices)

  app$set_inputs(`sc-form-cluster-selected_cluster` = '1-vs-2')
  expect_equal(
    app$wait_for_value(input = 'sc-form-cluster-selected_cluster', timeout = timeout, ignore = list('1')),
    '1-vs-2'
  )
  app$set_inputs(`sc-form-cluster-selected_cluster` = '1')
  app$click("sc-form-cluster-show_contrasts")

  marker_colors <- app$wait_for_value(export = 'sc-marker_plot_cluster-colors', timeout = timeout)
  expect_length(unique(marker_colors), 94)

  # can set custom boolean metric
  app$click("sc-form-gene_clusters-show_custom_metric")
  app$wait_for_idle()
  app$set_inputs(`sc-form-gene_clusters-custom_metric` = "EFNB1>0")

  marker_colors <- app$wait_for_value(export = 'sc-marker_plot_cluster-colors', timeout = timeout, ignore = list(marker_colors))
  expect_length(unique(marker_colors), 2)

  # can view expression of gene using custom metric
  app$set_inputs(`sc-form-gene_clusters-custom_metric` = "EFNB1")
  marker_colors <- app$wait_for_value(export = 'sc-marker_plot_cluster-colors', timeout = timeout, ignore = list(marker_colors))
  ncols_efnb1 <- length(unique(marker_colors))
  expect_gt(ncols_efnb1, 2)

  # can use cluster for metric
  app$set_inputs(`sc-form-gene_clusters-custom_metric` = "cluster==1")
  marker_colors <- app$wait_for_value(export = 'sc-marker_plot_cluster-colors', timeout = timeout, ignore = list(marker_colors))
  expect_length(unique(marker_colors), 2)

  # closing custom metric plots previously selected gene
  app$click("sc-form-gene_clusters-show_custom_metric")
  marker_colors <- app$wait_for_value(export = 'sc-marker_plot_cluster-colors', timeout = timeout, ignore = list(marker_colors))
  expect_length(unique(marker_colors), ncols_efnb1)

  # re-opening custom metric plots previous specified metric
  app$click("sc-form-gene_clusters-show_custom_metric")
  marker_colors <- app$wait_for_value(export = 'sc-marker_plot_cluster-colors', timeout = timeout, ignore = list(marker_colors))
  expect_length(unique(marker_colors), 2)

  # can add custom metric
  prev_files <- list_files()
  app$click("sc-form-gene_clusters-save_custom_metric")
  expect_equal(app$get_value(input = "sc-form-gene_clusters-custom_metric"), "")
  saved_metric_files <- setdiff(list_files(), prev_files)
  expect_snapshot(saved_metric_files)

  # custom metric was saved
  saved_metrics <- qs::qread(saved_metric_files)
  expect_equal(colnames(saved_metrics), "cluster==1")

  # closing custom metric after save shows no marker plot
  app$click("sc-form-gene_clusters-show_custom_metric")
  marker_colors <- app$wait_for_value(export = 'sc-marker_plot_cluster-colors', timeout = timeout, ignore = list(marker_colors))
  expect_null(marker_colors)

  # can change resolution of a dataset
  app$click("sc-form-dataset-show_label_resoln")
  prev_files <- list_files()
  app$set_inputs("sc-form-resolution-resoln" = 2, wait_ = FALSE)
  app$wait_for_idle()
  change_resoln_files <- setdiff(list_files(), prev_files)
  expect_snapshot(change_resoln_files)
  app$click("sc-form-dataset-show_label_resoln")

  # can subset a dataset
  dataset_name <- app$get_value(export = "sc-form-dataset-dataset_name")
  expect_equal(dataset_name, "mock_10x")

  dataset_names <- app$get_value(export = "sc-form-dataset-dataset_names")
  app$click("sc-form-dataset-show_subset")
  app$set_inputs("sc-form-subset-subset_features" = c(1, 2, 3, 4))
  app$click("sc-form-subset-toggle_exclude")
  app$set_inputs("sc-form-subset-subset_name" = "1234")
  app$click("sc-form-subset-submit_subset")
  app$wait_for_idle()
  app$click("sc-form-subset-confirm_subset")
  app$click("sc-form-dataset-show_subset")

  # subset dataset has been added
  # dataset name duplicated as also in "previous"
  new_dataset_names <- app$wait_for_value(export = 'sc-form-dataset-dataset_names', timeout = timeout, ignore = list(dataset_names))
  expect_setequal(new_dataset_names, c('mock_10x', 'mock_10x_1234'))

  # check that previous dataset still selected
  expect_equal(dataset_name, "mock_10x")

  # check that new dataset has cells from selected clusters
  # load_scseq_qs so that attaches clusters from above resolution change
  mock10x <- load_scseq_qs('test_data_dir/test_user/default/single-cell/mock_10x')
  mock10x_ncells <- ncol(mock10x)
  expect_ncells <- sum(mock10x$cluster %in% 1:4)

  mock10x_1234 <- load_scseq_qs('test_data_dir/test_user/default/single-cell/mock_10x_1234')
  mock10x_1234_ncells <- ncol(mock10x_1234)
  expect_equal(mock10x_1234_ncells, expect_ncells)

  # run label transfer from initial dataset to subset dataset
  app$click("sc-form-dataset-show_label_resoln")
  app$set_inputs(`sc-form-dataset-selected_dataset` = 2, wait_ = FALSE)
  new_dataset_name <- app$get_value(export = "sc-form-dataset-dataset_name")
  expect_equal(new_dataset_name, 'mock_10x_1234')

  init_annot <- app$get_value(export = "sc-form-sample_clusters-annot")
  expect_equal(init_annot, as.character(1:4))

  app$set_inputs("sc-form-transfer-ref_name" = "mock_10x", wait_ = FALSE)
  pred_annot <- app$wait_for_value(export = "sc-form-transfer-pred_annot", timeout = timeout)
  expect_equal(pred_annot, c('CD14 Mono', 2:4))
  app$click("sc-form-transfer-overwrite_annot")
  app$wait_for_idle()
  app$click("sc-form-transfer-confirm_overwrite")

  new_annot <- app$wait_for_value(export = "sc-form-sample_clusters-annot", timeout = timeout, ignore = list(init_annot))
  expect_equal(new_annot, c('CD14 Mono', 2:4))

  # can integrate two datasets
  pre_integration_files <- list_files()
  app$click("datasets_dropdown")
  app$wait_for_idle()
  app$click("integrate_dataset")
  app$wait_for_idle()
  app$set_inputs("sc-form-integration-integration_datasets" = c("mock_10x", "mock_10x_1234"), wait_ = FALSE)
  app$set_inputs("sc-form-integration-integration_name" = c("mock_10x_integrated"), wait_ = FALSE)
  app$wait_for_idle()
  app$click("sc-form-integration-submit_integration")
  current_datasets <- app$get_value(export = 'sc-form-dataset-dataset_names')
  app$wait_for_value(export = 'sc-form-dataset-dataset_names', timeout = timeout, ignore = list(current_datasets))

  integrated_dataset_name <- setdiff(app$get_value(export = 'sc-form-dataset-dataset_names'), current_datasets)
  expect_snapshot(integrated_dataset_name)

  integrated_files <- setdiff(list_files(), pre_integration_files)
  expect_snapshot(integrated_files)

  # previous dataset still selected
  expect_equal(app$get_value(export = "sc-form-dataset-dataset_name"), "mock_10x_1234")

  dataset_names <- app$get_value(export = 'sc-form-dataset-dataset_names')
  is.int <- grep('mock_10x_integrated', dataset_names)[1]

  # clusters change after selecting new dataset
  current_clusters <- app$get_value(export = 'sc-form-sample_clusters-annot')
  app$set_inputs(`sc-form-dataset-selected_dataset` = as.character(is.int), wait_ = FALSE)
  app$wait_for_value(export = 'sc-form-sample_clusters-annot', timeout = timeout, ignore = list(current_clusters))

  integrated_clusters <- app$get_value(export = 'sc-form-sample_clusters-annot')
  expect_snapshot(integrated_clusters)

  expect_equal(app$get_value(export = "sc-form-dataset-dataset_name"), "mock_10x_integrated_harmony")

  # custom metric using either sample (batch) gives correct number of unique marker colors
  # fails if duplicate cell names not accounted for
  app$click("sc-form-gene_clusters-show_custom_metric")
  app$set_inputs(`sc-form-gene_clusters-custom_metric` = "batch=='mock_10x'")
  marker_colors <- app$wait_for_value(export = 'sc-marker_plot_cluster-colors', timeout = timeout, ignore = list(marker_colors))
  expect_setequal(table(marker_colors), c(mock10x_1234_ncells, mock10x_ncells))

  app$set_inputs(`sc-form-gene_clusters-custom_metric` = "batch=='mock_10x_1234'")
  marker_colors <- app$wait_for_value(export = 'sc-marker_plot_cluster-colors', timeout = timeout, ignore = list(marker_colors))
  expect_setequal(table(marker_colors), c(mock10x_1234_ncells, mock10x_ncells))

  # colors from unique expression values for EFNB1 equals unique colors in marker plot
  # fails if duplicate cell names not accounted for
  scseq <- load_scseq_qs('test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony', with_logs = TRUE)
  ft <- logcounts(scseq)['EFNB1', ]
  ft.scaled <- scales::rescale(ft)
  colors_expected <- get_expression_colors(ft.scaled)
  expect_ncolors <- length(unique(colors_expected))

  app$set_inputs(`sc-form-gene_clusters-custom_metric` = "EFNB1")
  marker_colors <- app$wait_for_value(export = 'sc-marker_plot_cluster-colors', timeout = timeout, ignore = list(marker_colors))
  ncolors <- length(unique(marker_colors))
  expect_equal(ncolors, expect_ncolors)
  expect_setequal(colors_expected, marker_colors)

})
