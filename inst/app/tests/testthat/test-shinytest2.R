suppressPackageStartupMessages({
  library(shinytest2)
  library(SingleCellExperiment)
})

# create single-cell dataset for upload
mock_10x_files <- function(dir_name) {

  set.seed(0)
  ngenes <- 10000
  sce <- scDblFinder::mockDoubletSCE(ncells = c(200, 300, 100, 200, 300), ngenes = ngenes)
  counts <- as(counts(sce), 'dgCMatrix')
  t2g <- dseqr.data::load_tx2gene()
  DropletUtils::write10xCounts(
    path = dir_name,
    x = counts,
    gene.id = t2g$gene_id[seq_len(ngenes)],
    gene.symbol = t2g$gene_name[seq_len(ngenes)]
  )

  for (file in list.files(dir_name, full.names = TRUE)) R.utils::gzip(file)
}

# 10 mins
timeout <- 1000*60*10

test_that("{shinytest2} recording: Single-Cell Tab", {
  suppressWarnings(
    app <- AppDriver$new(
      name = "import_sc",
      options = list(shiny.maxRequestSize = 30*1024*1024^2),
      width = 1619,
      height = 909,
      seed = 0)
  )

  app$wait_for_idle()

  sample_dir <- 'mock_10x'
  data_dir <- 'test_data_dir'

  list_files <- function()
    file.path(data_dir, list.files(data_dir, recursive = TRUE, all.files = TRUE, include.dirs = TRUE))

  # created expected files/folders
  init_files <- list_files()
  expect_snapshot(init_files)

  mock_10x_files(sample_dir)
  on.exit(unlink(sample_dir, recursive = TRUE), add = TRUE)
  on.exit(unlink(data_dir, recursive = TRUE), add = TRUE)


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
    as.character(1:6))

  # saved prev_dataset
  prev_file <- setdiff(list_files(), c(init_files, dataset_files))
  expect_equal(qs::qread(prev_file), sample_dir)

  # didn't auto-select cluster
  expect_equal(app$get_value(input = 'sc-form-cluster-selected_cluster'), '')
  app$set_inputs(`sc-form-cluster-selected_cluster` = "1")

  # created markers file
  app$wait_for_value(export = 'sc-form-cluster-have_selected_markers', ignore = list(FALSE))
  all_prev_files <- c(init_files, dataset_files, prev_file)
  markers_file <- setdiff(list_files(), all_prev_files)
  expect_equal(markers_file, "test_data_dir/test_user/default/single-cell/mock_10x/snn1/markers.qs")
  expect_equal(unname(tools::md5sum(markers_file)), "3f7e85bde918f9c04a73af1804986414")

  # clusters and markers plot have same coordinates
  app$set_inputs(`sc-form-gene_clusters-gene_table_rows_selected` = 1, allow_no_input_binding_ = TRUE)
  cluster_coords <- jsonlite::fromJSON(app$get_value(output = 'sc-cluster_plot-cluster_plot'))$x$coords
  markers_coords <- jsonlite::fromJSON(app$get_value(output = 'sc-marker_plot_cluster-marker_plot'))$x$coords

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
    c("CD14 Mono", as.character(2:6)))

  # TODO: add input values to exports
  # can compare one cluster vs one other cluster
  app$click("sc-form-cluster-show_contrasts")
  app$get_value(input = "sc-form-dataset-selected_dataset")
  app$click("sc-form-cluster-show_contrasts")

  # can set custom boolean metric
  app$click("sc-form-gene_clusters-show_custom_metric")
  app$set_inputs(`sc-form-gene_clusters-custom_metric` = "KCNT2>0")
  colors <- jsonlite::fromJSON(app$get_value(output = 'sc-marker_plot_cluster-marker_plot'))$x$colors
  expect_length(unique(colors), 2)

  # can view expression of gene using custom metric
  app$set_inputs(`sc-form-gene_clusters-custom_metric` = "KCNT2")
  colors <- jsonlite::fromJSON(app$get_value(output = 'sc-marker_plot_cluster-marker_plot'))$x$colors
  expect_length(unique(colors), 163)

  # can use cluster for metric
  app$set_inputs(`sc-form-gene_clusters-custom_metric` = "cluster==1")
  colors <- jsonlite::fromJSON(app$get_value(output = 'sc-marker_plot_cluster-marker_plot'))$x$colors
  expect_length(unique(colors), 2)

  # TODO: check that it was added to table
  # can add custom metric
  prev_files <- list_files()
  app$click("sc-form-gene_clusters-save_custom_metric")
  expect_equal(app$get_value(input = "sc-form-gene_clusters-custom_metric"), "")
  saved_metric_files <- setdiff(list_files(), prev_files)
  expect_snapshot(saved_metric_files)

  # TODO: check that closing custom metric input hides markers plot
  app$click("sc-form-gene_clusters-show_custom_metric")

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

  # run label transfer from new to current
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

  # clusters change after selecting new dataset
  current_clusters <- app$get_value(export = 'sc-form-sample_clusters-annot')
  app$set_inputs(`sc-form-dataset-selected_dataset` = "3", wait_ = FALSE)
  app$wait_for_value(export = 'sc-form-sample_clusters-annot', timeout = timeout, ignore = list(current_clusters))

  integrated_clusters <- app$get_value(export = 'sc-form-sample_clusters-annot')
  expect_snapshot(integrated_clusters)

  expect_equal(app$get_value(export = "sc-form-dataset-dataset_name"), "mock_10x_integrated_harmony")
})
