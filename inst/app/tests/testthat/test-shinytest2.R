library(shinytest2)
library(SingleCellExperiment)

# create single-cell dataset for upload
mock_10x_files <- function() {

  set.seed(0)
  ngenes <- 10000
  sce <- scDblFinder::mockDoubletSCE(ncells = c(200, 300, 100, 200, 300), ngenes = ngenes)
  counts <- as(counts(sce), 'dgCMatrix')
  t2g <- dseqr.data::load_tx2gene()
  DropletUtils::write10xCounts(
    path = 'mock_10x',
    x = counts,
    gene.id = t2g$gene_id[seq_len(ngenes)],
    gene.symbol = t2g$gene_name[seq_len(ngenes)]
  )

  for (file in list.files('mock_10x', full.names = TRUE)) R.utils::gzip(file)
}


test_that("{shinytest2} recording: import_sc", {
  suppressWarnings(
  app <- AppDriver$new(
    name = "import_sc",
    options = list(shiny.maxRequestSize = 30*1024*1024^2),
    width = 1619,
    height = 909,
    seed = 0)
  )

  app$wait_for_idle()

  mock_10x_files()
  upload_dir <- 'mock_10x'
  on.exit(unlink(upload_dir, recursive = TRUE), add = TRUE)
  on.exit(unlink('test_data_dir', recursive = TRUE), add = TRUE)

  app$upload_file(`sc-form-dataset-up_raw` = file.path(upload_dir, "barcodes.tsv.gz"))
  app$upload_file(`sc-form-dataset-up_raw` = file.path(upload_dir, "genes.tsv.gz"))
  app$upload_file(`sc-form-dataset-up_raw` = file.path(upload_dir, "matrix.mtx.gz"))
  app$set_inputs(`sc-form-dataset-up_table_rows_selected` = 1:3, allow_no_input_binding_ = TRUE)
  app$set_inputs(`sc-form-dataset-sample_name` = "mock_10x")
  app$click("sc-form-dataset-add_sample")
  app$click("sc-form-dataset-import_datasets")
  app$wait_for_idle()

  app$click("sc-form-dataset-confirm_import_datasets")
  app$wait_for_value(export = 'sc-form-dataset-dataset_names', timeout = 1000*60*5, ignore = list(character(0)))
  app$set_inputs(`sc-form-dataset-selected_dataset` = "1", wait_ = FALSE)
  app$wait_for_idle()
  app$expect_values()

})
