library(shinytest2)
library(SingleCellExperiment)
library(org.Hs.eg.db)

# create single-cell dataset for upload
mock_10x_files <- function() {

  set.seed(0)
  sce <- scDblFinder::mockDoubletSCE(ncells = c(200, 300, 100, 200, 300), ngenes = 10000)
  counts <- as(counts(sce), 'dgCMatrix')
  ensids <- keys(org.Hs.eg.db, keytype = 'ENSEMBL')[1:10000]
  genes <- mapIds(org.Hs.eg.db, ensids, column = 'SYMBOL', keytype = 'ENSEMBL')
  DropletUtils::write10xCounts(
    path = 'mock_10x',
    x = counts,
    gene.id = ensids,
    gene.symbol = genes
  )

  for (file in list.files('mock_10x', full.names = TRUE)) R.utils::gzip(file)
}


test_that("{shinytest2} recording: import_sc", {
  app <- AppDriver$new(
    name = "import_sc",
    options = list(shiny.maxRequestSize = 30*1024*1024^2),
    width = 1619,
    height = 909,
    seed = 0)

  app$wait_for_idle()

  mock_10x_files()
  upload_dir <- 'mock_10x'
  on.exit(unlink(uploads_dir, recursive = TRUE), add = TRUE)
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
