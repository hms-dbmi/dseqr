context("load_seq helpers")

test_that("setup_fdata expands entrez ids and makes characters", {

  # list of entrezids for single transcript
  enids <- c(5554, 11272, 100533464)
  tx2gene <- data.frame(tx_id = 'ENST00000228811',
                        gene_name = 'PRR4',
                        entrezid = I(list(ENST00000228811 = enids)))

  fdata <- setup_fdata(tx2gene)

  expect_equal(as.character(enids), fdata$ENTREZID_HS)
  expect_equal(class(fdata$ENTREZID_HS), 'character')

})

test_that("setup_fdata removes rows with duplicate gene_name and entrezid", {

  # list of entrezids for two different transcripts with same gene
  enids <- c(5554, 11272, 100533464)
  tx2gene <- data.frame(tx_id = c('ENST00000228811', 'ENST00000228812'),
                        gene_name = c('PRR4', 'PRR4'),
                        entrezid = I(list(ENST00000228811 = enids, ENST00000228812 = enids)))

  # just the first row
  tx2gene_row1 <- tx2gene[1,]

  fdata <- setup_fdata(tx2gene)
  fdata_row1 <- setup_fdata(tx2gene_row1)

  expect_equal(fdata, fdata_row1)
})


test_that("get_fastq_id1s returns sequence ids (starts with @)", {

  fastq_paths <- list.files(system.file('extdata', 'IBD', package='drugseqr'), '.fastq.gz$', full.names = TRUE)
  fastq_id1s <- get_fastq_id1s(fastq_paths)
  expect_true(all(grepl('^@', fastq_id1s)))

})

test_that("add_norms fails if it can't match file names with quant folders", {
  quants <- data.frame(file1 = 1, file2 = 2)
  pdata <- data.frame('File Name' = c('file-1.fastq.gz', 'file-2.fastq.gz'), check.names = FALSE, stringsAsFactors = FALSE)

  expect_error(add_norms(quants, pdata), 'failed to match')

})
