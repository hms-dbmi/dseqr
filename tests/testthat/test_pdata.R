context("Pdata file matching")

test_that("match_pdata will find exact matches", {

  pdata <- data.frame(id=c('SRA_1', 'SRA_2'), stringsAsFactors = FALSE)
  qdirs <- c('SRA_1', 'SRA_2')
  pdata <- match_pdata(pdata, qdirs)

  expect_equal(pdata$id, pdata$quants_dir)
})

test_that("match_pdata adds qdirs in correct order", {

  pdata <- data.frame(id=c('SRA_1', 'SRA_2'), stringsAsFactors = FALSE)
  # qdirs in reverse order
  qdirs <- c('SRA_2', 'SRA_1')
  pdata <- match_pdata(pdata, qdirs)

  expect_equal(pdata$id, pdata$quants_dir)
})

test_that("match_pdata can ignore '_' and '-' as word seperators", {

  pdata <- data.frame(id=c('SRA_1', 'SRA_2'))
  qdirs <- c('SRA-1', 'SRA-2')
  expect_error(match_pdata(pdata, qdirs), NA)

})

test_that("match_pdata ignores '.fastq.gz' suffix in qdirs while finding match", {

  pdata <- data.frame(id=c('SRA_1', 'SRA_2'), stringsAsFactors = FALSE)
  qdirs <- c('SRA_1.fastq.gz', 'SRA_2.fastq.gz')
  pdata <- match_pdata(pdata, qdirs)

  expect_equal(pdata$quants_dir, qdirs)
})



