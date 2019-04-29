context("Pdata file matching")

test_that("match_pdata will find exact matches", {

  pdata <- data.frame(id=c('SRA_1', 'SRA_2'), stringsAsFactors = FALSE)
  fnames <- c('SRA_1', 'SRA_2')
  pdata <- match_pdata(pdata, fnames)

  expect_equal(pdata$id, pdata$`File Name`)
})

test_that("match_pdata adds file names in correct order", {

  pdata <- data.frame(id=c('SRA_1', 'SRA_2'), stringsAsFactors = FALSE)
  # fnames in reverse order
  fnames <- c('SRA_2', 'SRA_1')
  pdata <- match_pdata(pdata, fnames)

  expect_equal(pdata$id, pdata$`File Name`)
})

test_that("match_pdata can ignore '_' and '-' as word seperators", {

  pdata <- data.frame(id=c('SRA_1', 'SRA_2'))
  fnames <- c('SRA-1', 'SRA-2')
  expect_error(match_pdata(pdata, fnames), NA)

})

test_that("match_pdata ignores '.fastq.gz' suffix in file names while finding match", {

  pdata <- data.frame(id=c('SRA_1', 'SRA_2'), stringsAsFactors = FALSE)
  fnames <- c('SRA_1.fastq.gz', 'SRA_2.fastq.gz')
  pdata <- match_pdata(pdata, fnames)

  expect_equal(pdata$`File Name`, fnames)
})



