context("Pdata file matching")

test_that("match_pdata will find exact matches", {

  pdata <- data.frame(id=c('SRA_1', 'SRA_2'))
  files <- c('SRA_1', 'SRA_2')
  pdata <- match_pdata(pdata, files)

  expect_equal(pdata$id, pdata$file)
})

test_that("match_pdata adds files in correct order", {

  pdata <- data.frame(id=c('SRA_1', 'SRA_2'))
  # files in reverse order
  files <- c('SRA_2', 'SRA_1')
  pdata <- match_pdata(pdata, files)

  expect_equal(pdata$id, pdata$file)
})

test_that("match_pdata can ignore '_' and '-' as word seperators", {

  pdata <- data.frame(id=c('SRA_1', 'SRA_2'))
  files <- c('SRA-1', 'SRA-2')
  expect_error(match_pdata(pdata, files), NA)

})

test_that("match_pdata ignores '.fastq.gz' suffix in files while finding match", {

  pdata <- data.frame(id=c('SRA_1', 'SRA_2'))
  files <- c('SRA_1.fastq.gz', 'SRA_2.fastq.gz')
  pdata <- match_pdata(pdata, files)

  expect_equal(pdata$file, files)
})



