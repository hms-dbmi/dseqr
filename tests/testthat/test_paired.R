context("detect if samples are paired")

test_that("detect_paired identifies SRA fastq files", {

  fastq_id1s <- c('@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=72',
                  '@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=72')

  expect_error(detect_paired(fastq_id1s), regexp = 'not yet implemented')

})


test_that("detect_paired handles older illumina sequence identifiers", {

  fastq_id1s <- c('@HWUSI-EAS100R:6:73:941:1973#0/1',
                  '@HWUSI-EAS100R:6:73:941:1973#0/1')

  expect_false(detect_paired(fastq_id1s))

})

test_that("detect_paired handles newer illumina sequence identifiers", {

  fastq_id1s <- c('@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG',
                  '@EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG')

  expect_true(detect_paired(fastq_id1s))

})

