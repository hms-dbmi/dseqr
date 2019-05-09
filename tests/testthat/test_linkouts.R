context("linkouts can be added and merged")
library(tibble)

query_res <- tibble('Pubchem CID' = c('1234', '5678', NA),
                    'SIDER' = c('1234', NA, '5678'))

pub_img <- 'https://pubchem.ncbi.nlm.nih.gov/pcfe/favicon/favicon.ico'
pub_pre <- 'https://pubchem.ncbi.nlm.nih.gov/compound/'

sidder_img <- 'http://sideeffects.embl.de/media/images/EMBL_Logo.png'
sidder_pre <- 'http://sideeffects.embl.de/drugs/'

with_linkouts <- add_linkout(query_res, 'Pubchem CID', pub_img, pub_pre, '/', title = 'Pubchem')
with_linkouts <- add_linkout(with_linkouts, 'SIDER', sidder_img, sidder_pre)


test_that("add_linkout skips NA entries", {
  expect_true(is.na(with_linkout$`Pubchem CID`[3]))
})

test_that("add_linkout uses pre and post url", {
  expect_match(with_linkout$`Pubchem CID`[1], 'https://pubchem.ncbi.nlm.nih.gov/compound/1234/')
})

test_that("add_linkout adds title with 'Go to' prepended", {
  expect_match(with_linkout$`Pubchem CID`[1], 'title="Go to Pubchem"')
})

test_that("add_linkout uses column name as default title", {
  default_title <- add_linkout(query_res,
                              'Pubchem CID',
                              'https://pubchem.ncbi.nlm.nih.gov/pcfe/favicon/favicon.ico',
                              'https://pubchem.ncbi.nlm.nih.gov/compound/',
                              '/post_stuff')

  expect_match(default_title$`Pubchem CID`[1], 'title="Go to Pubchem CID"')
})

test_that("merge_linkouts can merge non-NA linkouts from multiple columns", {
  with_external <- merge_linkouts(with_linkouts, c('Pubchem CID', 'SIDER'))

  # second entry only has Pubchem
  expect_equal(with_external$`External Links`[2], with_linkouts$`Pubchem CID`[2])

  # third entry only has SIDER
  expect_equal(with_external$`External Links`[3], with_linkouts$SIDER[3])

  # first entry has both
  expect_equal(with_external$`External Links`[1],
               paste0(with_linkouts$`Pubchem CID`[1], with_linkouts$SIDER[1]))
})
