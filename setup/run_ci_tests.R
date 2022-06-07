Sys.setenv('NOT_CRAN' = 'true')
args <- c('--no-manual', '--no-build-vignettes')
rcmdcheck::rcmdcheck(error_on = 'error', build_args = args, args = args)
testthat::test_file('inst/app/tests/testthat/test-shinytest2.R')
