# things loaded in here are loaded once (even if multiple users)

suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(shinyBS)
  library(shinyWidgets)
  library(rlang)
  library(shinypanel)
  library(magrittr)
  dseqr_path <- find.package('dseqr')
  lazyLoad(file.path(dseqr_path, 'R/dseqr'))
  lazyLoad(file.path(dseqr_path, 'R/sysdata'))
})



# setup Drugs table annotation
# variable get updated when they are first needed
cmap_annot <- NULL
l1000_drugs_annot <- NULL
l1000_genes_annot <- NULL
