# things loaded in here are loaded once (even if multiple users)

suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(rlang)
})



# setup Drugs table annotation
# variable get updated when they are first needed
cmap_annot <- NULL
l1000_drugs_annot <- NULL
l1000_genes_annot <- NULL
