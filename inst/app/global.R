# things loaded in here are loaded once (even if multiple users)

suppressPackageStartupMessages({
  library(drugseqr)
  library(shiny)
  library(shinyBS)
  library(shinydlplot)
  library(shinyjs)
  library(shinyWidgets)
  library(dplyr)
  # Seurat fails to load during cell-type deconvolution
})



# setup Drugs table annotation
# variable get updated when they are first needed
cmap_annot <- NULL
l1000_drugs_annot <- NULL
l1000_genes_annot <- NULL
