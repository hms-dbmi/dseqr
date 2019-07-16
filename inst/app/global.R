# things loaded in here are loaded once (even if multiple users)
library(drugseqr)
library(shiny)
library(shinyjs)
library(shinyWidgets)


# setup Drugs table annotation
cmap_annot <- get_drugs_table('CMAP02')
l1000_annot <- get_drugs_table('L1000')
