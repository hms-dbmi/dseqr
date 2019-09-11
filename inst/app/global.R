# things loaded in here are loaded once (even if multiple users)
library(drugseqr)
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(dplyr)


# setup Drugs table annotation
# variable get updated when they are first needed
cmap_annot <- NULL
l1000_drugs_annot <- NULL
l1000_genes_annot <- NULL

# TODO make not hardcoded
pert_query_dir <- '/home/alex/Documents/Batcave/zaklab/drugseqr/data-raw/drug_gene_queries/data'
