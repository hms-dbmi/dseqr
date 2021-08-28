# things loaded in here are loaded once (even if multiple users)

suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(shinyBS)
  library(shinyWidgets)
  library(rlang)
  library(shinypanel)
  library(magrittr)
})

# attach if loaded otherwise lazy
if (isNamespaceLoaded('dseqr')) {
  require(dseqr)

} else {
  print('lazy loading dseqr')
  pdir <- find.package('dseqr')
  lazyLoad(file.path(pdir, 'R/dseqr'))
  lazyLoad(file.path(pdir, 'R/sysdata'))
}

# attach if loaded otherwise lazy
if (isNamespaceLoaded('SingleCellExperiment')) {
  require(SingleCellExperiment)

} else {
  print('lazy loading SingleCellExperiment')
  pdir <- find.package('SingleCellExperiment')
  lazyLoad(file.path(pdir, 'R/SingleCellExperiment'))
  lazyLoad(file.path(pdir, 'R/SingleCellExperiment'))
}

# setup Drugs table annotation
# variable get updated when they are first needed
cmap_annot <- NULL
l1000_drugs_annot <- NULL
l1000_genes_annot <- NULL
