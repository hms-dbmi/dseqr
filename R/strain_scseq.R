
#' Remove ambient gene expression from sc-RNAseq counts
#'
#' Uses \code{SoupX} to estimate and adjust counts for ambient contamination.
#' Currently only works with kallisto quantification results (and needs some work).
#'
#' @param counts sparse dgTMatrix returned by \code{\link{load_kallisto}}.
#' @param empty Boolean indicating columns in \code{counts} that are empty droplets. Return by \code{\link{get_empty}}
#'
#' @return \code{counts} corrected for ambient contamination.
#' @export
strain_scseq <- function(counts, empty) {

  hgGenes = c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE",
              "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")

  # construct soup
  toc <- counts[, !empty]
  scl <- SoupX::SoupChannelList(list(SoupX::SoupChannel(tod=counts, toc=toc, 'Channel1')))

  # infer useful candidate genes
  scl <- SoupX::inferNonExpressedGenes(scl)
  # SoupX::plotCandidateMarkerGenes(scl, 'Channel1')

  nonExpressedGenes <- scl$channels$Channel1$nonExpressedGenes
  isUseful <- row.names(nonExpressedGenes)[nonExpressedGenes$isUseful]

  # check if any in hg/ig subset
  hgUseful <- hgGenes[hgGenes %in% isUseful]
  igUseful <- igGenes[igGenes %in% isUseful]

  nonExpressedGeneList <- list()
  if (length(hgUseful)) nonExpressedGeneList[['HG']] <- hgUseful
  if (length(igUseful)) nonExpressedGeneList[['IG']] <- igUseful

  if (!length(nonExpressedGeneList)) {
    message("Immunoglobulin and haemoglobin genes not useful. Failed to strain soup.")
    return(counts)
  }

  scl <- SoupX::calculateContaminationFraction(scl, "Channel1", nonExpressedGeneList, tgtSoupCntsPerGroup = 1000)
  # SoupX::plotChannelContamination(scl, "Channel1")

  # interpolate and adjust
  scl <- SoupX::interpolateCellContamination(scl, "Channel1")
  scl <- SoupX::adjustCounts(scl)


  colnames(scl$atoc) <- gsub('^Channel1___', '', colnames(scl$atoc))
  return(scl$atoc)
}

get_empty <- function(counts) {

  out <- DropletUtils::emptyDrops(counts)
  empty <- out$FDR > 0.05 | is.na(out$FDR)

  return(empty)
}


