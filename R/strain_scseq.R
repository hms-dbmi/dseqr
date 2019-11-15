
#' Remove ambient gene expression from sc-RNAseq counts
#'
#' Uses \code{SoupX} to estimate and adjust counts for ambient contamination.
#' Currently only works with kallisto quantification results (and needs some work).
#'
#' @param counts sparse dgTMatrix of counts.
#' @param empty Boolean indicating columns in \code{counts} that are empty droplets.
#'
#' @return \code{counts} corrected for ambient contamination.
#' @export
strain_scseq <- function(counts, empty) {

   hgGenes = c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
   igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE",
               "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")
   nkGenes = c('GNLY', 'NKG7', 'GZMB')
   sfGenes = c('SFTPA2', 'SFTPA1', 'SFTPB')
   epGenes = c('KRT19', 'CAV1', 'EMP2')

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
  nkUseful <- nkGenes[nkGenes %in% isUseful]
  sfUseful <- sfGenes[sfGenes %in% isUseful]
  epUseful <- epGenes[epGenes %in% isUseful]

  nonExpressedGeneList <- list()
  if (length(hgUseful)) nonExpressedGeneList[['HG']] <- hgUseful
  if (length(igUseful)) nonExpressedGeneList[['IG']] <- igUseful
  if (length(nkUseful)) nonExpressedGeneList[['NK']] <- nkUseful
  if (length(sfUseful)) nonExpressedGeneList[['SF']] <- sfUseful
  if (length(epUseful)) nonExpressedGeneList[['EP']] <- epUseful

  if (!length(nonExpressedGeneList)) {
    message("Immunoglobulin and haemoglobin genes not useful. Failed to strain soup.")
    return(counts)
  }

  scl <- SoupX::calculateContaminationFraction(scl, "Channel1", nonExpressedGeneList)
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


