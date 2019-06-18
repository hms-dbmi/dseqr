#' Get empty droplets
#'
#' Used for calling empty droplets for testing kallisto and SoupX. Recommendation is currently to use alevin cell
#' filtering instead.
#'
#'
#' @param counts sparse dgTMatrix returned by \code{\link{load_kallisto}}.
#'
#' @return Boolean with length \code{ncol(counts)} indicating if the droplet is empty (\code{TRUE}) or a cell(\code{FALSE}).
#' @export
#'
#' @examples
get_empty <- function(counts) {
  out <- DropletUtils::emptyDrops(counts)

  # same thresholds as SoupX paper
  exp_genes <- tabulate(counts@j + 1)
  is.cell <- !(out$FDR > 0.05 | is.na(out$FDR) | exp_genes < 50)
  return(!is.cell)
}

#' Remove ambient gene expression from sc-RNAseq counts
#'
#' Uses \code{SoupX} to estimate and adjust counts for ambient contamination.
#' Currently only works with kallisto quantification results.
#'
#' @param counts sparse dgTMatrix returned by \code{\link{load_kallisto}}.
#' @param empty Boolean with indicating columns in \code{counts} that are empty droplets. Return by \code{\link{get_empty}}
#' @param project String. Passed to \code{CreateSeuratObject}.
#'
#' @return \code{Seurat} object with counts corrected for ambient contamination.
#' @export
#'
#' @examples
strain_scseq <- function(counts, empty, project) {

  hgGenes = c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE",
              "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")

  # construct soup
  toc <- counts[, !empty]
  scl <- SoupX::SoupChannelList(list(SoupX::SoupChannel(tod=counts, toc=toc, 'Channel1')))

  # infer useful candidate genes
  scl <- SoupX::inferNonExpressedGenes(scl)
  SoupX::plotCandidateMarkerGenes(scl, 'Channel1')

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
    return(Seurat::CreateSeuratObject(toc, project))
  }


  scl <- SoupX::calculateContaminationFraction(scl, "Channel1", nonExpressedGeneList, tgtSoupCntsPerGroup = 100)
  SoupX::plotChannelContamination(scl, "Channel1")

  # interpolate and adjust
  scl <- SoupX::interpolateCellContamination(scl, "Channel1")
  scl <- SoupX::adjustCounts(scl)

  strained_scseq <- Seurat::CreateSeuratObject(scl$atoc, project)
  return(strained_scseq)

}

#' Filter cells for scRNA-seq based on quality metrics
#'
#' This is used for filtering low quality cells from kallisto quantification after running \code{get_empty} and
#' either \code{strain_scseq} or just creating a Seurat object from the counts with empty droplets removed.
#'
#' Removes cells that are three median absolute deviations above or below the median for:
#' percentage mitochondrial counts, percentage ribosomal counts, number of expressed features, and
#' number of total counts.
#'
#' @param scseq Seurat object.
#'
#' @return \code{scseq} with low quality cells removed.
#' @export
#'
#' @examples
qc_scseq <- function(scseq) {
  # simple for now

  # ribosomal and mitochondrial genes
  rrna_path <- system.file('extdata', 'rrna.csv', package = 'drugseqr')
  mrna_path <- system.file('extdata', 'mrna.csv', package = 'drugseqr')

  rrna <- read.delim1(rrna_path)
  mrna <- read.delim1(mrna_path)

  rrna <- rrna[rrna %in% row.names(scseq)]
  mrna <- mrna[mrna %in% row.names(scseq)]

  scseq$percent.mito <- Seurat::PercentageFeatureSet(scseq, features = mrna)
  scseq$percent.ribo <- Seurat::PercentageFeatureSet(scseq, features = rrna)

  # exclude 3 mads below and above median for each
  scseq <- subset_scseq(scseq, 'percent.mito')
  scseq <- subset_scseq(scseq, 'percent.ribo')
  scseq <- subset_scseq(scseq, 'nFeature_RNA')
  scseq <- subset_scseq(scseq, 'nCount_RNA')
  return(scseq)

}

#' Remove cells above/below three mads from the median. Used by \code{qc_scseq}.
#'
#' @param scseq Seurat object
#' @param by String, a numeric column in \code{scseq@meta.data} to subset by.
#'
#' @return \code{scseq} with cells that are within three mads of the median of \code{by}.
#' @export
#' @keywords internal
#'
#' @examples
subset_scseq <- function(scseq, by) {
  # keep within 3 mads of median
  x <- scseq[[by]][[1]]
  med.x <- median(x)
  mad.x <- mad(x, med.x)
  keep <- which((x > med.x - 3*mad.x) & (x < med.x + 3*mad.x))
  return(scseq[, keep])
}
