# based on https://www.huber.embl.de/users/klaus/stat_methods_bioinf/graphics_bioinf.html
# starts in drugseqr::diff_expr @ drugseqr::diff_anal

library(factoextra)
library(tibble)
library(MASS)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

create_dist_mat <- function(tdata, anno, method = "pearson"){

  dists <- as.matrix(get_dist(tdata, method = method))
  diag(dists) <- NA

  rownames(dists) <-  anno
  colnames(dists) <- anno

  return(dists)
}

# with sva ----

exprs_sva <- dups$exprs_sva
lib.size <- pData(eset)$lib.size * pData(eset)$norm.factors
exprs_sva <- exprs_sva[edgeR::filterByExpr(exprs_sva), ]
exprs_sva <- limma::voom(exprs_sva, modsv, lib.size)$E

ibd_log_counts <- t(exprs_sva)
ibd_log_counts[1:5,1:5]

dist_ibd <- get_dist(ibd_log_counts, method = "spearman")

scaling_ibd <- as_tibble(isoMDS(dist_ibd, k = 2)$points)
colnames(scaling_ibd) <- c("MDS_dimension_1", "MDS_dimension_2")


scaling_ibd <- add_column(scaling_ibd, cell_id = labels(dist_ibd), group = pData(eset)$group, .before = "MDS_dimension_1")


linecolors <- c("#714C02", "#01587A")
fillcolors <- c("#9D6C06", "#077DAA")

mds_plot_ibd <- ggplot(scaling_ibd, aes(x = MDS_dimension_1,
                                        MDS_dimension_2,
                                        color = group, fill = group)) +
  geom_point(size = 3, shape=21, alpha = 0.5, position=position_jitter(h=0.003, w=0.003)) +
  ggtitle("Kruskal MDS of the IBD data (with sva)") +
  scale_color_manual(values=linecolors) +
  scale_fill_manual(values=fillcolors) +
  coord_equal() +
  theme_bw()

mds_plot_ibd

hmcol <- colorRampPalette(brewer.pal(9, "RdPu"))(255)

pheatmap(create_dist_mat(ibd_log_counts,
                         anno = pData(eset)$group,
                         method = "euclidian"),
         trace="none", col = rev(hmcol))


# without sva -----

exprs_nosva <- exprs(eset)
lib.size <- pData(eset)$lib.size * pData(eset)$norm.factors
exprs_nosva <- limma::voom(exprs_nosva, modsv, lib.size)$E

ibd_log_counts <- t(exprs_nosva)
ibd_log_counts[1:5,1:5]

dist_ibd <- get_dist(ibd_log_counts, method = "spearman")

scaling_ibd <- as_tibble(isoMDS(dist_ibd, k = 2)$points)
colnames(scaling_ibd) <- c("MDS_dimension_1", "MDS_dimension_2")


scaling_ibd <- add_column(scaling_ibd, cell_id = labels(dist_ibd), group = pData(eset)$group, .before = "MDS_dimension_1")

linecolors <- c("#714C02", "#01587A")
fillcolors <- c("#9D6C06", "#077DAA")

mds_plot_ibd <- ggplot(scaling_ibd, aes(x = MDS_dimension_1,
                                        MDS_dimension_2,
                                        color = group, fill = group)) +
  geom_point(size = 3, shape=21, alpha = 0.5, position=position_jitter(h=0.003, w=0.003)) +
  ggtitle("Kruskal MDS of the IBD data (without sva)") +
  scale_color_manual(values=linecolors) +
  scale_fill_manual(values=fillcolors) +
  coord_equal() +
  theme_bw()

mds_plot_ibd

hmcol <- colorRampPalette(brewer.pal(9, "RdPu"))(255)

pheatmap(create_dist_mat(ibd_log_counts,
                         anno = pData(eset)$group,
                         method = "euclidian"),
         trace="none", col = rev(hmcol))


