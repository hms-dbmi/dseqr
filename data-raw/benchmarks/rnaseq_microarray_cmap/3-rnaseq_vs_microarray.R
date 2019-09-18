library(drugseqr)

# load RNA-Seq results
seq_name <- 'GSE55347'
mic_name <- 'GSE47875'
seq_dir <- file.path('/mnt/shared', seq_name)
mic_dir <- file.path('~/Documents/Batcave/zaklab/drugseqr/data-raw/benchmarks/rnaseq_microarray_cmap', mic_name)

# load CMAP02 data
cmap_path <- system.file('extdata', 'cmap_es_ind.rds', package = 'drugseqr.data', mustWork = TRUE)
cmap_es <- readRDS(cmap_path)

# load differential expression results
kal_version <- get_pkg_version('kallisto')
kal_prefix <- paste('kallisto', kal_version, sep = '_')
mic_prefix <- 'microarray'

gse_drugs <- c('BETA-ESTRADIOL', 'BEZAFIBRATE', 'CLOFIBRIC ACID', 'CLOTRIMAZOLE', 'ECONAZOLE', 'GEMFIBROZIL', 'IFOSFAMIDE',
               'LEFLUNOMIDE', 'LOVASTATIN', 'MICONAZOLE', 'PIRINIXIC ACID', 'ROSIGLITAZONE', 'SIMVASTATIN')

load_query_genes <- function(drug, prefix, data_dir) {
  anal_path <- file.path(data_dir, paste0('diff_expr_symbol_', prefix, '_', drug, '.rds'))
  if (!file.exists(anal_path)) browser()
  anal <- readRDS(anal_path)

  dprimes <- get_dprimes(anal)
  return(dprimes)
}

seq_query_genes <- lapply(gse_drugs, load_query_genes, kal_prefix, seq_dir)
mic_query_genes <- lapply(gse_drugs, load_query_genes, mic_prefix, mic_dir)

seq_query_res <- lapply(seq_query_genes, query_drugs, drug_es = cmap_es)
mic_query_res <- lapply(mic_query_genes, query_drugs, drug_es = cmap_es)

# equivalent drug names in cmap database
cmap_drugs <- c('estradiol', 'bezafibrate', 'clofibrate', 'clotrimazole', 'econazole', 'gemfibrozil', 'ifosfamide',
                'leflunomide', 'lovastatin', 'miconazole', 'pirinixic acid', 'rosiglitazone', 'simvastatin')

# get ranks ordered by similarity
get_query_ranks <- function(query_res, cmap_drug) {

  # sort result by decreasing similarity
  query_res <- sort(query_res, decreasing = TRUE)

  # get ranks of correct drug
  grep(paste0('^', cmap_drug, '_'),names(query_res))
}

seq_ranks <- lapply(seq_along(seq_query_res), function(i) get_query_ranks(seq_query_res[[i]], cmap_drugs[i]))
mic_ranks <- lapply(seq_along(mic_query_res), function(i) get_query_ranks(mic_query_res[[i]], cmap_drugs[i]))
seq_ranks <- unlist(seq_query_ranks)
mic_ranks <- unlist(mic_query_ranks)



# ROCR
get_rates <- function(res) {


  # for each posible position, add to tpr and fpr
  tpr <- c(0)
  fpr <- c(0)

  for (i in 1:ncol(cmap_es)) {

    # add to tpr num results that are correct at this position
    tpr <- c(tpr, tail(tpr, 1) + (sum(res == i)))

    # add to fpr num results that are incorrect at this position
    fpr <- c(fpr, tail(fpr, 1) + (sum(res != i)))
  }

  # devide by total number of true/false positives
  return(list(tpr = tpr / tail(tpr, 1), fpr = fpr / tail(fpr, 1)))
}


mic_rates <- get_rates(mic_ranks)
mic_rates_df <- data.frame(mic_rates, Approach = "Microarray")
MESS::auc(x = mic_rates_df$fpr, y = mic_rates_df$tpr)
# [1] 0.6242003

seq_rates <- get_rates(seq_ranks)
seq_rates_df <- data.frame(seq_rates, Approach = "RNA-Seq")
MESS::auc(x = seq_rates_df$fpr, y = seq_rates_df$tpr)
# [1] 0.6101993

# plot together
library(ggplot2)
rtdf <- rbind(mic_rates_df, seq_rates_df)
scaleFUN <- function(x) sprintf("%.1f", x)
rtpl <- ggplot(rtdf) +
  geom_line(aes(y = tpr, x = fpr, colour = Approach),
            size=1, show.legend = TRUE) +
  scale_colour_manual(values=c('#4582ec', '#d9534f')) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  ggtitle('ROCR for 13 CMAP Drugs Assayed with RNA-Seq and Microarray') +
  geom_abline(slope=1, colour = 'gray') +
  theme(
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_line(colour = "#dddddd", size = 0.6),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
    axis.text.x = element_text(margin=margin(5,5,10,5,"pt")),
    axis.text.y = element_text(margin=margin(5,5,10,5,"pt")),
    axis.title = element_text(size=16),
    axis.text = element_text(size=14))


rtpl

df <- data.frame(type = rep(c('RNA-seq', 'Microarray'), each = length(cmap_drugs)),
                 rank = c(seq_query_ranks, mic_query_ranks))
