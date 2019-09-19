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

# sort by decreasing similarity and take first instance of each drug (current app display)
get_query_ranks <- function(query_genes, drug_es, cmap_drug, sort_by = c('min', 'avg')) {
  res <- query_drugs(query_genes, drug_es)
  res <- sort(res, decreasing = TRUE)
  drugs <- gsub('^([^_]+)_.+?$', '\\1', names(res))

  if (sort_by[1] == 'min') {
    rank <- which(drugs[!duplicated(drugs)] == cmap_drug)

  } else if (sort_by[1] == 'avg') {
    tb <- tibble::tibble(cor = res, drug = drugs)
    avg_drugs <- tb %>%
      dplyr::group_by(drug) %>%
      dplyr::summarise(mean_cor = mean(cor)) %>%
      dplyr::arrange(-mean_cor) %>%
      dplyr::pull(drug)

    rank <- which(avg_drugs == cmap_drug)
  }

  return(rank)
}

# equivalent drug names in cmap database
cmap_drugs <- c('estradiol', 'bezafibrate', 'clofibrate', 'clotrimazole', 'econazole', 'gemfibrozil', 'ifosfamide',
                'leflunomide', 'lovastatin', 'miconazole', 'pirinixic acid', 'rosiglitazone', 'simvastatin')

seq_ranks <- lapply(seq_along(seq_query_genes), function(i) get_query_ranks(seq_query_genes[[i]], cmap_es, cmap_drugs[i]))
mic_ranks <- lapply(seq_along(mic_query_genes), function(i) get_query_ranks(mic_query_genes[[i]], cmap_es, cmap_drugs[i]))

seq_ranks_avg <- lapply(seq_along(seq_query_genes), function(i) get_query_ranks(seq_query_genes[[i]], cmap_es, cmap_drugs[i], sort_by = 'avg'))
mic_ranks_avg <- lapply(seq_along(mic_query_genes), function(i) get_query_ranks(mic_query_genes[[i]], cmap_es, cmap_drugs[i], sort_by = 'avg'))

# ROCR seq vs mic
get_rates <- function(res) {


  # for each posible position, add to tpr and fpr
  tpr <- c(0)
  fpr <- c(0)

  for (i in 1:1309) {

    # add to tpr num results that are correct at this position
    tpr <- c(tpr, tail(tpr, 1) + (sum(res == i)))

    # add to fpr num results that are incorrect at this position
    fpr <- c(fpr, tail(fpr, 1) + (sum(res != i)))
  }

  # devide by total number of true/false positives
  return(list(tpr = tpr / tail(tpr, 1), fpr = fpr / tail(fpr, 1)))
}

mic_rates <- get_rates(mic_ranks)
mic_rates_df <- data.frame(mic_rates, Approach = "Microarray (min)")
MESS::auc(x = mic_rates_df$fpr, y = mic_rates_df$tpr)
# [1] 0.7592919

mic_rates_avg <- get_rates(mic_ranks_avg)
mic_rates_avg_df <- data.frame(mic_rates_avg, Approach = "Microarray (avg)")
MESS::auc(x = mic_rates_avg_df$fpr, y = mic_rates_avg_df$tpr)
# [1] 0.6724888

seq_rates <- get_rates(seq_ranks)
seq_rates_df <- data.frame(seq_rates, Approach = "RNA-Seq (min)")
MESS::auc(x = seq_rates_df$fpr, y = seq_rates_df$tpr)
# [1] 0.727123

seq_rates_avg <- get_rates(seq_ranks_avg)
seq_rates_avg_df <- data.frame(seq_rates_avg, Approach = "RNA-Seq (avg)")
MESS::auc(x = seq_rates_avg_df$fpr, y = seq_rates_avg_df$tpr)
# [1] 0.6638438

# plot together
library(ggplot2)
rtdf <- rbind(mic_rates_df, seq_rates_df)
rtdf <- rbind(mic_rates_avg_df, mic_rates_df)
rtdf <- rbind(seq_rates_avg_df, seq_rates_df)

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
