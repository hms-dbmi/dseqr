library(dseqr)

gse_name <- 'GSE55347'
data_dir <- file.path('/mnt/shared', gse_name)

# load CMAP02 data
cmap_path <- system.file('extdata', 'cmap_es_ind.rds', package = 'dseqr.data', mustWork = TRUE)
cmap_es <- readRDS(cmap_path)

# load differential expression results
kal_version <- get_pkg_version('kallisto')
sal_version <- get_pkg_version('salmon')

kal_prefix <- paste('kallisto', kal_version, sep = '_')
sal_prefix <- paste('salmon', sal_version, sep = '_')

gse_drugs <- c('BETA-ESTRADIOL', 'BEZAFIBRATE', 'CLOFIBRIC ACID', 'CLOTRIMAZOLE', 'ECONAZOLE', 'GEMFIBROZIL', 'IFOSFAMIDE',
                'LEFLUNOMIDE', 'LOVASTATIN', 'MICONAZOLE', 'PIRINIXIC ACID', 'ROSIGLITAZONE', 'SIMVASTATIN')

load_query_genes <- function(drug, prefix) {
  anal_path <- file.path(data_dir, paste0('diff_expr_symbol_', prefix, '_', drug, '.rds'))
  anal <- readRDS(anal_path)

  dprimes <- get_dprimes(anal)
  return(dprimes)
}

kal_query_genes <- lapply(gse_drugs, load_query_genes, kal_prefix)
sal_query_genes <- lapply(gse_drugs, load_query_genes, sal_prefix)

# sort by decreasing similarity and take first instance of each drug (current app display)
get_query_ranks <- function(query_genes, drug_es, cmap_drug) {
  res <- query_drugs(query_genes, drug_es)
  res <- sort(res, decreasing = TRUE)
  drugs <- gsub('^([^_]+)_.+?$', '\\1', names(res))
  res <- res[!duplicated(drugs)]

  # minimum rank of correct drug
  grep(paste0('^', cmap_drug, '_'), names(res))
}

# equivalent drug names in cmap database
cmap_drugs <- c('estradiol', 'bezafibrate', 'clofibrate', 'clotrimazole', 'econazole', 'gemfibrozil', 'ifosfamide',
                'leflunomide', 'lovastatin', 'miconazole', 'pirinixic acid', 'rosiglitazone', 'simvastatin')


kal_ranks <- sapply(seq_along(kal_query_genes), function(i) get_query_ranks(kal_query_genes[[i]], cmap_es, cmap_drugs[i]))
sal_ranks <- sapply(seq_along(sal_query_genes), function(i) get_query_ranks(sal_query_genes[[i]], cmap_es, cmap_drugs[i]))

df <- data.frame(type = rep(c('salmon', 'kallisto'), each = length(sal_ranks)),
                 rank = c(sal_ranks, kal_ranks))

library(ggplot2)
ggplot(df, aes(x = type, y = rank)) + geom_boxplot() + geom_jitter(width = 0.15, height = 0.05)


# ROCR
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

kal_rates <- get_rates(kal_ranks)
kal_rates_df <- data.frame(kal_rates, Approach = "Kallisto")
MESS::auc(x = kal_rates_df$fpr, y = kal_rates_df$tpr)
# [1] 0.727123

sal_rates <- get_rates(sal_ranks)
sal_rates_df <- data.frame(sal_rates, Approach = "Salmon")
MESS::auc(x = sal_rates_df$fpr, y = sal_rates_df$tpr)
# [1] 0.6972477

# plot together
library(ggplot2)
rtdf <- rbind(kal_rates_df, sal_rates_df)
scaleFUN <- function(x) sprintf("%.1f", x)
rtpl <- ggplot(rtdf) +
  geom_line(aes(y = tpr, x = fpr, colour = Approach),
            size=1, show.legend = TRUE) +
  scale_colour_manual(values=c('#4582ec', '#d9534f')) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  ggtitle('ROCR for 13 CMAP Drugs Quantified with Kallisto and Salmon') +
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
