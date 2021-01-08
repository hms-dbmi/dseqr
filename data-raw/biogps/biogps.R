# This script annotates the BioGPS Human U133A Gene Atlas for use by explore_scseq_clusters
# the BioGPS data was downloaded from http://plugins.biogps.org/download/gnf1h-gcrma.zip

library(data.table)
library(crossmeta)

# load BioGPS Human U133A Gene Atlas
atlas <- fread('data-raw/biogps/U133AGNF1B.gcrma.avg.csv')
colnames(atlas)[1] <- 'PROBE'

# annotate to SYMBOL from U133A platform ----

# annotation for U133A
gse_name <- 'GSE1133'
data_dir <- 'data-raw/single-cell/example-anals'
eset <- load_raw(gse_name, data_dir)[[1]]
fdat <- fData(eset)[, c('PROBE', 'ENTREZID', 'SYMBOL')]

# join with atlas
table(fdat$PROBE %in% atlas$PROBE)
atlas <- atlas[PROBE %in% fdat$PROBE]
atlas <- atlas[fdat, on = 'PROBE']

# use max IQR to resolve duplicates by symbol
expr <- atlas[, -c('SYMBOL', 'ENTREZID', 'PROBE')]
fdat <- AnnotatedDataFrame(atlas[, .(SYMBOL, ENTREZID, PROBE)])
eset <- ExpressionSet(as.matrix(expr), featureData = fdat)

max_iqr <- which_max_iqr(eset, 'SYMBOL')
atlas <- atlas[max_iqr,  -'PROBE']

# used SYMBOL for key column
setkey(atlas, SYMBOL)

# save expression values
saveRDS(atlas, 'data-raw/biogps/biogps.rds')
