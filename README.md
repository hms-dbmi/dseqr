# drugseqr

<!-- badges: start -->
<!-- badges: end -->

The goal of drugseqr (*drug-seek-R*) is to find CMAP02/L1000 compounds that oppose RNA-seq gene expression signatures.

## Installation on Debian

[Salmon](https://combine-lab.github.io/salmon/) is required for quantifying the expression of transcripts using RNA-seq data. To install and make `salmon` command available:

```bash
# get salmon
cd /tmp
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.13.1/salmon-0.13.1_linux_x86_64.tar.gz

# extract
mkdir salmon
tar -xvf salmon-0.13.1_linux_x86_64.tar.gz -C salmon --strip-components=1

# make salmon command available
cd salmon
sudo mv bin/salmon /usr/local/bin
sudo mv lib/* /usr/local/bin

# check
salmon --version # salmon 0.13.1
```

To install drugseqr clone the repo then install from source code:

```R
install.packages('path/to/drugseqr', repos=NULL)
```

## Download CMAP02 and LINCS L1000 Data

The CMAP02 and LINCS L1000 data will be used to search for small molecules that oppose your query gene expression signature. To download this data (only has to be performed once):

```R
library(drugseqr)
dl_drug_es()
```

## Build Salmon Index and Ensembl Annotation Package

[Salmon](https://combine-lab.github.io/salmon/) requires an index built on a target transcriptome for quantification. For mapping counts from transcripts to genes, `drugseqr` requires an Ensembl annotation package. This step only has to be performed once.

To download a transcriptome, build an index, and build an Ensembl annotation package for homo sapiens (default):

```R
library(drugseqr)
build_index()
build_ensdb()
```

## Run Salmon, Load ExpressionSet, and Annotate

After building an index and ensembldb annotation package, you are ready to run salmon quantification and load/annotate the results. To do so:

```R
library(drugseqr)
# replace with path to folder with your raw fastq.gz files
data_dir <- system.file('extdata',  'IBD', package='drugseqr')

# replace with path to text file with sample annotations
# see pdata_path argument in ?run_salmon for specifications
pdata_path <- system.file('extdata',  'IBD', 'Phenotypes.csv', package='drugseqr')

# run transcript quantification using salmon
# this will load a GUI and prompt for various checks/annotations
run_salmon(data_dir, pdata_path)

# load and annotate RNA-seq quants
eset <- load_seq(data_dir)
```

## Run Differential Expression Analysis


After loading and annotating the raw RNA-seq data, you are ready to run differential expression analysis. To do so:

```R
library(drugseqr)

# replace with path to folder with your raw fastq.gz files
data_dir <- system.file('extdata',  'IBD', package='drugseqr')

# load eset saved from previous call to load_seq
eset <- readRDS(file.path(data_dir, 'eset.rds'))

# run differential expression analysis
# this will load a GUI and prompt you to select samples (rows) belonging to the control and test groups
anal <- diff_expr(eset, data_dir)
```

## Query CMAP02/L1000 Data

After running the differential expression analysis, you are ready to query against CMAP02 and L1000 signatures. To do so:

```R
library(drugseqr)

# load previous differential expression analysis
data_dir <- system.file('extdata', 'IBD', package='drugseqr')
anal <- readRDS(file.path(data_dir, 'diff_expr_symbol.rds'))

# extract effect size values used for query signature
dprimes <- get_dprimes(anal)

# load CMAP02 and L1000 datasets
cmap_es <- readRDS(system.file('extdata', 'cmap_es_ind.rds', package = 'drugseqr'))
l1000_es <- readRDS(system.file('extdata', 'l1000_es.rds', package = 'drugseqr'))

# run queries
cmap_res <- query_drugs(dprimes, cmap_es)
l1000_res <- query_drugs(dprimes, l1000_es)

# explore results
explore_results(cmap_res, l1000_res)
```


