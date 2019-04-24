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

## Build Salmon Index and Ensembl Annotation Package

Salmon requires an index built on a target transcriptome for quantification. For mapping counts from transcripts to genes, `drugseqr` requires an Ensembl annotation package. 

To download a transcriptome, build an index, and build an Ensembl annotation package for homo sapiens (default):

```R
library(drugseqr)
build_index()
build_ensdb()
```

## Run Salmon, Load ExpressionSet, and Annotate

After building and index and ensembldb annotation package, you are ready to run salmon quantification and load/annotate the results. To do so:

```R
# replace with path to folder with your raw fastq.gz files
data_dir <- system.file('extdata',  'IBD', package='drugseqr')

# replace with path to text file with sample annotations
# see pdata_path argument in ?load_seq for specifications
pdata_path <- system.file('extdata',  'IBD', 'Phenotypes.csv', package='drugseqr')

# run transcript quantification using salmon
run_salmon(data_dir)

# load and annotate RNA-seq quants
eset <- load_seq(data_dir, pdata_path)
```

