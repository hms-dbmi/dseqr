# drugseqr

<!-- badges: start -->
<!-- badges: end -->

The goal of drugseqr (*drug-seek-R*) is to find CMAP02/L1000 compounds that oppose RNA-seq gene expression signatures. The following instructions are to set up an Amazon EC2 spot instance to host and share the `drugseqr` web app.

## EC2 setup

Launch an Ubuntu Server 18.04 AMI instance with sufficient resources to meet your requirements. For example, I will launch a r4.large spot instance with a 50GiB SSD. I prefer to host a local copy of `drugseqr` to run all quantification and then transfer the saved data to the server. If you plan to upload raw RNA-Seq data to the server and run quantification there, you may need more resources.

Make sure that port 3838 is open to all inbound traffic so that the shiny server can be accessed.

## Kallisto installation

Kallisto is only required for quantification and, as such, may not be required on the server if quantification will be performed locally (see EC2 setup).

First get miniconda:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
sh Miniconda2-latest-Linux-x86_64.sh
```

Now exit the instance and re-ssh for changes to take effect.

[Kallisto](https://pachterlab.github.io/kallisto/download) is required for quantifying the expression of transcripts using RNA-seq data. To install `kallisto` with conda:

```bash
# setup bioconda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# install kallisto
conda install kallisto

# check
kallisto
```

## R and shiny server installation

Get R 3.6:

```bash
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
sudo apt install r-base r-base-dev
```

Get shiny and shiny-server:

```bash
sudo R
install.packages("shiny")
q('no')

sudo apt-get install gdebi-core
wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.9.923-amd64.deb
sudo gdebi shiny-server-1.5.9.923-amd64.deb
```

## drugseqr installation

To install drugseqr, first install system dependencies and `remotes`:

```bash
sudo apt install libcurl4-openssl-dev libv8-dev libxml2-dev libssl-dev
sudo R
install.packages('remotes')

# TODO: once public and release, remove access token and specify release
remotes::install_github("hms-dbmi/drugseqr", auth_token = 'e13eb58d90b1a6a62798c995485ad437be5e008f')
remotes::install_github("hms-dbmi/drugseqr.data", auth_token = 'e13eb58d90b1a6a62798c995485ad437be5e008f')
q('no')
# devtools::install_github("hms-dbmi/drugseqr@*release", auth_token = 'e13eb58d90b1a6a62798c995485ad437be5e008f')
```

## Download CMAP02 and LINCS L1000 data

The CMAP02 and LINCS L1000 data will be used to search for small molecules that oppose your query gene expression signature. To download this data into the `drugseqr.data` package:

```R
library(drugseqr.data)
dl_drug_es()
```

## Build kalisto index

[Kallisto](https://pachterlab.github.io/kallisto/) requires an index built on a target transcriptome for quantification. You may not need to do this step on the server if quantification will be run locally (see EC2 setup).

To download a transcriptome and build an index for homo sapiens (default):

```R
library(drugseqr.data)
build_kallisto_index()
```

## Setup the server

This sets up the file structure needed to initialize the example app:

```bash
sudo R
library(drugseqr)
init_drugseqr('example')
```

You should not be able to navigate your browser to  [EC2 Public DNS]:3838/drugseqr/example/ where EC2 Public DNS can be found in the EC2 instance description.


## Adding datasets
