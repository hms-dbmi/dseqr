#!/bin/bash

# get miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
echo 'export PATH=$HOME/miniconda/bin:$PATH' >> ~/.bashrc

# add miniconda to PATH
source ~/.bashrc

# add miniconda to PATH for Rstudio session (local)
echo 'PATH=${HOME}/miniconda/bin:${PATH}' >> ~/.Renviron

# setup bioconda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# install kallisto
conda install kallisto=0.46.0 -y
conda install -c bioconda bustools=0.39.3 -y

# get R 3.6
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
sudo DEBIAN_FRONTEND=noninteractive apt install -y r-base r-base-dev

# install shiny and shiny-server
echo "install.packages(\"shiny\", repos=\"https://cran.rstudio.com\")" | sudo R --no-save

sudo apt install gdebi-core -y
wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.9.923-amd64.deb
sudo gdebi shiny-server-1.5.9.923-amd64.deb --n

# install drugseqr
sudo apt install libcurl4-openssl-dev libv8-dev libxml2-dev libssl-dev -y
echo "install.packages(\"remotes\", repos=\"https://cran.rstudio.com\")" | sudo R --no-save

# this takes a while to download drugseqr dependencies
# TODO: once public and release, remove access token and specify release
echo "remotes::install_github(\"hms-dbmi/drugseqr\", auth_token = \"e13eb58d90b1a6a62798c995485ad437be5e008f\", upgrade = FALSE)" | sudo R --no-save
echo "remotes::install_github(\"hms-dbmi/drugseqr.data\", auth_token = \"e13eb58d90b1a6a62798c995485ad437be5e008f\")" | sudo R --no-save

# download CMAP02/LINCS data and build kallisto index
echo "drugseqr.data::dl_drug_es()" | sudo R --no-save
echo "drugseqr.data::build_kallisto_index()" | sudo R --no-save

# setup example app
echo "drugseqr::init_drugseqr(\"example\")" | sudo R --no-save
