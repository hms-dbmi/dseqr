FROM rocker/shiny:3.6.1

# install Ubuntu packages
RUN apt-get update && apt-get install -y \
    sudo \
    libcurl4-openssl-dev \
    libv8-dev \
    libxml2-dev \
    libssl-dev \
    git \
    wget



# Download miniconda and kallisto/bustools
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh -O ~/miniconda.sh && \
    bash ~/miniconda.sh -b -p $HOME/miniconda

ENV PATH="/root/miniconda/bin:$PATH"


RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install kallisto=0.46.0 -y && \
    conda install -c bioconda bustools=0.39.3 -y


# install drugseqr
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@0.7.1-20')"
RUN R -e "remotes::install_github('hms-dbmi/drugseqr.data', auth_token = 'e13eb58d90b1a6a62798c995485ad437be5e008f')"

# clone the code base
RUN git clone https://e13eb58d90b1a6a62798c995485ad437be5e008f@github.com/hms-dbmi/drugseqr.git

# restore the package environment
RUN R -e 'setwd("./drugseqr"); renv::restore()'
