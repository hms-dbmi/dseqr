FROM rocker/shiny:3.6.1

# install Ubuntu packages
RUN apt-get update && apt-get install -y \
    sudo \
    libcurl4-openssl-dev \
    libv8-dev \
    libxml2-dev \
    libssl-dev \
    libbz2-dev \
    liblzma-dev \
    git \
    wget && rm -rf /var/lib/apt



# Download miniconda and kallisto/bustools
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh -O ~/miniconda.sh && \
    bash ~/miniconda.sh -b -p $HOME/miniconda

ENV PATH="/root/miniconda/bin:$PATH"


RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install kallisto=0.46.0 -y && \
    conda install -c bioconda bustools=0.39.3 -y


# install drugseqr dependencies from renv.lock file
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -e "remotes::install_github('rstudio/renv@0.7.1-20')"

COPY renv.lock .
COPY .Renviron .

# restore the package environment
RUN R -e 'options(renv.consent = TRUE); renv::restore()'

# clone the code base
RUN R -e "remotes::install_github('hms-dbmi/drugseqr@0.1.0', dependencies = FALSE, upgrade = FALSE)"

RUN R -e "drugseqr::init_drugseqr('example')"
