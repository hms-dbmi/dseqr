FROM rocker/r-ver:4.0.2

# install Ubuntu packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    pkg-config \
    libcurl4-openssl-dev \
    libv8-dev \
    libxml2-dev \
    libssl-dev \
    libbz2-dev \
    liblzma-dev \
    libhdf5-dev \
    zlib1g-dev \
    libpng-dev \
    libjpeg-dev \
    git \
    wget && rm -rf /var/lib/apt/lists/*


# install drugseqr dependencies from renv.lock file
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -e "remotes::install_github('rstudio/renv@0.9.2-20')"


COPY renv.lock .

# restore the package environment
RUN R -e 'options(renv.consent = TRUE); renv::restore()'

# Download miniconda and kallisto/bustools
# install in system-wide location
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh -O ~/miniconda.sh && \
    bash ~/miniconda.sh -b -p /opt/miniconda

# this sets path for current (root) user
ENV PATH="/opt/miniconda/bin:$PATH"

RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install kallisto=0.46.0 -y && \
    conda install -c bioconda bustools=0.39.3 -y

# download drug effect size data
RUN R -e "drugseqr.data::dl_drug_es()"

# install drugseqr last as will have to redo often
RUN R -e "remotes::install_github('hms-dbmi/drugseqr@0.1.76', dependencies = FALSE, upgrade = FALSE)"


# save image to a tar.gz file and upload to s3
# sudo docker build -t drugseqr .
# sudo docker save drugseqr:latest | gzip > drugseqr_latest.tar.gz
# aws s3 cp drugseqr_latest.tar.gz s3://drugseqr/drugseqr_latest.tar.gz
