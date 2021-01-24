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
    R -e "remotes::install_github('rstudio/renv@0.12.5')"


COPY ./renv.lock .

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
RUN R -e "remotes::install_github('hms-dbmi/drugseqr@0.4.7', dependencies = FALSE, upgrade = FALSE)"
