FROM rocker/r-ver:4.0.5

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

WORKDIR /src/dseqr

# install dseqr dependencies from renv.lock file
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -e "remotes::install_github('rstudio/renv@0.14.0')" && \
    R -e "renv::init(bare = TRUE, settings = list(use.cache = FALSE))"

# initial lockfile: sync periodically
COPY ./renv.lock.init .
RUN R -e 'renv::restore(lockfile="renv.lock.init")'

# Download miniconda and kallisto/bustools
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh -O /src/miniconda.sh && \
bash /src/miniconda.sh -b -p /src/miniconda

# this sets path for current (root) user
ENV PATH="/src/miniconda/bin:$PATH"

RUN conda config --add channels bioconda && \
conda config --add channels conda-forge && \
conda install kallisto=0.46.0 -y && \
conda install -c bioconda bustools=0.39.3 -y

# download drug effect size data
RUN R -e "dseqr.data::dl_data()"

# lockfile: use this until slow
COPY ./renv.lock .
RUN R -e 'renv::restore(lockfile="renv.lock", clean = TRUE)'

# move library and delete renv
RUN mv /src/dseqr/renv/library/R-4.0/x86_64-pc-linux-gnu /src/library && \
    rm -rf /src/dseqr && \
    mkdir /src/dseqr && \
    echo ".libPaths(c('/src/library', .libPaths()))" >> $R_HOME/etc/Rprofile.site

# -------
FROM rocker/r-ver:4.0.5
WORKDIR /src/dseqr

# get source code and R packages
COPY --from=0 /src /src
COPY --from=0 $R_HOME/etc/Rprofile.site $R_HOME/etc/Rprofile.site

# add conda to path
ENV PATH="/src/miniconda/bin:$PATH"

# set temporary directory for R
# need dseqr in libPaths
ENV TMP_DIR=/srv/dseqr/tmp
RUN mkdir -p $TMP_DIR && \
    echo "TMPDIR = $TMP_DIR" > ${HOME}/.Renviron && \
    apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev libhdf5-dev && \
    rm -rf /var/lib/apt/lists/* && \
    R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -e "remotes::install_github('hms-dbmi/dseqr@0.21.58', dependencies = FALSE, upgrade = FALSE)" && \
    R -e "remove.packages('remotes')"

# add source files
COPY inst/run.R .
COPY R/* R/

# docker build -t alexvpickering/dseqr:latest .
# docker push alexvpickering/dseqr
