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

# sets cache location for renv
ENV RENV_PATHS_ROOT=/src/.cache/R/renv

# install dseqr dependencies from renv.lock file
RUN touch /src/dseqr/.Renviron && \
    echo "RENV_PATHS_ROOT = $RENV_PATHS_ROOT" > /src/dseqr/.Renviron && \
    R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -e "remotes::install_github('rstudio/renv@0.14.0')"

# initial lockfile: sync periodically
COPY ./renv.lock.init .
RUN R -e 'renv::restore(lockfile="renv.lock.init")'

# lockfile: use this until slow
COPY ./renv.lock .
RUN R -e 'renv::restore(lockfile="renv.lock", clean = TRUE)'

# move packages to system library and clear out renv
RUN mv /src/.cache/R/renv/cache/v5/R-4.0/x86_64-pc-linux-gnu/* /usr/local/lib/R/library/ && \
    rm -rf /src && \
    mkdir -p /src/dseqr

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

COPY inst/run.R .
COPY R/* R/

# -------
FROM rhub/r-minimal:4.0.5
WORKDIR /src/dseqr

# get source code and R packages
COPY --from=0 /src /src
COPY --from=0 /usr/local/lib/R/library /usr/local/lib/R/library/

# set cache location for renv & add miniconda to PATH
ENV TMP_DIR=/srv/dseqr/tmp PATH="/src/miniconda/bin:$PATH"
RUN mkdir -p $TMP_DIR && \
    echo "TMPDIR = $TMP_DIR" > ${HOME}/.Renviron

# install dseqr last as will have to redo often
RUN R -e "remotes::install_github('hms-dbmi/dseqr@0.17.3', dependencies = FALSE, upgrade = FALSE)"

# docker build -t alexvpickering/dseqr:latest .
# docker push alexvpickering/dseqr
