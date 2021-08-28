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
    R -e "remotes::install_github('rstudio/renv@0.14.0')"

# sets cache location for renv
ENV RENV_PATHS_ROOT=/src/.cache/R/renv

# initial lockfile: sync periodically
COPY ./renv.lock.init .
RUN R -e 'renv::restore(lockfile="renv.lock.init")'

# lockfile: use this until slow
COPY ./renv.lock .
RUN R -e 'renv::restore(lockfile="renv.lock", clean = TRUE)'

RUN rm -rf $RENV_PATHS_ROOT/source

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

COPY --from=0 /src/* /src/
WORKDIR /src/dseqr

# set tmp directory on EFS (for file uploads)
ENV TMP_DIR=/srv/dseqr/tmp

# install dseqr last as will have to redo often
RUN R -e "renv::install('hms-dbmi/dseqr@0.17.1')" && \
    mkdir -p $TMP_DIR && \
    echo "TMPDIR = $TMP_DIR" > ${HOME}/.Renviron

# docker build -t alexvpickering/dseqr:latest .
# docker push alexvpickering/dseqr
