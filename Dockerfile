FROM rocker/r-ver:4.4.0 AS build
WORKDIR /src/dseqr

# install required debian packages to install R packages
COPY setup/install_debian_packages.sh .
COPY setup/sysdeps_build_debian.txt .
RUN cat sysdeps_build_debian.txt | xargs ./install_debian_packages.sh

# add renv library to .libPaths
ENV RENV_LIB=/src/lib
RUN echo ".libPaths(c('$RENV_LIB', .libPaths()))" >> $(R RHOME)/etc/Rprofile.site

# install dseqr dependencies from renv.lock file
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -e "remotes::install_github('rstudio/renv@v1.0.7')" && \
    R -e "renv::init(bare = TRUE, settings = list(use.cache = FALSE))"

# initial lockfile: sync periodically
# delete renv cache
# strip debug from shared libraries
# see http://dirk.eddelbuettel.com/blog/2017/08/20/#010_stripping_shared_libraries
COPY ./renv.lock.init .
RUN R -e "renv::restore(lockfile='renv.lock.init', library = '$RENV_LIB')" && \
    R -e 'root <- renv::paths$root(); unlink(root, recursive = TRUE)' && \
    strip --strip-debug $RENV_LIB/*/libs/*.so

RUN R -e "renv::deactivate()"

# this sets path for current (root) user
ENV PATH="/src/miniconda/bin:$PATH"

# Download miniconda and kallisto/bustools
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh -O /src/miniconda.sh && \
  bash /src/miniconda.sh -b -p /src/miniconda && \
  conda config --add channels bioconda && \
  conda config --add channels conda-forge && \
  conda install kallisto=0.46.0 -y && \
  conda clean --force-pkgs-dirs -y

# current lockfile: use this until slow
COPY ./renv.lock .
RUN R -e "renv::restore(lockfile='renv.lock', library = '$RENV_LIB', clean = TRUE)" && \
    R -e 'root <- renv::paths$root(); unlink(root, recursive = TRUE)' && \
    strip --strip-debug $RENV_LIB/*/libs/*.so


# remove unecessary R packages
RUN R -e "remove.packages(c('remotes', 'renv'), .libPaths())"

# determine system run-time deps
COPY setup/get_sysdeps_run.R .
RUN Rscript get_sysdeps_run.R

# ----------
# COMMON
#-----------
FROM rocker/r-ver:4.4.0 AS common
WORKDIR /src/dseqr

# add conda to path
ENV PATH="/src/miniconda/bin:$PATH"

# get source code and R packages
COPY --from=build /src /src

# add renv library to .libPaths
ENV RENV_LIB=/src/lib
RUN echo ".libPaths(c('$RENV_LIB', .libPaths()))" >> $(R RHOME)/etc/Rprofile.site

# install runtime system deps
RUN cat sysdeps_run.txt | xargs ./install_debian_packages.sh

# set temporary directory for R
ENV TMP_DIR=/srv/dseqr/tmp
RUN mkdir -p $TMP_DIR && \
    echo "TMPDIR = $TMP_DIR" > ${HOME}/.Renviron

# install dseqr for callr::r_bg
# delete all except R/ directory
ADD R ./R
ADD inst ./inst
COPY DESCRIPTION NAMESPACE ./
RUN R -e "install.packages(repos=NULL, '.')" && \
    find . -maxdepth 1 ! -name R -exec rm -r "{}" \;

# ----------
# PRODUCTION
#-----------
from common AS production

# ----------
# TESTING
#-----------
from common AS testing

# google-chrome: to run shinytest2
# textlive-*: to build docs for rcmdcheck
RUN apt-get update && \
    apt-get -y install --no-install-recommends wget && \
    wget https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb && \
    apt-get -y install --no-install-recommends ./google-chrome-stable_current_amd64.deb && \
    rm ./google-chrome-stable_current_amd64.deb && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ADD tests ./tests
ADD man ./man
ADD inst ./inst
COPY DESCRIPTION NAMESPACE LICENSE setup/run_ci_tests.R ./

CMD ["Rscript", "run_ci_tests.R"]


# docker build -t alexvpickering/dseqr:latest .
