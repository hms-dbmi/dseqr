FROM rocker/r-ver:4.2.0 AS build
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
    R -e "remotes::install_github('rstudio/renv@0.15.5')" && \
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

# -------
FROM rocker/r-ver:4.2.0 AS deploy
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
# need dseqr in libPaths
ENV TMP_DIR=/srv/dseqr/tmp
RUN mkdir -p $TMP_DIR && \
    echo "TMPDIR = $TMP_DIR" > ${HOME}/.Renviron

# need dseqr installed for callr::r_bg
ADD R ./R
COPY inst/run.R .
ADD inst ./inst
COPY DESCRIPTION .
COPY NAMESPACE .
RUN R -e "install.packages(repos=NULL, '.')"

# docker build -t alexvpickering/dseqr:latest .
# docker push alexvpickering/dseqr
