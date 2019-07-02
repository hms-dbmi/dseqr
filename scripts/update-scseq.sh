#!/bin/bash

# assumes on remote ssh
cd drugseqr
git pull
cd ..
R CMD build drugseqr
sudo R CMD INSTALL drugseqr_0.1.0.tar.gz

sudo rm -rf /srv/shiny-server/drugseqr/test/*
cp -r /usr/local/lib/R/site-library/drugseqr/shiny-apps/scseq/* /srv/shiny-server/drugseqr/test/

sudo systemctl restart shiny-server