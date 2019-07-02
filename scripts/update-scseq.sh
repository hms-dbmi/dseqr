#!/bin/bash
# run by remote-update.sh

# get latest package and build/install
cd drugseqr
git pull https://$1:$2@github.com/hms-dbmi/drugseqr
cd ..
R CMD build drugseqr
sudo R CMD INSTALL drugseqr_0.1.0.tar.gz

# replace the shiny app
sudo rm -rf /srv/shiny-server/drugseqr/test/*
cp -r /usr/local/lib/R/site-library/drugseqr/shiny-apps/scseq/* /srv/shiny-server/drugseqr/test/

# permission to write to data directory
chmod -R 0777 /srv/shiny-server/drugseqr/scseq/sjia/

# restart the server
sudo systemctl restart shiny-server