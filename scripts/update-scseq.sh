#!/bin/bash
# run by remote-update.sh

# get latest package and build/install
cd drugseqr
git pull https://$1:$2@github.com/hms-dbmi/drugseqr
cd ..
R CMD build drugseqr
sudo R CMD INSTALL drugseqr_0.1.0.tar.gz

# replace the shiny app
sudo rsync -a drugseqr/inst/app/ /srv/shiny-server/drugseqr/sjia/
sudo rsync -a drugseqr/inst/app/ /srv/shiny-server/drugseqr/amnon/

# permission to write to data directory
sudo chmod -R 0777 /srv/shiny-server/drugseqr/sjia/data_dir
sudo chmod -R 0777 /srv/shiny-server/drugseqr/amnon/data_dir

# restart the server
sudo systemctl restart shiny-server