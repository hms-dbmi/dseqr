# install latest R
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
sudo apt install r-base r-base-dev

# install shiny
sudo R
install.packages("shiny")

# install latest shiny server
sudo apt-get install gdebi-core
wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.9.923-amd64.deb
sudo gdebi shiny-server-1.5.9.923-amd64.deb

# install app
git clone https://github.com/hms-dbmi/drugseqr.git

# resolve dependencies until works
R CMD build drugseqr
sudo R CMD INSTALL drugseqr_0.1.0.tar.gz


# move over app


