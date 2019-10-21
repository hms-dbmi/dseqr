# drugseqr

<!-- badges: start -->
<!-- badges: end -->

The goal of drugseqr (*drug-seek-R*) is to find CMAP02/L1000 compounds that oppose RNA-seq gene expression signatures. The following instructions set up an Amazon EC2 spot instance to host and share the `drugseqr` web app.

## EC2 setup

Launch an instance with sufficient resources to meet your requirements. For example, I will launch a r5.large spot instance with a 50GiB SSD. I prefer to host a local copy of `drugseqr` to run all quantification and then transfer the saved data to the server. If you plan to upload raw RNA-Seq data to the server and run quantification there, you will likely need more resources.

Make sure that port 80 is open to all inbound traffic so that the shiny server can be accessed.

## Setup the server


ssh into your instance and follow instructions to [install docker](https://docs.docker.com/install/).

Next, download the `drugseqr` image, load it, and initialize an empty `drugseqr` app:

```bash
# retrieve pre-built drugseqr docker image
wget https://drugseqr.s3.us-east-2.amazonaws.com/drugseqr_latest.tar.gz
sudo docker load < drugseqr_latest.tar.gz
rm drugseqr_latest.tar.gz

# permission change needed so that can init new app with user shiny inside the container
sudo mkdir -p /srv/shiny-server/
sudo chmod -R 0777 /srv/shiny-server/

# init new example app by running container
# host:container mounted volume in order to persist example app folders that are created inside the container
sudo docker run --user shiny --rm \
  -v /srv/shiny-server:/srv/shiny-server \
  drugseqr R -e "drugseqr::init_drugseqr('example')"
```


Then download example data and sync with previously initialized app:

```bash
wget https://drugseqr.s3.us-east-2.amazonaws.com/example_data.tar.gz
tar -xzvf example_data.tar.gz
rm example_data.tar.gz
sudo rsync -av example/ /srv/shiny-server/drugseqr/example/data_dir/
sudo chmod -R 0777 /srv/shiny-server/
sudo chmod -R 0777 /var/log/shiny-server/
```

Build kallisto index (optional - if will quantify bulk/sc fastq files on the server):

```bash
sudo docker run --user shiny --rm \
  -v /srv/shiny-server:/srv/shiny-server \
  drugseqr R -e "drugseqr.data::build_kallisto_index('/srv/shiny-server/drugseqr/indices')"
```

Now run a container to host the example app:

```bash
sudo docker run -d --user shiny --rm -p 80:3838 \
  -v /srv/shiny-server/:/srv/shiny-server/ \
  -v /var/log/shiny-server/:/var/log/shiny-server/ \
  drugseqr
```

You should now be able to navigate your browser to  [EC2 Public DNS]/drugseqr/example/ where EC2 Public DNS can be found in the EC2 instance description.


## Adding datasets

Use `rsync` to sync local and server data:

```bash
rsync -av --progress -e "ssh -i /path/to/mykeypair.pem" \
       ~/path/to/local/data_dir/ \ 
       ubuntu@[EC2 Public DNS]:/srv/shiny-server/drugseqr/example/data_dir/
```
