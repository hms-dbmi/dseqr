# drugseqr

<!-- badges: start -->
<!-- badges: end -->

The goal of drugseqr (*drug-seek-R*) is to find CMAP02/L1000 compounds that oppose RNA-seq gene expression signatures. The following instructions set up an Amazon EC2 spot instance to host and share the `drugseqr` web app.

## EC2 setup

Launch an instance with sufficient resources to meet your requirements. For example, I will launch a r5.large spot instance with a 50GiB SSD. I prefer to host a local copy of `drugseqr` to run all quantification and then transfer the saved data to the server. If you plan to upload raw RNA-Seq data to the server and run quantification there, you will likely need more resources.

Make sure that port 8080 is open to all inbound traffic so that the web app can be accessed.

## Setup the server

The basic setup is going to be a docker container running ShinyProxy which will orchestrate starting docker containers running the app.

ssh into your instance and follow instructions to [install docker](https://docs.docker.com/install/). You likely also want to [configure](https://docs.docker.com/install/linux/linux-postinstall/#configure-docker-to-start-on-boot) docker to start on boot.

Next, create a docker network that ShinyProxy will use to communicate with the Shiny containers and build the ShinyProxy image. To do so, follow these [instructions](https://github.com/hms-dbmi/drugseqr.sp).

The `drugseqr` app won't work yet. To get it working, download the `drugseqr` image, load it, and initialize the example app:

```bash
# retrieve pre-built drugseqr docker image
wget https://drugseqr.s3.us-east-2.amazonaws.com/drugseqr_latest.tar.gz
sudo docker load < drugseqr_latest.tar.gz
rm drugseqr_latest.tar.gz


# we mount host:container volume in order to persist example app files/folders that are created inside the container
sudo docker run --rm \
  -v /srv/drugseqr:/srv/drugseqr \
  drugseqr R -e "drugseqr::init_drugseqr('example')"
```

Then download example data and sync with previously initialized app:

```bash
wget https://drugseqr.s3.us-east-2.amazonaws.com/example_data.tar.gz
tar -xzvf example_data.tar.gz
rm example_data.tar.gz
sudo rsync -av example/ /srv/drugseqr/example/
```

Build kallisto index (optional - if will quantify bulk/sc fastq files on the server):

```bash
sudo docker run --rm \
  -v /srv/drugseqr:/srv/drugseqr \
  drugseqr R -e "drugseqr.data::build_kallisto_index('/srv/drugseqr')"
```

## Run the app

Run a ShinyProxy container in detached mode `-d` and set the policy to always [restart](https://docs.docker.com/config/containers/start-containers-automatically/#use-a-restart-policy):

```bash
sudo docker run -d --restart always -v /var/run/docker.sock:/var/run/docker.sock --net sp-example-net -p 8080:8080 drugseqr.sp
```

You should now be able to navigate your browser to  [EC2 Public DNS]:8080/app/example where EC2 Public DNS can be found in the EC2 instance description.

## Adding single-cell datasets

Add single cell fastq.gz or CellRanger files to a directory in `single-cell`. For example, download CellRanger files:

```bash
# directory to store single cell sample data in  
cd /srv/drugseqr/example/single-cell
mkdir GSM2560249_pbmc_ifnb_full
cd GSM2560249_pbmc_ifnb_full

# get CellRanger files
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560249/suppl/GSM2560249%5F2%2E2%2Emtx%2Egz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560249/suppl/GSM2560249%5Fbarcodes%2Etsv%2Egz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl/GSE96583%5Fbatch2%2Egenes%2Etsv%2Egz

```

Run the app as before, create a new dataset, and select the created folder. Single cell 10X fastq.gz files can be added similarly. For example, using [bamtofastq](https://support.10xgenomics.com/docs/bamtofastq) to convert from 10XBAMs back to FASTQ:

```bash
# download bamtofastq into directory on the path
wget http://cf.10xgenomics.com/misc/bamtofastq -O ~/bin

# directory to store single cell sample data in  
cd /srv/drugseqr/example/single-cell
mkdir GSM3304014_lung_healthy
cd GSM3304014_lung_healthy

# download 10XBAM
wget https://sra-pub-src-1.s3.amazonaws.com/SRR7586091/P4_Normal_possorted_genome_bam.bam.1

# convert to fastq
bamtofastq P4_Normal_possorted_genome_bam.bam.1 ./fastqs
```

Run the app again, create a new dataset, and select the folder with the fastq files.

## Adding bulk datasets

Adding bulk datasets is similar to adding single-cell datasets. [GEOfastq](https://github.com/alexvpickering/GEOfastq) has a couple of utilities that make adding public bulk datasets particularly easy. For example, install `GEOfastq` then:

```R
# download bulk fastqs to appropriate directory for example app
data_dir <- '/srv/drugseqr/example/bulk'
gse_name <- 'GSE35296'

# first four samples for demonstration
srp_meta <- GEOfastq::get_srp_meta(gse_name, data_dir)
GEOfastq::get_fastqs(gse_name, srp_meta[1:4, ], data_dir)
```

## Sync local data with server

One way to this is is with `rsync`. For example:

```bash
rsync -av --progress -e "ssh -i /path/to/mykeypair.pem" \
       ~/path/to/local/example/ \ 
       ubuntu@[EC2 Public DNS]:/srv/drugseqr/example/
```
