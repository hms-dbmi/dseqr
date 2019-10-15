# drugseqr

<!-- badges: start -->
<!-- badges: end -->

The goal of drugseqr (*drug-seek-R*) is to find CMAP02/L1000 compounds that oppose RNA-seq gene expression signatures. The following instructions set up an Amazon EC2 spot instance to host and share the `drugseqr` web app.

## EC2 setup

Launch an Ubuntu Server 18.04 AMI instance with sufficient resources to meet your requirements. For example, I will launch a r4.large spot instance with a 50GiB SSD. I prefer to host a local copy of `drugseqr` to run all quantification and then transfer the saved data to the server. If you plan to upload raw RNA-Seq data to the server and run quantification there, you will likely need more resources.

Make sure that port 3838 is open to all inbound traffic so that the shiny server can be accessed.

## Setup the server

ssh into your instance and run [setup-aws.sh](scripts/setup-aws.sh). This will take 20 minutes or so:

```bash
curl -s https://e13eb58d90b1a6a62798c995485ad437be5e008f@raw.githubusercontent.com/hms-dbmi/drugseqr/master/scripts/setup-aws.sh | sudo bash
```

You should now be able to navigate your browser to  [EC2 Public DNS]:3838/drugseqr/example/ where EC2 Public DNS can be found in the EC2 instance description.


## Adding datasets

Use `rsync` to sync local and server data:

```bash
rsync -av --progress -e "ssh -i /path/to/mykeypair.pem" \
       ~/path/to/local/data_dir/ \ 
       ubuntu@[EC2 Public DNS]:/srv/shiny-server/drugseqr/example/data_dir/
```
