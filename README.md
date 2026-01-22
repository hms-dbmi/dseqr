<!-- badges: start -->
[![CI](https://github.com/hms-dbmi/dseqr/actions/workflows/ci.yml/badge.svg)](https://github.com/hms-dbmi/dseqr/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/182834359.svg)](https://zenodo.org/badge/latestdoi/182834359)
<!-- badges: end -->

## Dseqr
#### **End-to-End RNA-Seq Analysis**

Dseqr is a web application that helps you run 10X single-cell and bulk RNA-seq analyses from fastq → pathways → drug candidates.

💡 [Read the Docs →](https://docs.dseqr.com)


### Local setup with Docker

```bash
# pull image
docker pull alexvpickering/dseqr --platform linux/amd64

# make directory to store data
mkdir dseqr_data

# run at http://0.0.0.0:3838/ and keep data on exit
docker run -v $(pwd)/dseqr_data:/srv/dseqr \
-p 3838:3838 \
alexvpickering/dseqr R -e 'library(dseqr); run_dseqr("example", "/srv/dseqr")'
```

<h2></h2>
  <a href="https://docs.dseqr.com">
    <img src="https://user-images.githubusercontent.com/15719520/136054436-77ba2a23-1b0c-475e-a1d5-da5983edf2fd.gif"/>
  </a>
<h2></h2>
  


