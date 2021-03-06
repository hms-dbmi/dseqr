## Dseqr
#### **End-to-End RNA-Seq Analysis**

Dseqr is a web application that helps you run 10X single-cell and bulk RNA-seq analyses from fastq → pathways → drug candidates.

💡 [Read the Docs and Deploy →](https://docs.dseqr.com)

### Local setup

```R
# install
install.packages('remotes')
remotes::install_github('hms-dbmi/dseqr')

# initialize and run new app
library(dseqr)
app_name <- 'example'
data_dir <- 'path/to/app_dir'
init_dseqr(app_name, data_dir)
run_dseqr(app_name, data_dir)
```

If using fastq.gz files, build a `kallisto` index for quantification. To do so run:

```R
rkal::build_kallisto_index('/srv/dseqr/indices')
```

### Prefer docker?

```bash
# pull image
docker pull alexvpickering/dseqr

# run at http://0.0.0.0:3838/ and keep data on exit
docker run -v path/to/app_dir:/srv/dseqr \
-p 3838:3838 \
alexvpickering/dseqr R -e 'dseqr::run_dseqr("example")'
```


### Host it

To spin up your own AWS infrastructure to host `dseqr`, see [dseqr.aws →](https://github.com/hms-dbmi/dseqr.aws)
