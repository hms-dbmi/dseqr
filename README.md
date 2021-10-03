## Dseqr
#### **End-to-End RNA-Seq Analysis**

Dseqr is a web application that helps you run 10X single-cell and bulk RNA-seq analyses from fastq → pathways → drug candidates.

💡 [Read the Docs and Deploy →](https://docs.dseqr.com)

### Local setup

```R
# install
install.packages('remotes')
remotes::install_github('hms-dbmi/dseqr')

# initialize and run new project
library(dseqr)
project_name <- 'example'

# directory to store application and project files in
data_dir <- './dseqr'

run_dseqr(project_name, data_dir)
```

If using fastq.gz files, build a `kallisto` index for quantification. To do so run:

```R
# default as used by run_dseqr
indices_dir <- file.path(data_dir, '.indices_dir')

rkal::build_kallisto_index(indices_dir)
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
