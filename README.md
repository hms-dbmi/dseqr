<!-- badges: start -->
[![CI](https://github.com/hms-dbmi/dseqr/actions/workflows/ci.yml/badge.svg)](https://github.com/hms-dbmi/dseqr/actions/workflows/ci.yml)
<!-- badges: end -->

## Dseqr
#### **End-to-End RNA-Seq Analysis**

Dseqr is a web application that helps you run 10X single-cell and bulk RNA-seq analyses from fastq → pathways → drug candidates.

💡 [Read the Docs and Open Dseqr →](https://docs.dseqr.com)


<h2></h2>
  <a href="https://docs.dseqr.com">
    <img src="https://user-images.githubusercontent.com/15719520/136054436-77ba2a23-1b0c-475e-a1d5-da5983edf2fd.gif"/>
  </a>
<h2></h2>
  

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

To enable bulk fastq.gz import, first build a `kallisto` index for quantification. To do so run:

```R
# default as used by run_dseqr
indices_dir <- file.path(data_dir, '.indices_dir')

rkal::build_kallisto_index(indices_dir)
```

### scRNA-seq fastqs
`dseqr` can directly import `cellranger` formatted count matrices. If you are starting
from fastq files, first install `kb-python`:

```console
# install kallisto|bustools wrapper (required)
pip install kb-python
```

Then run pseudo-quantification:

```R
# download pre-built index (mouse or human)
dseqr::download_kb_index(indices_dir, species = 'human')

# run pseudo-quantification
data_dir <- 'path/to/folder/with/fastqs'
dseqr::run_kb_scseq(indices_dir, data_dir, species = 'human')

# clean intermediate files produced by kb
dseqr::clean_kb_scseq(data_dir)
```

The resulting `cellranger` formatted count matrix files will be in the `data_dir`
subdirectory `bus_output/counts_unfiltered/cellranger`.


### Prefer docker?

```bash
# pull image
docker pull alexvpickering/dseqr

# run at http://0.0.0.0:3838/ and keep data on exit
docker run -v /full/path/to/data_dir:/srv/dseqr \
-p 3838:3838 \
alexvpickering/dseqr R -e 'library(dseqr); run_dseqr("example", "/srv/dseqr")'
```


### Host it

To spin up your own AWS infrastructure to host `dseqr`, see [dseqr.aws →](https://github.com/hms-dbmi/dseqr.aws)
