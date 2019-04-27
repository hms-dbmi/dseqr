# this script is used to generate small example (first 1000 sequences) from fastq.gz files
# to run it requires the IBD fastq.gz files, Phenotypes.csv, and eset.rds generated from full IBD dataset
# example data is stored in inst/extdata/IBD (see http://r-pkgs.had.co.nz/data.html#data-extdata)

# full IBD fastq.gz files
fastq_files <- list.files('data-raw/example-data', '*.fastq.gz$')

extdata_dir <- file.path('inst', 'extdata', 'IBD')
dir.create(extdata_dir, recursive = TRUE)

for (inf in fastq_files) {

  # read 4 lines per sequence * 1000 sequences
  inf_path <- file.path('data-raw', 'example-data', inf)
  incon <- gzfile(inf_path)
  lines <- readLines(incon, n = 1000*4)

  # write to extdata
  outf <- gsub('IBD', 'IBD_1000', inf)
  outf_path <- file.path(extdata_dir, outf)
  outcon <- gzfile(outf_path, 'w')
  writeLines(lines, outcon)

  close(incon)
  close(outcon)
}

# also move over phenotype annotation and ExpressionSet from full IBD dataset
file.copy('data-raw/example-data/Phenotypes.csv', 'inst/extdata/IBD/Phenotypes.csv')
file.copy('data-raw/example-data/eset.rds', 'inst/extdata/IBD/eset.rds', overwrite = TRUE)
file.copy('data-raw/example-data/dgel.rds', 'inst/extdata/IBD/dgel.rds', overwrite = TRUE)
