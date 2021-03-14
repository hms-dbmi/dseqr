# this script moves everything from rds to qs
# it is necessary to transition data created prior to version 0.10.0

# fill in path to app data
app_dir <- "/home/alex/patient_data/example"

# get all .rds files
rds_files <- list.files(app_dir, '.rds$', full.names = TRUE, recursive = TRUE)

for (rds_file in rds_files) {
  item <- readRDS(rds_file)
  qs_file <- gsub('.rds$', '.qs', rds_file)
  qs::qsave(item, qs_file)
  unlink(rds_file)
}
