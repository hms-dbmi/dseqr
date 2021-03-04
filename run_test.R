# for (f in list.files('R', '.R$', full.names = TRUE)) source(f)
# load('R/sysdata.rda')
# library(rlang)


# override default data_dir etc for development
app_name <- 'test'
app_dir <- 'inst/app'
data_dir <- '~/patient_data'
pert_query_dir <- '~/pert_query_dir'
pert_signature_dir <- '~/pert_signature_dir'
indices_dir <- '/srv/dseqr/indices'
run_dseqr(app_name, data_dir, app_dir, pert_query_dir, pert_signature_dir, indices_dir, test=TRUE, test_data=TRUE)
