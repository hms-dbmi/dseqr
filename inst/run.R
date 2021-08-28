
for (f in list.files('R', '.R$', full.names = TRUE)) source(f)
load('R/sysdata.rda')

args <- commandArgs(trailingOnly = TRUE)

app_name <- args[1]
is_example <- app_name == 'example'
logout_url <- sprintf('https://%s/logout', args[2])

message('app_name: ', app_name)
message('is_example: ', is_example)
message('logout_url: ', logout_url)

run_dseqr(app_name,
          app_dir = 'inst/app',
          data_dir = '~/patient_data',
          pert_query_dir = '~/dseqr/pert_query_dir',
          pert_signature_dir = '~/dseqr/pert_signature_dir',
          gs_dir = '~/dseqr/gs_dir',
          indices_dir = '~/dseqr/indices',
          tx2gene_dir = '~/dseqr/tx2gene',
          logout_url = logout_url,
          is_example = is_example)
