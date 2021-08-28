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
          logout_url = logout_url,
          is_example = is_example)
