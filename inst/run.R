for (f in list.files('R', '.R$', full.names = TRUE)) source(f)
load('R/sysdata.rda')

args <- commandArgs(trailingOnly = TRUE)

project_name <- args[1]
is_example <- project_name == 'example'
if (project_name == 'private') project_name <- Sys.getenv('SHINYPROXY_USERNAME')

logout_url <- sprintf('https://%s/logout', args[2])

message('project_name: ', project_name)
message('is_example: ', is_example)
message('logout_url: ', logout_url)

run_dseqr(project_name,
          logout_url = logout_url,
          is_example = is_example)
