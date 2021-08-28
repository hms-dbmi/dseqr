dseqr_path <- find.package('dseqr')
invisible({
    lazyLoad(file.path(find.package('SingleCellExperiment'), 'R/SingleCellExperiment'))
    lazyLoad(file.path(dseqr_path, 'R/dseqr'))
    lazyLoad(file.path(dseqr_path, 'R/sysdata'))
})

args <- commandArgs(trailingOnly = TRUE)

app_name <- args[1]
is_example <- app_name == 'example'
logout_url <- sprintf('https://%s/logout', args[2])

message('app_name: ', app_name)
message('is_example: ', is_example)
message('logout_url: ', logout_url)

run_dseqr(app_name, logout_url = logout_url, is_example = is_example)
