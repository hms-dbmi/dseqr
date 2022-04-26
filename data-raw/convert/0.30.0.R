# add default project dirs
# remove sample metrics from cdata

data_dir <- '/home/alex/Documents/Data/dseqr'
app_dirs <- list.files(data_dir)
app_dirs <- setdiff(app_dirs,
                    c('indices', 'node_modules', 'pert_signature_dir',
                      'pert_query_dir', 'gs_dir', 'example_data.tar.gz', 'tx2gene'))


for  (app_dir in app_dirs) {
  # get dirs
  move_dirs <- list.files(file.path(data_dir, app_dir))

  # make default
  default_dir <- file.path(data_dir, app_dir, 'default')
  if (dir.exists(default_dir)) next()
  dir.create(default_dir)

  file.move(from = file.path(data_dir, app_dir, move_dirs),
            to = file.path(default_dir, move_dirs))
}
