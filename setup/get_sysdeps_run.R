install.packages('maketools')

lock <- renv:::renv_lockfile_read('renv.lock')
pkgs <- names(lock$Packages)

sysdeps_run <- NULL
for (pkg in pkgs) {
  deps <- maketools::package_sysdeps(pkg)
  if (!nrow(deps)) next

  sysdeps_run <- c(sysdeps_run, deps$package)
}

sysdeps_run <- unique(sysdeps_run)
writeLines(sysdeps_run, 'sysdeps_run.txt')

remove.packages('maketools')

