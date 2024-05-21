# update packages when R version updates
# assumes new R version is installed

# upgrade renv
renv::upgrade(reload = TRUE)

# start with clean renv package cache
# only keep renv
renv_lib <- grep('renv', .libPaths(), value = TRUE)[1]

del_pkgs <- list.files(renv_lib)
del_pkgs <- del_pkgs[del_pkgs != 'renv']

unlink(file.path(renv_lib, del_pkgs), recursive = TRUE)

# Bioconductor packages ----

# read lockfile records
lockfile_old <- renv:::renv_lockfile_read("renv.lock")
records <- lockfile_old$Packages

# Bioconductor packages tied to specific R version so they must be updated
# make sure we have latest version of BiocManager
renv::install("BiocManager")


# keep only packages from Bioconductor
bioc <- Filter(function(record) {
  record$Source == "Bioconductor"
}, records)

bioc <- names(bioc)

# install those packages and their dependencies
options(repos = BiocManager::repositories())
bioc.installed <- renv::install(bioc)

# snapshot without removing prior package records
renv::snapshot(packages = bioc, update = TRUE)

# Other packages ----

# read updated lockfile records
lockfile_new <- renv:::renv_lockfile_read("renv.lock")
records <- lockfile_new$Packages

# lockfile packages that were not updated above
all.pkgs <- names(records)
not.updated <- setdiff(all.pkgs, names(bioc.installed))

# restore those packages to recorded version
# update any that fail
done <- FALSE
updated <- NULL

while (!done) {

  tryCatch({
    renv::restore(packages = not.updated, prompt = FALSE)
    done <- TRUE

  }, error = function(e) {
    message(e$message)
    message('Updating failed package')
    failed.pkg <- gsub("^install of package '(.+?)' failed .+?$", "\\1", e$message)
    updated <<- c(updated, failed.pkg)

    renv::install(failed.pkg, prompt = FALSE)
    renv::snapshot(packages = failed.pkg, update = TRUE, prompt = FALSE)
  })
}

if (length(updated)) {
  warning('The following packages were updated:\n', paste(updated, collapse = '\n'))
}

# should be up to date
renv::status()
