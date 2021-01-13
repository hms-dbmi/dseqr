#! /bin/bash
git merge --no-commit --no-ff master
git reset HEAD .gitignore
git reset HEAD .dockerignore
git reset HEAD .renvignore
git reset HEAD Dockerfile
git reset HEAD renv/.gitignore
git reset HEAD renv/activate.R
git reset HEAD renv/settings.dcf
git reset HEAD .Rprofile
git reset HEAD renv.lock
git restore .
