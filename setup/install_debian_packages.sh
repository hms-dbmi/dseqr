#!/bin/bash

# Bash "strict mode", to help catch problems and bugs in the shell
# script. Every bash script you write should include this. See
# http://redsymbol.net/articles/unofficial-bash-strict-mode/ for
# details.
set -euo pipefail

export DEBIAN_FRONTEND=noninteractive

apt-get update

# Install required Debian packages
apt-get -y install --no-install-recommends "$@"

# Delete cached files we don't need anymore
apt-get clean
rm -rf /var/lib/apt/lists/*
