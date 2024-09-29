#!/bin/bash

# Bash "strict mode", to help catch problems and bugs in the shell
# script. Every bash script you write should include this. See
# http://redsymbol.net/articles/unofficial-bash-strict-mode/ for
# details.
set -euo pipefail

export DEBIAN_FRONTEND=noninteractive

apt update
apt install -y gnupg2 wget

wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/e6ff8e9c-ee28-47fb-abd7-5c524c983e1c/l_BaseKit_p_2024.2.1.100_offline.sh
sh ./l_BaseKit_p_2024.2.1.100_offline.sh -a --silent --cli --eula accept
rm l_BaseKit_p_2024.2.1.100_offline.sh

wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/d461a695-6481-426f-a22f-b5644cd1fa8b/l_HPCKit_p_2024.2.1.79_offline.sh
sh ./l_HPCKit_p_2024.2.1.79_offline.sh -a --silent --cli --eula accept
rm l_HPCKit_p_2024.2.1.79_offline.sh

# Delete cached files we don't need anymore
apt-get clean
rm -rf /var/lib/apt/lists/*


# link to MKL blas
source /opt/intel/oneapi/setvars.sh
cd /usr/local/lib/R/lib
mv libR.so libR.so.keep

ln -s $MKLROOT/lib/libmkl_rt.so libR.so
export MKL_INTERFACE_LAYER=GNU,LP64
export MKL_THREADING_LAYER=GNU

source $MKLROOT/env/vars.sh
