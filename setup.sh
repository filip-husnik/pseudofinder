#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'
PATH_TO_PSEUDOFINDER=`dirname $0`

printf "\n    ${GREEN}Setting up conda environment...${NC}\n\n"

chmod +x $PATH_TO_PSEUDOFINDER/pseudofinder.py

## creating environment and installing dependencies
conda env create --name pseudofinder --file $PATH_TO_PSEUDOFINDER/modules/environment.yml

## activating environment
source activate pseudofinder

## creating directory for conda-env-specific source files
mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d

## adding codeml-2.ctl file path:
echo '#!/bin/sh'" \

export PATH=\"$(pwd):"'$PATH'\"" \

export ctl=\"$(pwd)/codeml-2.ctl\"" >> ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

# re-activating environment so variable and PATH changes take effect
source activate pseudofinder

printf "\n        ${GREEN}DONE!${NC}\n\n"

# to reset:
# conda env remove --name pseudofinder
