#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

printf "\n    ${GREEN}Setting up conda environment...${NC}\n\n"

chmod +x pseudofinder.py

## adding conda channels (THIS PART DEPRACATED)
# conda config --add channels defaults 2> /dev/null
# conda config --add channels bioconda 2> /dev/null
# conda config --add channels conda-forge 2> /dev/null
# conda config --add channels au-eoed 2> /dev/null

## creating GToTree environment and installing dependencies
conda create -n gtotree -c conda-forge -c bioconda -c defaults biopython hmmer=3.2.1 muscle=3.8.1551 trimal=1.4.1 fasttree=2.1.10 iqtree=1.6.9 prodigal=2.6.3 taxonkit=0.3.0 parallel=20190722 --yes

## activating environment
source activate pseudo

## creating directory for conda-env-specific source files
mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d

## adding codeml-2.ctl file path:
echo '#!/bin/sh'" \
export PATH=\"$(pwd):"'$PATH'\"" \
export ctl=\"$(pwd)/codeml-2.ctl\"" >> ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

# re-activating environment so variable and PATH changes take effect
source activate pseudo


printf "\n        ${GREEN}DONE!${NC}\n\n"
