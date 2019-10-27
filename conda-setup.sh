#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

printf "\n    ${GREEN}Setting up conda environment...${NC}\n\n"

chmod +x pseudofinder.py

## adding conda channels
conda config --add channels defaults 2> /dev/null
conda config --add channels bioconda 2> /dev/null
conda config --add channels conda-forge 2> /dev/null
conda config --add channels au-eoed 2> /dev/null

## creating GToTree environment and installing dependencies
conda create -n pseudo blast diamond pal2nal muscle paml prodigal biopython plotly pandas numpy reportlab --yes

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
