#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

printf "\n    ${GREEN}Setting up conda environment...${NC}\n\n"

chmod +x pseudofinder.py

## creating environment and installing dependencies
conda create \
-n pseudofinder \
-c conda-forge \
-c bioconda \
-c defaults \
python=3.10.2 \
biopython=1.79 \
blast=2.12.0 \
diamond=2.0.14 \
pal2nal=14.1 \
paml=4.9 \
plotly=5.6.0 \
pandas=1.4.1 \
numpy=1.22.2 \
reportlab=3.5.68 \
hmmer=3.2.1 \
muscle=3.8.1551 \
prodigal=2.6.3 \
parallel=20190722 \
--yes # Do not ask for confirmation

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
