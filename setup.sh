#! /usr/bin/env bash 

### Install Mamba in Base Environment
conda install -c conda-forge micromamba

### Create SPTCR_ENV
mamba env create -f SPTCR_ENV.yml

## Activate ENV
eval "$(conda shell.bash hook)"
conda activate SPTCR_ENV

## Setup Databases for PyIR to Annotate TCRs
pyir setup

### Download and Build Scripts needed for the Pipeline
cd ./TOOLS

##sc Tagger
git clone https://github.com/vpc-ccg/scTagger.git

## Rattle
git clone --recurse-submodules https://github.com/comprna/RATTLE 
cd RATTLE 
./build.sh
