#! /usr/bin/env bash 

### Create SPTCR_ENV
micromamba create -f SPTCR_ENV.yml

## Activate ENV
micromamba activate SPTCR_ENV

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
