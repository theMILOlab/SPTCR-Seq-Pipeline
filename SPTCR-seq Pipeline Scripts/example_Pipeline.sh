#! /usr/bin/env bash

demux=./1_Demultiplex_UMI_Extraction.sh
preproc=./2_Preprocess_Reads.sh
clusco=./3_Cluster_Correct.sh

##############################

${demux}     \
    -i ./test_demux.fastq 

${preproc} \
    -i ./test_demux.fastq 

${clusco} \
    -i ./test_demux.fastq \
    -b ./PreProcessing/IGB_Trimmed/test_demux_preprocessed_IGB.tsv