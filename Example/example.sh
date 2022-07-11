#! /usr/bin/env bash 

full="../SPTCR-seq Pipeline Scripts/SPTCR_FULL_Pipeline.sh"

sample="./Example_TCR.fastq"

SAMPLE_NAME="$(basename "${sample}")"
SAMPLE_NAME="$(cut -d'_' -f1 <<<${SAMPLE_NAME})"

OUT="../Example/${SAMPLE_NAME}"
echo "$SAMPLE_NAME"

############################
"${full}" \
    -n "$SAMPLE_NAME" \
    -i"${sample}" \
    -t 8 \
    -mem 16 \
    -rep ".." \
    -cln False 
