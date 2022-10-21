#! /usr/bin/env bash 

### EXAMPLE Pipeline

### PreProcess, corrects, demultiplexes, umicorrects all Fastqs found in a specified Folder. Names The Outfolder like the basname of the Fastq it finds.
### Cleanup is set to False for Debugging
### to apply on your own samples, just specify the path SPTCR_FULL_PIPELINE.sh,the repository, Number of Threads, Memory as well as the path to your sample folder

full="../SPTCR-seq Pipeline Scripts/SPTCR_FULL_Pipeline.sh"

sample="./Example_TCR.fastq"
THREADS=12
MEMORY=32
SAMPLE_NAME="$(basename "${sample}")"
SAMPLE_NAME="$(cut -d'_' -f1 <<<${SAMPLE_NAME})"

OUT="../Example/${SAMPLE_NAME}"
echo "$SAMPLE_NAME"
cd "${OUT}"

############################
"${full}" \
    -n "$SAMPLE_NAME" \
    -i"${sample}" \
    -t ${THREADS} \
    -mem ${MEMORY} \
    -rep ".." \
    -cln False \
    -o "${OUT}"
