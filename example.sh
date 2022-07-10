#! /usr/bin/env bash 

### EXAMPLE Pipeline

### PreProcess, corrects, demultiplexes, umicorrects all Fastqs found in a specified Folder. Names The Outfolder like the basname of the Fastq it finds.
### Cleanup is set to False for Debugging
### to apply on your own samples, just specify the path SPTCR_FULL_PIPELINE.sh,the repository, Number of Threads, Memory as well as the path to your sample folder

THREADS=12
MEMORY=16

full="/PATH/TO/GITHUB/REPOSITORY/SPTCR-seq Pipeline Scripts/SPTCR_FULL_Pipeline.sh"

for sample in "/YOUR/SAMPLE/FOLDER/WITH/*.fastq"
do  
    SAMPLE_NAME="$(basename "${sample}")"
    SAMPLE_NAME="$(cut -d'_' -f1 <<<${SAMPLE_NAME})"

    OUT="/YOUR/DESIRED/OUTPAH/${SAMPLE_NAME}"
    echo "$SAMPLE_NAME"
    mkdir "${OUT}"
    cd "${OUT}"

    ############################
    "${full}" \
    -n "$SAMPLE_NAME" \
    -i"${sample}" \
    -t ${THREADS} \
    -mem ${MEMORY} \
    -rep "/PATH/TO/GITHUB/REPOSITORY/SPTCR-Seq-Pipeline" \
    -cln False 
    
    
done
