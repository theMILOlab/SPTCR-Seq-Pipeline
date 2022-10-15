#! /usr/bin/env bash 

### EXAMPLE Pipeline

### PreProcess, corrects, demultiplexes, umicorrects all Fastqs found in a specified Folder. Names The Outfolder like the basname of the Fastq it finds.
### Cleanup is set to False for Debugging
### to apply on your own samples, just specify the path SPTCR_FULL_PIPELINE.sh,the repository, Number of Threads, Memory as well as the path to your sample folder

THREADS=12
MEMORY=64

#full="/PATH/TO/GITHUB/REPOSITORY/SPTCR-seq Pipeline Scripts/SPTCR_FULL_Pipeline.sh"
full="/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Github SPTCR-seq Pipeline/SPTCR-Seq-Pipeline/SPTCR-seq Pipeline Scripts/SPTCR_FULL_Pipeline.sh"
for sample in "/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Github SPTCR-seq Pipeline/SPTCR-Seq-Pipeline/Example/"*.fastq
do  
    SAMPLE_NAME="$(basename "${sample}")"
    SAMPLE_NAME="$(cut -d'_' -f1 <<<${SAMPLE_NAME})"

    OUT="/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Github SPTCR-seq Pipeline/SPTCR-Seq-Pipeline/Example/"
    "${SAMPLE_NAME}"
    #echo "$SAMPLE_NAME"
    #mkdir "${OUT}"
    cd "${OUT}"

    ############################
    "${full}" \
        -n "$SAMPLE_NAME" \
        -i"${sample}" \
        -t ${THREADS} \
        -mem ${MEMORY} \
        -rep "/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Github SPTCR-seq Pipeline/SPTCR-Seq-Pipeline" \
        -cln False \
        -o "${OUT}"
    #-rep "/PATH/TO/GITHUB/REPOSITORY" \
    
done
