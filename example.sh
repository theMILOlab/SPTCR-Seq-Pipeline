#! /usr/bin/env bash 

### EXAMPLE Pipeline

### PreProcess, corrects, demultiplexes, umicorrects all Fastqs found in a specified Folder. Names The Outfolder like the basname of the Fastq it finds.
### Cleanup is set to False for Debugging
### to apply on your own samples, just specify the path SPTCR_FULL_PIPELINE.sh,the repository, Number of Threads, Memory as well as the path to your sample folder

THREADS=12
MEMORY=64


full="/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Github SPTCR-seq Pipeline/SPTCR-Seq-Pipeline/SPTCR-seq Pipeline Scripts/SPTCR_FULL_Pipeline.sh"
sample="/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Github SPTCR-seq Pipeline/SPTCR-Seq-Pipeline/Example/Example_TCR.fastq"

SAMPLE_NAME="$(basename "${sample}")"
SAMPLE_NAME="$(cut -d'_' -f1 <<<${SAMPLE_NAME})"

OUT="/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Github SPTCR-seq Pipeline/SPTCR-Seq-Pipeline/${SAMPLE_NAME}"
echo "$SAMPLE_NAME"
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

exit


full="PATH/TO/GITHUB_REPO/SPTCR-Seq-Pipeline/SPTCR-seq Pipeline Scripts/SPTCR_FULL_Pipeline.sh"
sample="PATH/TO/GITHUB_REPO/SPTCR-Seq-Pipeline/Example/Example_TCR.fastq"

SAMPLE_NAME="$(basename "${sample}")"
SAMPLE_NAME="$(cut -d'_' -f1 <<<${SAMPLE_NAME})"

OUT="PATH/TO/GITHUB_REPO/SPTCR-Seq-Pipeline/Example/${SAMPLE_NAME}"
echo "$SAMPLE_NAME"
cd "${OUT}"

############################

"${full}" \
    -n "$SAMPLE_NAME" \
    -i"${sample}" \
    -t ${THREADS} \
    -mem ${MEMORY} \
    -rep "PATH/TO/GITHUB_REPO/SPTCR-Seq-Pipeline" \
    -cln False \
    -o "${OUT}"