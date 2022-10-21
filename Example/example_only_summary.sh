#! /usr/bin/env bash 

## Example Pipeline to Demultiplex & Summarize given IGB Table

REPOSITORY="./SPTCR-Seq-Pipeline"
OUT="./SPTCR-Seq-Pipeline/Example/Only_Summary"
IN="./SPTCR-Seq-Pipeline/Example/Example_TCR.fastq_sana" 
INGB="./SPTCR-Seq-Pipeline/Example/Example/ClusterCorrect/Example_corrected_IGB.tsv"

THREADS=12
MEMORY=24


SAMPLE_NAME="$(basename "${IN}")"
SAMPLE_NAME="$(cut -d'_' -f1 <<<${SAMPLE_NAME})"


bash "./SPTCR-Seq-Pipeline/SPTCR-seq Pipeline Scripts/Demultiplex_IGB.sh" \
    -i "${IN}" \
    -n "$SAMPLE_NAME" \
    -igb "${INGB}" \
    -o "${OUT}" \
    -t  ${THREADS} \
    -mem ${MEMORY} \
    -rep "${REPOSITORY}"

