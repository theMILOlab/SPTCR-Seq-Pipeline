#!/usr/bin/env bash
##### Argparse Scripts #####

# Use python's argparse module in shell scripts
#
# The function `argparse` parses its arguments using
# argparse.ArgumentParser; the parser is defined in the function's
# stdin.
#
# Executing ``argparse.bash`` (as opposed to sourcing it) prints a
# script template.
#
# https://github.com/nhoffman/argparse-bash
# MIT License - Copyright (c) 2015 Noah Hoffman

argparse(){
    argparser=$(mktemp 2>/dev/null || mktemp -t argparser)
    cat > "$argparser" <<EOF
from __future__ import print_function
import sys 
import argparse
import os


class MyArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        """Print help and exit with error"""
        super(MyArgumentParser, self).print_help(file=file)
        sys.exit(1)

parser = MyArgumentParser(prog=os.path.basename("$0"),
            description="""$ARGPARSE_DESCRIPTION""")
EOF

    # stdin to this function should contain the parser definition
    cat >> "$argparser"

    cat >> "$argparser" <<EOF
args = parser.parse_args()
for arg in [a for a in dir(args) if not a.startswith('_')]:
    key = arg.upper()
    value = getattr(args, arg, None)

    if isinstance(value, bool) or value is None:
        print('{0}="{1}";'.format(key, 'yes' if value else ''))
    elif isinstance(value, list):
        print('{0}=({1});'.format(key, ' '.join('"{0}"'.format(s) for s in value)))
    else:
        print('{0}="{1}";'.format(key, value))
EOF

    # Define variables corresponding to the options if the args can be
    # parsed without errors; otherwise, print the text of the error
    # message.
    if python "$argparser" "$@" &> /dev/null; then
        eval $(python "$argparser" "$@")
        retval=0
    else
        python "$argparser" "$@"
        retval=1
    fi

    rm "$argparser"
    return $retval
}
##########################################
##### Argparse Options #####
ARGPARSE_DESCRIPTION="Pipeline to group & Correct T-Cell Receptor Reads generated by Oxford Nanopore Reads for Libraries prepared for 10X Genomics"
argparse "$@" <<EOF || exit 1

parser.add_argument('-n','--NAME',help="Name of Output Folder",default="-")
parser.add_argument('-i', '--INPUT_FASTQ', help="Specify the Path (preprocessed) Fastq File",required=True)
parser.add_argument('-b', '--INPUT_IGB', help="Specify the Path to the Input IgBlast Result",required=True)
parser.add_argument('-o', '--OUTFOLDER', help="Specify the Directory for the Outputfolder", default="PWD")

parser.add_argument('-t','--THREADS',help="Number of Threads", default="2")

parser.add_argument('-rep', '--REPOSITORY', help="Specify the Location of the Repositroy Folder holding all References and scripts for SPTCR Seq",default="./")

parser.add_argument('-igb', '--IGBLAST', help="If True, the corrected Fastq is aligned with IgBLAST Following Correction.",default="True")

parser.add_argument('-cln', '--CLEANUP', help="If True, created intermediate Files and Folders will be deleted. For Debugging you can set this to False.",default="True")

parser.add_argument('-bc', '--BCUMI', help="Specify the path to the SAMPLENAME_barcode_umi.csv generated by the preprocessing pipeline. Defaults to OUTFOLDER/PreProcessing/SAMPLENAME_barcode_umi.csv",default="True")

EOF

################################################################
################## Define Variables BLOCK###########################
################################################################
if [ "${OUTFOLDER}" = "PWD" ];then
    OUTFOLDER=$PWD
    mkdir ${PWD}/ClusterCorrect
    out="${OUTFOLDER}"
    OUTFOLDER=${PWD}/ClusterCorrect
else
    mkdir "${OUTFOLDER}"/ClusterCorrect
    out="${OUTFOLDER}"
    OUTFOLDER="${OUTFOLDER}"/ClusterCorrect
fi 

if [ "${NAME}" = "-" ];then
    SAMPLE_NAME="$(basename "${INPUT_FASTQ}")"
    SAMPLE_NAME="$(cut -d'.' -f1 <<<"${SAMPLE_NAME}")"_$(date +%d_%Y)
else
    SAMPLE_NAME="${NAME}"
fi 


### Log Folder ####
mkdir "${OUTFOLDER}"/LOGS
LOGS="${OUTFOLDER}"/LOGS

### Timestamp Function

timestamp() {
    date +"%Y-%m-%d : %H-%M-%S" # current time
    }

STARTTIME=$(date +%s)


### Rattle Path
RATTLE_PATH="${REPOSITORY}"/TOOLS/RATTLE

### For TCR Annotation Summay
DEMUX_SUMMARY="${REPOSITORY}"/SCRIPTS/demultiplex_summarize.py
BARCODE_UMI_FILE="${out}"/PreProcessing/"${SAMPLE_NAME}"_barcode_umi.csv

################################################################
################## Cluster BLOCK #######################
################################################################

echo " :::: Clustering Reads by VJ Arrangement ::::"

mkdir ""${OUTFOLDER}""/IGB_CLUSTERS

python "${REPOSITORY}"/SCRIPTS/TCR_GROUPING_IGB.py \
    --IGB "${INPUT_IGB}" \
    --INPUT_FASTQ "${INPUT_FASTQ}" \
    --OUT "${OUTFOLDER}"/IGB_CLUSTERS 2> "${LOGS}"/tcr_clustering_log.txt

IGB_CLUSTERS="${OUTFOLDER}"/IGB_CLUSTERS



################################################################
################## CORRECTION BLOCK #######################
################################################################

echo " :::: Parallel Correction of all Clusters ::::"

mkdir "${OUTFOLDER}"/CORRECTION
mkdir "${OUTFOLDER}"/CORRECTION/RATTLE_CLUSTERS
RATTLE_CLUSTERS_OUT="${OUTFOLDER}"/CORRECTION/RATTLE_CLUSTERS

mkdir "${OUTFOLDER}"/CORRECTION/RATTLE_CORRECT
RATTLE_CORRECT_OUT="${OUTFOLDER}"/CORRECTION/RATTLE_CORRECT

mkdir "${OUTFOLDER}"/CORRECTION/CORRECTED_MERGE
CORRECTED_MERGE="${OUTFOLDER}"/CORRECTION/CORRECTED_MERGE

mkdir "${LOGS}"/Rattle_Clustering 
mkdir "${LOGS}"/Rattle_Correction
#--fastq \

find "${IGB_CLUSTERS}" -type f |parallel --jobs 0 \
        "
        mkdir '${RATTLE_CLUSTERS_OUT}'/{/.}
        '${RATTLE_PATH}'/rattle cluster \
            -i {} \
            -t ${THREADS} \
            -o '${RATTLE_CLUSTERS_OUT}'/{/.} \
            --iso \
            2> '${LOGS}'/Rattle_Clustering/{/.}_clust_stderr.txt

        mkdir '${RATTLE_CORRECT_OUT}'/{/.}
        '${RATTLE_PATH}'/rattle correct \
            -i {} \
            -c '${RATTLE_CLUSTERS_OUT}'/{/.}/clusters.out \
            -t ${THREADS} \
            -o '${RATTLE_CORRECT_OUT}'/{/.} \
            2> '${LOGS}'/Rattle_Correction/{/.}_corr_stderr.txt
        "
echo " :::: Parallel Merge Corrected Fastqs ::::"
find "${RATTLE_CORRECT_OUT}" -type f -name "corrected.fq"|parallel --jobs 0 "cat {} >> "${CORRECTED_MERGE}"/"${SAMPLE_NAME}"_corrected_merged.fastq"

echo " :::: Repairing possibly broken Merge Corrected Fastqs ::::"

seqkit sana "${CORRECTED_MERGE}"/"${SAMPLE_NAME}"_corrected_merged.fastq -o "${CORRECTED_MERGE}"/"${SAMPLE_NAME}"_corrected_merged_sana.fastq 2> "${LOGS}"/seqkit_sana.txt

rm "${CORRECTED_MERGE}"/"${SAMPLE_NAME}"_corrected_merged.fastq
mv "${CORRECTED_MERGE}"/"${SAMPLE_NAME}"_corrected_merged_sana.fastq "${CORRECTED_MERGE}"/"${SAMPLE_NAME}"_corrected_merged.fastq

CORRECTED_MERGED_FASTQ="${CORRECTED_MERGE}"/"${SAMPLE_NAME}"_corrected_merged.fastq

################################################################
################## Final IGB Query BLOCK #######################
################################################################

if [ ${IGBLAST} = True ]; then
    echo " :::: Quering corrected Fastq to IGB :::: "
    mkdir "${OUTFOLDER}"/IGB_Corrected
    mkdir "${OUTFOLDER}"/IGB_Corrected/TEMP_"${SAMPLE_NAME}"
    cd "${OUTFOLDER}"/IGB_Corrected

    pyir \
        -t fastq \
        --outfmt tsv \
        -m ${THREADS} \
        --pretty \
        --debug \
        -r TCR \
        -s human \
        --numV 1 \
        --numD 1 \
        --numJ 1 \
        --tmp_dir "${OUTFOLDER}"/IGB_Corrected/TEMP_"${SAMPLE_NAME}" \
        -o "${SAMPLE_NAME}"_corrected_IGB \
        ${CORRECTED_MERGED_FASTQ} 2> "${LOGS}"/IGB_Annotation.txt

    gunzip "${SAMPLE_NAME}"_corrected_IGB.tsv.gz

    OUTPUT_IGB="${OUTFOLDER}"/IGB_Corrected/"${SAMPLE_NAME}"_corrected_IGB.tsv

else 
    echo " :::: Not quering IgBlast as indicated ::::"
fi


if [ ${CLEANUP} = True ]; then
    echo " :::: Removing Intermediate Files ::::"
    ### Move Output Files to the Front
    mv "${CORRECTED_MERGED_FASTQ}" "${OUTFOLDER}"/"${SAMPLE_NAME}"_corrected_merged.fastq
    mv "${OUTFOLDER}"/IGB_Corrected/"${SAMPLE_NAME}"_corrected_IGB.tsv "${OUTFOLDER}"/"${SAMPLE_NAME}"_corrected_IGB.tsv
    IGBLAST="${OUTFOLDER}"/"${SAMPLE_NAME}"_corrected_IGB.tsv

    ### Remove Created Working Directories
    rm -r "${OUTFOLDER}"/IGB_Corrected/TEMP_"${SAMPLE_NAME}"
    rm -r "${OUTFOLDER}"/CORRECTION
    
    rm -r "${IGB_CLUSTERS}"
    rm -r "${OUTFOLDER}"/CORRECTION/RATTLE_CORRECT
    rm -r "${OUTFOLDER}"/CORRECTION/RATTLE_CLUSTERS

else echo " :::: Keeping all Intermediate Files Intermediate Files ::::"
fi

echo " :::: Generate the VDj Annotation Summary ::::"

DEMUX_SUMMARY="/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Github SPTCR-seq Pipeline/SPTCR-Seq-Pipeline/SCRIPTS/demultiplex_summarize.py"
OUTFOLDER="/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/LONG_TCR/Testing_Zone/28.7/ClusterCorrect"
IGBLAST="/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/LONG_TCR/Testing_Zone/28.7/ClusterCorrect/test_demux_corrected_IGB.tsv"
BARCODE_UMI_FILE="/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/LONG_TCR/Testing_Zone/27.07/Demultiplexing_test_demux/test_demux_barcode_umi.csv"

python "${DEMUX_SUMMARY}" \
    -igb "${IGBLAST}" \
    -bc  "${BARCODE_UMI_FILE}" \
    --MOD False \
    --OUTN _corr_igb \
    -n "${SAMPLE_NAME}" \
    -o "${OUTFOLDER}" 
    
echo "Done with SPTCR-seq Correction Pipeline. Corrected IgBlast File is in ${OUTDIR}/"${SAMPLE_NAME}"_corrected_IGB.tsv"

exit

### Missing: Demultiplex Final IGB & UMI Correct

############### For umi Correction
mkdir "${OUTFOLDER}"/UMI_Correction
mkdir "${OUTFOLDER}"/UMI_Correction/UMI_CLUSTERING
mkdir "${OUTFOLDER}"/UMI_Correction/UMI_CORRECTED

/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/LONG_TCR/FINAL_21_2/Repository/SCRIPTS/umi_correction.py \
    -igb "${OUTPUT_IGB} "\
    -n "${SAMPLE_NAME}" \
    -O "${OUTFOLDER}"/UMI_Correction
