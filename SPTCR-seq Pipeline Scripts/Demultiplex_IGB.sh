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
ARGPARSE_DESCRIPTION="Pipeline to demultiplex & UMI Correct given IGB Table with raw Fastq & Barcodes"
argparse "$@" <<EOF || exit 1

parser.add_argument('-n','--NAME',help="Chosen Name for the Pipeline. Will Be Used as Name for the Pipelines Outfolder.",default="-")
parser.add_argument('-i', '--INPUT_FASTQ', help="Specify the Path to the raw Input Fastq File",required=True)
parser.add_argument('-igb', '--IGBLAST', help="Provide IGB for demultiplexing & UMI Correction",required=True)
parser.add_argument('-o', '--OUTFOLDER', help="Specify the Directory for the Output Folder. If not specified, will use current working Directory.", default="PWD")

parser.add_argument('-t','--THREADS',help="Number of Threads to use.", default="2")
parser.add_argument('-mem','--MEMORY',help="Gigabytes of RAM to use for Demultiplexing", default="8")

parser.add_argument('-rep', '--REPOSITORY', help="Specify the Location of the Github Repository Folder for SPTCR Seq",default="./")

parser.add_argument('-bc', '--BARCODES', default='visium_bc.tsv', 
                    help='Specify the Path to the Barcode Whitelist extracted from tissue_positions_list.csv from the Spaceranger Output. If left unfilled, all possible Visium Barcodes will be matched.')

EOF

################################################################
################## Define Variables BLOCK###########################
################################################################
if [ "${OUTFOLDER}" = "PWD" ];then
    mkdir "${PWD}"/PreProcessing
    OUTFOLDER="${PWD}"/PreProcessing
else
    OUTFOLDER="${OUTFOLDER}"
fi 

if [ ${NAME} = "-" ];then
    SAMPLE_NAME="$(basename ${INPUT_FASTQ})"
    SAMPLE_NAME="$(cut -d'.' -f1 <<<"${SAMPLE_NAME}")"_$(date +%d_%Y)
else
    SAMPLE_NAME=${NAME}
fi 


### Log Folder ####
mkdir "${OUTFOLDER}"/LOGS
LOGS="${OUTFOLDER}"/LOGS

### Timestamp Function

timestamp() {
    date +"%Y-%m-%d : %H-%M-%S" # current time
    }

STARTTIME=$(date +%s)

### For Barcode Demultiplexing

if [ "${BARCODES}" = "visium_bc.tsv" ];then
    BARCODES="${REPOSITORY}"/Reference/Barcodes/visium_bc.tsv
    
else
    BARCODES="${BARCODES}"
fi 

##Demultiplex
Demultiplex_UMI_Extraction="${REPOSITORY}/SPTCR-seq Pipeline Scripts/1_Demultiplex_UMI_Extraction.sh"

################################################################
################## Demultiplexing #######################
################################################################

echo " ::::: Demultiplexing Input IGB with given Barcodes & Fastq :::::"
cd "${OUTFOLDER}"

"${Demultiplex_UMI_Extraction}" \
    -i "${INPUT_FASTQ}" \
    -igb "${IGBLAST}" \
    -n "${SAMPLE_NAME}" \
    -t ${THREADS} \
    -mem ${MEMORY} \
    -bc "${BARCODES}" \
    -rep "${REPOSITORY}"


echo " :::: Performing UMI Correction on Uncorrected Summary ::::"#

"${REPOSITORY}/SCRIPTS/umi_correct_output.py" \
    -igb "${OUTFOLDER}/Demultiplexing_${SAMPLE_NAME}/${SAMPLE_NAME}_vdj_umi_barcode_uncorrected_df.csv" \
    -n "${SAMPLE_NAME}" \
    -outn 'UNCORRECTED_umi_corrected_count_table' \
    --BARCODES "${BARCODES}" \
    -O "${OUTFOLDER}" 


