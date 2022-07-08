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
ARGPARSE_DESCRIPTION="Pipeline to barcode and extract umir regions of Nanopore Reads for Libraries prepared for 10X Genomics"
argparse "$@" <<EOF || exit 1

parser.add_argument('-n','--NAME',help="Name of Output Folder",default="-")
parser.add_argument('-i', '--INPUT_FASTQ', help="Specify the Path to the merged Input Fastq File",required=True)
parser.add_argument('-igb', '--INPUT_IGB', help="Specify the Path to IgBlast File to be demultiplexed with demultiplex_summarize.py. If not specified will only output table for later demultiplexing.",required=False)

parser.add_argument('-o', '--OUTFOLDER', help="Specify the Directory for the Outputfolder", default="PWD")

parser.add_argument('-t','--THREADS',help="Number of Threads", default="2")
parser.add_argument('-mem','--MEMORY',help="RAM to user", default="16")

parser.add_argument('-rep', '--REPOSITORY', help="Specify the Location of the Repositroy Folder holding all References and scripts for SPTCR Seq",default="../")
parser.add_argument('-a', '--ADAPTER', default='CTACACGACGCTCTTCCGATCT', 
                    help='Specify Read 1 Sequence')

EOF

################################################################
###################### Variable BLOCK###########################
################################################################

#### Outdir ####
if [ "${OUTFOLDER}" = PWD ];then
    #mkdir ${PWD}/Demultiplexing
    OUTFOLDER=${PWD}
else
    OUTFOLDER="${OUTFOLDER}"
fi 

if [ ${NAME} = "-" ];then
    SAMPLE_NAME="$(basename "${INPUT_FASTQ}")"
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


BARCODE_LIST="${REPOSITORY}"/Reference/Barcodes/visium_bc.tsv
DEMUX_GENERATOR="${REPOSITORY}"/SCRIPTS/generate_demux_umi_dfs.sh
DEMUXXER="${REPOSITORY}"/SCRIPTS/demultiplex_extract_umi_region.py
SCTAGGER="${REPOSITORY}"/TOOLS/scTagger/py
DEMUX_SUMMARY="${REPOSITORY}"/SCRIPTS/demultiplex_summarize.py

################################################################
################## UMI Extraction BLOCK###########################
################################################################

echo "$(timestamp)"
echo " :::: Demultiplexing & Extracting UMIs with scTagger::::"
mkdir "${OUTFOLDER}"/Demultiplexing_"${SAMPLE_NAME}"
OUTFOLDER="${OUTFOLDER}"/Demultiplexing_"${SAMPLE_NAME}"

echo ":::::Demultiplexing and UMI Extraction ${samplename}:::::"
mkdir "${OUTFOLDER}"/"${SAMPLE_NAME}"_UMI_Extraction

echo "::: Extracting Adapter and Barcode :::"
"${SCTAGGER}"/extract_lr_bc.py \
    -r "${INPUT_FASTQ}" \
    -t ${THREADS} \
    -sa ${ADAPTER} \
    --num-bp-after 16 \
    -o "${OUTFOLDER}"/"${SAMPLE_NAME}"_UMI_Extraction/"${SAMPLE_NAME}"_Extracted_Adapter_Barcode.tsv.gz 2>> "${LOGS}"/"${SAMPLE_NAME}"_Demux_stderr.txt


echo "::: Extracting Adapter, Barcode and UMI:::"
"${SCTAGGER}"/extract_lr_bc.py \
    -r "${INPUT_FASTQ}" \
    -t ${THREADS} \
    -sa ${ADAPTER} \
    --num-bp-after 28 \
    -o "${OUTFOLDER}"/"${SAMPLE_NAME}"_UMI_Extraction/"${SAMPLE_NAME}"_Extracted_Adapter_Barcode_UMI.tsv.gz 2>> "${LOGS}"/"${SAMPLE_NAME}"_Demux_stderr.txt

################################################################
##################Demultiplexing BLOCK###########################
################################################################
echo "::::: Extracting Long Read Adapter Segments :::::"
mkdir "${OUTFOLDER}"/"${SAMPLE_NAME}"_DEMUX

"${SCTAGGER}"/extract_lr_bc.py \
    -r "${INPUT_FASTQ}" \
    -t ${THREADS} \
    -sa ${ADAPTER} \
    -o "${OUTFOLDER}"/"${SAMPLE_NAME}"_DEMUX/"${SAMPLE_NAME}"_lr_bc.tsv.gz 2>> "${LOGS}"/"${SAMPLE_NAME}"_Demux_stderr.txt

echo "::::: Matching Barcodes to Extracted Segments :::::"
"${SCTAGGER}"/match_lr_bc-trie.py \
    -lr "${OUTFOLDER}"/"${SAMPLE_NAME}"_DEMUX/"${SAMPLE_NAME}"_lr_bc.tsv.gz \
    -sr "${BARCODE_LIST}"  \
    -o "${OUTFOLDER}"/"${SAMPLE_NAME}"_DEMUX/"${SAMPLE_NAME}"_scTagger_DEMUX.tsv.gz \
    -t ${THREADS} \
    -m ${MEMORY} 2>> "${LOGS}"/"${SAMPLE_NAME}"_Demux_stderr.txt



gunzip "${OUTFOLDER}"/"${SAMPLE_NAME}"_DEMUX/*
#"${SAMPLE_NAME}"_scTagger_DEMUX.tsv.gz

gunzip "${OUTFOLDER}"/"${SAMPLE_NAME}"_UMI_Extraction/*

DEMUXED="${OUTFOLDER}"/"${SAMPLE_NAME}"_DEMUX/"${SAMPLE_NAME}"_scTagger_DEMUX.tsv 
umi_folder="${OUTFOLDER}"/"${SAMPLE_NAME}"_UMI_Extraction

################################################################
##########Generate Demultiplexing Table, Extract UMIs#########
################################################################
echo "::::: Generating Final Output Table :::::"

python "${DEMUXXER}" \
    -i ${DEMUXED} \
    -u ${umi_folder} \
    -o "${OUTFOLDER}" \
    -n "${SAMPLE_NAME}" \
    -b "${BARCODE_LIST}" >> "${LOGS}"/"${SAMPLE_NAME}"_Demux_stderr.txt 2>> "${LOGS}"/"${SAMPLE_NAME}"_Demux_stderr.txt

echo ":::::: Removing Intermediate Folders :::::"
rm -r ${umi_folder}
rm -r "${OUTFOLDER}"/"${SAMPLE_NAME}"_DEMUX

if [ -f "${INPUT_IGB}" ]; then
    echo ":::::: Generating IGB with UMI :::::"

    python "${DEMUX_SUMMARY}" \
        -igb "${INPUT_IGB}" \
        -bc "${OUTFOLDER}"/"${SAMPLE_NAME}"_barcode_umi.csv \
        -o "${OUTFOLDER}"
else
    echo " Not Input IgBlast File given"
fi
