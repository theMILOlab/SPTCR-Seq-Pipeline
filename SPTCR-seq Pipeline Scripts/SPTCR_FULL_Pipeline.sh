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
ARGPARSE_DESCRIPTION="Full Pipeline to Demultiplex, PreProcess & Correct T-Cell Receptor Reads optained by Oxford Nanopore Long Read Sequencing"
argparse "$@" <<EOF || exit 1

parser.add_argument('-n','--NAME',help="Sample Name/Name of Output Folder",default="-")
parser.add_argument('-i', '--INPUT_FASTQ', help="Specify the Path (preprocessed) Fastq File",required=True)
parser.add_argument('-o', '--OUTFOLDER', help="Specify the Directory for the Outputfolder", default=".")

parser.add_argument('-t','--THREADS',help="Number of Threads", default="2")
parser.add_argument('-mem','--MEMORY',help="RAM to user", default="8")

parser.add_argument('-rep', '--REPOSITORY', help="Specify the Location of the Github Repository Folder for SPTCR Seq",default="./Github SPTCR-seq Pipeline/SPTCR-Seq-Pipeline")
parser.add_argument('-cln', '--CLEANUP', help="If True, created intermediate Files and Folders will be deleted. For Debugging you can set this to False.",default="True")

EOF
################################################################
################## Define Variables BLOCK#######################
################################################################

preproc="${REPOSITORY}/SPTCR-seq Pipeline Scripts/2_Preprocess_Reads.sh"
cluscorr="${REPOSITORY}/SPTCR-seq Pipeline Scripts/3_Cluster_Correct.sh"

OUT="${OUTFOLDER}/${NAME}"
mkdir "${OUT}"

################################################################
################## Scrip Execution BLOCK########################
################################################################

echo "################## PreProcessing $NAME ########################"
cd "${OUT}"
bash "${preproc}" \
        -n ${NAME} \
        -i "${INPUT_FASTQ}" \
        -t ${THREADS} \
        -mem ${MEMORY} \
        -rep "${REPOSITORY}"

echo "################## Correcting $NAME ########################"
bash "${cluscorr}" \
        -i "./PreProcessing/${NAME}_Cutadapt_trimmed_sana.fastq" \
        -b "./PreProcessing/${NAME}_preprocessed_IGB.tsv" \
        -n ${NAME} \
        -t ${THREADS} \
        --GROUPER "locus,v_family" \
        -rep "${REPOSITORY}" \
        -cln ${CLEANUP}

        #-b "./PreProcessing/Demultiplexing_${NAME}/${NAME}_vdj_umi_barcode_uncorrected_df.csv" \