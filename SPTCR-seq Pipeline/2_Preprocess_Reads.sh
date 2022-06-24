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
parser.add_argument('-o', '--OUTFOLDER', help="Specify the Directory for the Outputfolder", default="PWD")

parser.add_argument('-t','--THREADS',help="Number of Threads", default="2")
parser.add_argument('-mem','--MEMORY',help="RAM to user", default="8")

parser.add_argument('-rep', '--REPOSITORY', help="Specify the Location of the Repositroy Folder holding all References and scripts for SPTCR Seq",default="/mnt/681ABFBA1ABF839A/Dropbox/KBJasim/Projects/Capture_Sequencing/LONG_TCR/FINAL/Repository")

parser.add_argument('-pri', '--PRIMER', help="Specify Custom Primers if not having used either 10X Visium or Single Cell for the reconstruction of Full Reads by Pychopper.",default="10X")
parser.add_argument('-conf', '--CONFIGURATION', help="Specify the possible Configurations of the Custom Primers for Pychopper if not having used either 10X Visium or Single Cell for the reconstruction of Full Reads by Pychopper.")

parser.add_argument('-chop', '--PYCHOPPER', help="Specify if Reads should be made full length by Pychopper",default="True")
parser.add_argument('-trim', '--ADAPTER_TRIM', help="Specify if Reads should be from Adapters",default="True")
parser.add_argument('-igb', '--IGBLAST', help="If True, the preprocessed Fastq is aligned with IgBLAST Following Processing.",default="True")

EOF

################################################################
################## Define Variables BLOCK###########################
################################################################
if [ ${OUTFOLDER} = "PWD" ];then
    OUTFOLDER=$PWD
else
    OUTFOLDER=${OUTFOLDER}
fi 

if [ ${NAME} = "-" ];then
    SAMPLE_NAME="$(basename ${INPUT_FASTQ})"
    SAMPLE_NAME="$(cut -d'.' -f1 <<<${SAMPLE_NAME})"_$(date +%d_%Y)
else
    SAMPLE_NAME=${NAME}
fi 


### Log Folder ####
mkdir ${OUTFOLDER}/LOGS
LOGS=${OUTFOLDER}/LOGS

### Timestamp Function

timestamp() {
    date +"%Y-%m-%d : %H-%M-%S" # current time
    }

STARTTIME=$(date +%s)


### PYCHOPPER
if [ ${PRIMER} = 10X ];then
    PRIMERS=${REPOSITORY}/REFERENCE/Primer/Pychopper/10XPrimers_pychopper.fa
    PRIMER_CONFIGURATION=${REPOSITORY}/REFERENCE/Primer/Pychopper/10XPrimers_pychopper_configuration.txt
else
    PRIMER=${PRIMER}
    PRIMER_CONFIGURATION=${CONFIGURATION}
fi

##Cutadapt
DUAL_ADAPTER_10X=${REPOSITORY}/REFERENCE/Primer/Cutadapt/10X_Dual_Adapter.fa
Adapter_5_3_10X=${REPOSITORY}/REFERENCE/Primer/Cutadapt/5_3_10X_Adapter.fa

################################################################
################## Trim & Reorient BLOCK #######################
################################################################



if [ ${PYCHOPPER} = True ]; then
    echo "$(timestamp)" 
    echo " :::: Pychopper ::::"
    mkdir ${OUTFOLDER}/PYCHOPPER
    pychop_dir=${OUTFOLDER}/PYCHOPPER

    cdna_classifier.py \
        -m edlib \
        -b ${PRIMERS} \
        -c ${PRIMER_CONFIGURATION} \
        -r ${LOGS}/${SAMPLE_NAME}_Pychopper_report.pdf \
        -S ${LOGS}/${SAMPLE_NAME}_Pychopper_report.tsv \
        -t ${THREADS} \
        -p \
        ${INPUT_FASTQ} \
        ${pychop_dir}/${SAMPLE_NAME}_pychopped.fastq \
        2> ${LOGS}/${SAMPLE_NAME}_Pychopper_stderr.txt

    PYCHOPPED=${pychop_dir}/${SAMPLE_NAME}_pychopped.fastq

else echo " :::: Skipping Pychopper ::::"
    PYCHOPPED=${INPUT_FASTQ}
fi
##################
if [ ${ADAPTER_TRIM} = True ]; then
    echo "$(timestamp)"
    echo " :::: Trimming Adapters ::::"
    mkdir ${OUTFOLDER}/CUTADAPT

    echo "Extracting 10X Amplicon (R1,TSO) of Reads"
    cutadapt \
        -g file:"${DUAL_ADAPTER_10X}" \
        --cores=0 \
        -e 0.2 \
        --action trim \
        --match-read-wildcards \
        -o ${OUTFOLDER}/CUTADAPT/${SAMPLE_NAME}_Cutadapt_dual_trim.fastq \
        ${PYCHOPPED} \
        > ${LOGS}/Cutadapt_Amplicon_report_${SAMPLE_NAME}.txt 

    DUAL_TRIMMED=${OUTFOLDER}/CUTADAPT/${SAMPLE_NAME}_Cutadapt_dual_trim.fastq

    echo "Trimming remaining R1 or TSO on either Side Reads for chimeric Reads/PCR Artifacts"
    cutadapt \
        -b file:"${Adapter_5_3_10X}" \
        --cores=0 \
        --action trim \
        --match-read-wildcards \
        -o ${OUTFOLDER}/CUTADAPT/${SAMPLE_NAME}_Cutadapt_trimmed.fastq \
        ${DUAL_TRIMMED} \
        > ${LOGS}/Cutadapt_trimmed_${SAMPLE_NAME}.txt

    TRIMMED=${OUTFOLDER}/CUTADAPT/${SAMPLE_NAME}_Cutadapt_trimmed.fastq

else echo " :::: Skipping Trimming Adapters ::::"
    TRIMMED=${INPUT_FASTQ}
fi

################################################################
################## IgBlast BLOCK #######################
################################################################

if [ ${IGBLAST} = True ]; then
    echo " :::: Quering (trimmed) Input to IgBlast for vdj Clustering :::: "
    mkdir ${OUTFOLDER}/IGB_Trimmed
    mkdir ${OUTFOLDER}/IGB_Trimmed/TEMP_${SAMPLE_NAME}
    cd ${OUTFOLDER}/IGB_Trimmed

    pyir \
        -t fastq \
        -m ${THREADS} \
        --outfmt tsv \
        --pretty \
        -r TCR \
        -s human \
        --numV 1 \
        --numD 1 \
        --numJ 1 \
        --tmp_dir ${OUTFOLDER}/IGB_Trimmed/TEMP_${SAMPLE_NAME} \
        -o ${SAMPLE_NAME}_preprocessed_IGB \
        ${TRIMMED}

    gunzip ${SAMPLE_NAME}_preprocessed_IGB.tsv.gz

    INPUT_IGB=${SAMPLE_NAME}_preprocessed_IGB.tsv

    rmdir ${OUTFOLDER}/IGB_Trimmed/TEMP_${SAMPLE_NAME}
else 
    echo " :::: Not quering IgBlast as indicated ::::"
fi
