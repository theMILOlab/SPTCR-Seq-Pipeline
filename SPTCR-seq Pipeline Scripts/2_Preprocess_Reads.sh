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
ARGPARSE_DESCRIPTION="Pipeline to preprocess, demultiplex and extract UMI regions of Nanopore Reads for Libraries prepared for 10X Genomics"
argparse "$@" <<EOF || exit 1

parser.add_argument('-n','--NAME',help="Chosen Name for the Pipeline. Will Be Used as Name for the Pipelines Outfolder.",default="-")
parser.add_argument('-i', '--INPUT_FASTQ', help="Specify the Path to the raw Input Fastq File",required=True)
parser.add_argument('-o', '--OUTFOLDER', help="Specify the Directory for the Output Folder. If not specified, will use current working Directory.", default="PWD")

parser.add_argument('-t','--THREADS',help="Number of Threads to use.", default="2")
parser.add_argument('-mem','--MEMORY',help="Gigabytes of RAM to use for Demultiplexing", default="8")

parser.add_argument('-rep', '--REPOSITORY', help="Specify the Location of the Github Repository Folder for SPTCR Seq",default="./")

parser.add_argument('-pri', '--PRIMER', help="Specify Custom Primers if not having used either 10X Visium or Single Cell for the reconstruction of Full Reads by Pychopper.",default="10X")
parser.add_argument('-conf', '--CONFIGURATION', help="Specify the possible Configurations of the Custom Primers for Pychopper if not having used either 10X Visium or Single Cell for the reconstruction of Full Reads by Pychopper. See PyChoppers Documentation (https://github.com/epi2me-labs/pychopper) for explanation")

parser.add_argument('-chop', '--PYCHOPPER', help="Specify if Pychopper should be performed on Input. Will use Input Fastq as Pychopped File.",default="True")
parser.add_argument('-trim', '--ADAPTER_TRIM', help="Specify if Reads should be trimmed from Adapters. If set to False will use PyChopper Output",default="True")
parser.add_argument('-igb', '--IGBLAST', help="If True, the preprocessed Fastq is aligned with IgBLAST Following Processing. If already done, use the Path to the IgBlast Output and skip",default="YES")
parser.add_argument('-demux', '--DEMULTIPLEX', help="If set to True, extracts Barcode and UMI Region of the Reads and updates the IgBlast Table. Form is default for downstream purposes, it is recommended to leave as default if you intend to correct the SPTCR-seq reads as well.",default="True")
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


### PYCHOPPER
if [ ${PRIMER} = 10X ];then
    PRIMERS="${REPOSITORY}"/Reference/Primer/Pychopper/10XPrimers_pychopper.fa
    PRIMER_CONFIGURATION="${REPOSITORY}"/Reference/Primer/Pychopper/10XPrimers_pychopper_configuration.txt
else
    PRIMER=${PRIMER}
    PRIMER_CONFIGURATION=${CONFIGURATION}
fi

##Cutadapt
DUAL_ADAPTER_10X="${REPOSITORY}"/Reference/Primer/Cutadapt/10X_Dual_Adapter.fa
Adapter_5_3_10X="${REPOSITORY}"/Reference/Primer/Cutadapt/5_3_10X_Adapter.fa

##Demultiplex
Demultiplex_UMI_Extraction="${REPOSITORY}/SPTCR-seq Pipeline Scripts/1_Demultiplex_UMI_Extraction.sh"

################################################################
################## Trim & Reorient BLOCK #######################
################################################################
echo " :::: Checking Fastq Integrity ::::"
seqkit sana "${INPUT_FASTQ}" -o "${INPUT_FASTQ}_sana" 2> "${LOGS}"/seqkit_sana.txt
INPUT_FASTQ="${INPUT_FASTQ}_sana"

if [ ${PYCHOPPER} = True ]; then
    echo "$(timestamp)" 
    echo " :::: Pychopper ::::"
    mkdir "${OUTFOLDER}"/PYCHOPPER
    pychop_dir="${OUTFOLDER}"/PYCHOPPER

    cdna_classifier.py \
        -m edlib \
        -b "${PRIMERS}" \
        -c "${PRIMER_CONFIGURATION}" \
        -r "${LOGS}"/"${SAMPLE_NAME}"_Pychopper_report.pdf \
        -S "${LOGS}"/"${SAMPLE_NAME}"_Pychopper_report.tsv \
        -t ${THREADS} \
        -p \
        "${INPUT_FASTQ}"  \
        "${pychop_dir}"/"${SAMPLE_NAME}"_pychopped.fastq \
        2> "${LOGS}"/"${SAMPLE_NAME}"_Pychopper_stderr.txt

    PYCHOPPED="${pychop_dir}"/"${SAMPLE_NAME}"_pychopped.fastq

else echo " :::: Skipping Pychopper ::::"
    PYCHOPPED="${INPUT_FASTQ}" 
fi

##################
if [ ${ADAPTER_TRIM} = True ]; then
    echo "$(timestamp)"
    echo " :::: Trimming Adapters with Cutadapt ::::"
    mkdir "${OUTFOLDER}"/CUTADAPT

    echo "Extracting 10X Amplicon (R1,TSO) of Reads"
    cutadapt \
        -g file:"${DUAL_ADAPTER_10X}" \
        --cores=0 \
        -e 0.2 \
        --action trim \
        --match-read-wildcards \
        -o "${OUTFOLDER}"/CUTADAPT/"${SAMPLE_NAME}"_Cutadapt_dual_trim.fastq \
        "${PYCHOPPED}" \
        > "${LOGS}"/Cutadapt_Amplicon_report_"${SAMPLE_NAME}".txt 

    DUAL_TRIMMED="${OUTFOLDER}"/CUTADAPT/"${SAMPLE_NAME}"_Cutadapt_dual_trim.fastq

    echo "Trimming remaining R1 or TSO on either Side Reads for chimeric Reads/PCR Artifacts"
    cutadapt \
        -b file:"${Adapter_5_3_10X}" \
        --cores=0 \
        --action trim \
        --match-read-wildcards \
        -o "${OUTFOLDER}"/CUTADAPT/"${SAMPLE_NAME}"_Cutadapt_trimmed.fastq \
        "${DUAL_TRIMMED}" \
        > "${LOGS}"/Cutadapt_trimmed_"${SAMPLE_NAME}".txt
    
    TRIMMED="${OUTFOLDER}"/CUTADAPT/"${SAMPLE_NAME}"_Cutadapt_trimmed.fastq

    echo " :::: Checking Fastq Integrity ::::"
    seqkit sana "${TRIMMED}" -o "${TRIMMED}_sana" 2> "${LOGS}"/seqkit_sana_trimmed.txt
    TRIMMED="${TRIMMED}_sana"

else echo " :::: Skipping Trimming Adapters ::::"
    TRIMMED="${PYCHOPPED}"
fi

################################################################
################## IgBlast BLOCK #######################
################################################################

if [ ${IGBLAST} = YES ]; then
    echo " :::: Quering (trimmed) Input to IgBlast for vdj Clustering :::: "
    mkdir "${OUTFOLDER}"/IGB_Trimmed
    mkdir "${OUTFOLDER}"/IGB_Trimmed/TEMP_"${SAMPLE_NAME}"
    cd "${OUTFOLDER}"/IGB_Trimmed

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
        --gzip False \
        --tmp_dir "${OUTFOLDER}"/IGB_Trimmed/TEMP_"${SAMPLE_NAME}" \
        -o "${SAMPLE_NAME}"_preprocessed_IGB \
        "${TRIMMED}"

    #gunzip "${OUTFOLDER}"/IGB_Trimmed/"${SAMPLE_NAME}"_preprocessed_IGB.tsv.gz

    #IGBLAST="${OUTFOLDER}"/IGB_Trimmed/"${SAMPLE_NAME}"_preprocessed_IGB.tsv

    ## Moving Output Files to the Front
    mv "${OUTFOLDER}"/IGB_Trimmed/"${SAMPLE_NAME}"_preprocessed_IGB.tsv "${OUTFOLDER}"/"${SAMPLE_NAME}"_preprocessed_IGB.tsv
    
    IGBLAST="${OUTFOLDER}"/"${SAMPLE_NAME}"_preprocessed_IGB.tsv

else 
    echo " :::: Not quering IgBlast as indicated ::::"

fi


echo " ::::: Cleaning Up :::::"
### Move Output Files to the Front
mv "${TRIMMED}" "${OUTFOLDER}"/"${SAMPLE_NAME}_Cutadapt_trimmed_sana.fastq"

### Remove Created Working Directories
rmdir "${OUTFOLDER}"/IGB_Trimmed/TEMP_"${SAMPLE_NAME}"
rmdir "${OUTFOLDER}"/IGB_Trimmed
rm -r "${OUTFOLDER}"/CUTADAPT
rm -r "${OUTFOLDER}"/PYCHOPPER




echo " ::::: Demultiplexing Output :::::"
if [ ${DEMULTIPLEX} = True ]; then
    cd "${OUTFOLDER}"
    
    "${Demultiplex_UMI_Extraction}" \
        -i "${INPUT_FASTQ}" \
        -igb "${IGBLAST}" \
        -n "${SAMPLE_NAME}" \
        -t ${THREADS} \
        -mem ${MEMORY} \
        -bc "${BARCODES}" \
        -rep "${REPOSITORY}"
else
    echo ":::: Demultiplexing set to False. If you use the Output for Downstream Applications, please match the Table Columns. ::::"
fi
