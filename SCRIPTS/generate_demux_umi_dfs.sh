#!/usr/bin/env bash


##########################################
##### Argparse Options #####
ARGPARSE_DESCRIPTION="Shell File to generate the Files needed for demultiplexing"      # this is optional
source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-i','--input',help='Input Raw Fastq',required=True)
parser.add_argument('-o','--output',help='Outfolder',required=True)
parser.add_argument('-n','--samplename',help='Sample Name',required=True)
parser.add_argument('-a', '--adapter', default='CTACACGACGCTCTTCCGATCT', 
                    help='Specify Read 1 Sequence')
parser.add_argument('-b', '--barcodes', default='../REFERENCE/Barcodes/visium_bc.tsv', 
                    help='Whitelist of used Barcodes')
parser.add_argument('-s', '--scTagger', default='', 
                    help='Path to scTagger Scripts')
parser.add_argument('-m', '--memory',default=16,help='Specify RAM')
parser.add_argument('-t', '--threads',default='2',
                    help='Specify Number of threads to use')
EOF


echo ":::::Demultiplexing and UMI Extraction ${samplename}:::::"
echo ${samplename}
echo ${input}
mkdir ${output}/${samplename}_UMI_Extraction

echo "::: Extracting Adapter and Barcode :::"
${scTagger}/extract_lr_bc.py \
    -r ${input} \
    -t ${threads} \
    -sa ${adapter} \
    --num-bp-after 16 \
    -o ${output}/${samplename}_UMI_Extraction/${samplename}_Extracted_Adapter_Barcode.tsv.gz

echo "::: Extracting Adapter, Barcode and UMI:::"
${scTagger}/extract_lr_bc.py \
    -r ${RAW_FASTQ} \
    -t ${threads} \
    -sa ${adapter} \
    --num-bp-after 28 \
    -o ${output}/${samplename}_UMI_Extraction/${samplename}_Extracted_Adapter_Barcode_UMI.tsv.gz

################
mkdir ${output}/${samplename}_DEMUX
echo "::::: Extracting Long Read Adapter Segments :::::"
${scTagger}/extract_lr_bc.py \
    -r ${RAW_FASTQ} \
    -t ${threads} \
    -sa ${adapter} \
    -o ${output}/${samplename}_DEMUX/${samplename}_lr_bc.tsv.gz

echo "::::: Matching Barcodes to Extracted Segments :::::"
${scTagger}/match_lr_bc-trie.py \
    -lr ${output}/${samplename}_DEMUX/${samplename}_lr_bc.tsv.gz \
    -sr ${VISIUM_BCs}  \
    -o ${output}/${samplename}_scTagger_DEMUX.tsv.gz \
    -t ${threads} \
    -m ${memory} 


