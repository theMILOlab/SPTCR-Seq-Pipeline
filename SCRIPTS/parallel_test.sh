#! /usr/bin/env bash 

CORRECTION_INPUT=/mnt/681ABFBA1ABF839A/Dropbox/KBJasim/Projects/Capture_Sequencing/LONG_TCR/FINAL/Parallel_Test/INPUT
RATTLE_PATH=/mnt/681ABFBA1ABF839A/Dropbox/KBJasim/Projects/Capture_Sequencing/LONG_TCR/FINAL/Repository/TOOLS/RATTLE
LOGS=/mnt/681ABFBA1ABF839A/Dropbox/KBJasim/Projects/Capture_Sequencing/LONG_TCR/FINAL/Parallel_Test/LOGS
OUTDIR=/mnt/681ABFBA1ABF839A/Dropbox/KBJasim/Projects/Capture_Sequencing/LONG_TCR/FINAL/Parallel_Test/OUT

THREADS=24

################################################################
echo "-------- First Round Cutadapt with Dual Adapters ----- "
find ${pychop_dir} -type f |parallel --jobs 0 "cutadapt \
        -g file:"${DUAL_ADAPTER_10X}" \
        --cores=0 \
        -e 0.2 \
        --action trim \
        --match-read-wildcards \
        -o ${TEMPORARY}/DUAL_TRIMMED/{/.}_trimmed_R1.fastq \
        {} \
        > ${LOGS}/CUTADAPT_LOGS/cutadapt_R1_dual_report_{/.}.txt" 

################################################################
################################################################

        #CLUSTERNAME=${FILE##*/}
        #CLUSTERNAME=${CLUSTERNAME%.fastq}
        #echo ${CLUSTERNAME}
         FILE={}
        echo ${FILE}
echo "-------- Parallel Correction of all Clusters ----- "

RATTLE_CLUSTERS=${OUTDIR}/CORRECTION/RATTLE_CLUSTERS
RATTLE_CORRECT_OUT=/mnt/681ABFBA1ABF839A/Dropbox/KBJasim/Projects/Capture_Sequencing/LONG_TCR/FINAL/Parallel_Test/OUT/CORRECTION/CORRECT


find ${CORRECTION_INPUT} -type f |parallel --jobs 0 \
        "
        mkdir ${RATTLE_CLUSTERS}/{/.}
        ${RATTLE_PATH}/rattle cluster \
            -i {} \
            -t ${THREADS} \
            -o ${RATTLE_CLUSTERS}/{/.} \
            --iso \
            --fastq \
            2> ${LOGS}/Rattle_Clustering/{/.}_clust_stderr.txt

        mkdir ${RATTLE_CORRECT_OUT}/{/.}
        ${RATTLE_PATH}/rattle correct \
            -i {} \
            -c ${RATTLE_CLUSTERS}/{/.}/clusters.out \
            -t ${THREADS} \
            -o ${RATTLE_CORRECT_OUT}/{/.} \
            2> ${LOGS}/Rattle_Correction/{/.}_corr_stderr.txt
        "

################################################################
################################################################
            



find ${CORRECTION_INPUT} -type f |parallel --jobs 0 \
        "
        ${RATTLE_PATH}/rattle cluster_summary \
            -i {} \
            -c ${RATTLE_CLUSTERS}/{/.}/clusters.out
            > ${RATTLE_CLUSTERS}/{/.}/{/.}_clust_summary.csv"


for cluster in ${RATTLE_CLUSTERS}/* ; do
        clustername="$(basename ${cluster})"
        clustername="$(cut -d'.' -f1 <<<${clustername})"
        ${RATTLE_PATH}/rattle cluster_summary \
            -i  ${CORRECTION_INPUT}/${clustername}.fastq \
            -c ${cluster}/clusters.out \
            > ${cluster}/${clustername}_clust_summary.csv
done

################################################################
find ${CORRECTION_INPUT} -type f |parallel --jobs 0 \
        "
        mkdir ${RATTLE_CLUSTERS}/{/.}_2
        ${RATTLE_PATH}/rattle cluster \
            -i {} \
            -o ${RATTLE_CLUSTERS}/{/.}_2 \
            -t ${THREADS} \
            --iso \
            --fastq \
            2> ${LOGS}/Rattle_Clustering/{/.}_clust_stderr2.txt"




find ${CORRECTION_INPUT} -type f |parallel --jobs 0 \
        "
        mkdir ${RATTLE_CORRECT_OUT}/{/.}
        ${RATTLE_PATH}/rattle correct \
            -i {} \
            -c ${RATTLE_CLUSTERS}/{/.}/clusters.out \
            -t ${THREADS} \
            -o ${RATTLE_CORRECT_OUT}/{/.} \
            2> ${LOGS}/Rattle_Correction/{/.}_corr_stderr.txt
        "


        
        mkdir ${RATTLE_CLUSTERS}/

        ${RATTLE_PATH}/rattle cluster \
            -i ${FILE} \
            -t ${THREADS} \
            -o ${RATTLE_CLUSTERS} \
            --iso \
            --fastq \
            > ${LOGS}/Rattle_Clustering/${CLUSTERNAME}.txt 2> ${LOGS}/stderr_Rattle_Clustering/${CLUSTERNAME}.txt

        ${RATTLE_PATH}/rattle correct \
            -i {} \
            -c ${RATTLE_CLUSTERS}/clusters.out \
            -t ${THREADS} \
            -o ${OUTDIR}/CORRECTION/RATTLE_Correction \
            > ${LOGS}/Rattle_Correction/${CLUSTERNAME}_correction_txt_output.txt 2> ${LOGS}/Rattle_Correction/${CLUSTERNAME}_correction_txt_output.txt




################################################################
    ${RATTLE_PATH}/rattle cluster \
        -i "${SAMPLE}" \
        -t ${THREADS} \
        -o "${OUTDIR}"/CORRECTION/RATTLE_CLUSTERS \
        --iso \
        --fastq \
        > ${LOGS}/Rattle_Clustering_${CLUSTERNAME}.txt 2> ${LOGS}/stderr_Rattle_Clustering_${CLUSTERNAME}.txt

    RATTLE_CLUSTERS="${OUTDIR}"/CORRECTION/RATTLE_CLUSTERS/clusters.out

    ${RATTLE_PATH}/rattle correct \
        -i "${SAMPLE}" \
        -c ${RATTLE_CLUSTERS} \
        -t ${THREADS} \
        -o "${OUTDIR}"/CORRECTION/RATTLE_Correction \
        > ${LOGS}/${CLUSTERNAME}_correction_txt_output.txt 2> ${LOGS}/${CLUSTERNAME}_correction_txt_output.txt


--colsep regexp Split input on regexp for positional replacements
{} {.} {/} {/.} {#} {%} {= perl code =} Replacement strings
{3} {3.} {3/} {3/.} {=3 perl code =}    Positional replacement strings
With --plus:    {} = {+/}/{/} = {.}.{+.} = {+/}/{/.}.{+.} = {..}.{+..} =
                {+/}/{/..}.{+..} = {...}.{+...} = {+/}/{/...}.{+...}
