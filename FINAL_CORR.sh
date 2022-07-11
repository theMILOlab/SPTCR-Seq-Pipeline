#! /usr/bin/env bash 

full="/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Github SPTCR-seq Pipeline/SPTCR-Seq-Pipeline/SPTCR-seq Pipeline Scripts/SPTCR_FULL_Pipeline.sh"

sample="/media/jkbuntu/JKB_500GB/Raw_Nanopore/Example_TCR.fastq"


SAMPLE_NAME="$(basename "${sample}")"
SAMPLE_NAME="$(cut -d'_' -f1 <<<${SAMPLE_NAME})"

OUT="/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Github SPTCR-seq Pipeline/SPTCR-Seq-Pipeline/Example/${SAMPLE_NAME}"
echo "$SAMPLE_NAME"
#mkdir "${OUT}"
#cd "${OUT}"

############################
"${full}" \
    -n "$SAMPLE_NAME" \
    -i"${sample}" \
    -t 24 \
    -mem 60 \
    -rep "/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Github SPTCR-seq Pipeline/SPTCR-Seq-Pipeline" \
    -cln False 


exit
for sample in "/media/jkbuntu/JKB_500GB/Raw_Nanopore/01.SPTCR12_splits_sana.fastq"
do  


    ## Copy Work to Harddrive and remove
    #cp -R "${OUT}" "/media/jkbuntu/WD12TB/FINAL_PUBLISH_CORRECTION"
    #rm -r "${OUT}"
    
done

exit


##############################################################
for sample in "/media/jkbuntu/JKB_500GB/Raw_Nanopore"/SPTCR16_raw.fastq
do  
    SAMPLE_NAME="$(basename "${sample}")"
    SAMPLE_NAME="$(cut -d'_' -f1 <<<${SAMPLE_NAME})"

    OUT="/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Testing_Zone/Final/${SAMPLE_NAME}"
    echo "$SAMPLE_NAME"
    mkdir "${OUT}"
    cd "${OUT}"

    ############################
    "${full}" \
    -n "$SAMPLE_NAME" \
    -i"${sample}" \
    -t 24 \
    -mem 60 \
    -rep "/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Github SPTCR-seq Pipeline/SPTCR-Seq-Pipeline" \
    -cln False 

    ## Copy Work to Harddrive and remove
    cp -R "${OUT}" "/media/jkbuntu/WD12TB/FINAL_PUBLISH_CORRECTION"
    rm -r "${OUT}"
    
done

###################################

declare -a nayrray=( "01.SPTCR12" "SPTCR16" )

for sample in "/media/jkbuntu/JKB_500GB/Raw_Nanopore"/*.fastq
do  
    SAMPLE_NAME="$(basename "${sample}")"
    SAMPLE_NAME="$(cut -d'_' -f1 <<<${SAMPLE_NAME})"

    OUT="/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Testing_Zone/Final/${SAMPLE_NAME}"
    echo "$SAMPLE_NAME"

        ############################
    if [[ ! " ${nayrray[*]} " =~ "${SAMPLE_NAME}" ]]; then
        mkdir "${OUT}"
        cd "${OUT}"

        "${full}" \
        -n "$SAMPLE_NAME" \
        -i"${sample}" \
        -t 24 \
        -mem 60 \
        -rep "/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Github SPTCR-seq Pipeline/SPTCR-Seq-Pipeline" \
        -cln False 

        ## Copy Work to Harddrive and remove
        cp -R "${OUT}" "/media/jkbuntu/WD12TB/FINAL_PUBLISH_CORRECTION"
        rm -r "${OUT}"
    else
            echo "$SAMPLE_NAME Did Correction already"
    
    fi
done
