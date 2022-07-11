# ONLY FOR INTERNAL PURPOSES (YET) REPOSITORY CONTAINS in ./TOOLS THIRD PARTY CODE


# SPTCR-seq
Explore spatially resolved T-Cell Infiltration at high resolution with Oxford Nanopore Sequencing. This Pipeline is used to demultiplex, correct , annotate and UMI correct target enriched T-Cell Receptor cDNA acquired with the SPTCR-seq Protocol.

#### Pipeline Overview
![image](https://user-images.githubusercontent.com/70334482/175873404-d5ef14b1-5be4-4789-8ae9-5214b924a89e.png)

#### TCR Infiltration Analysis in spatial Context
##### Sample 275 TRA/TRB Spatial T-Cell Deconvolution
<img src="https://github.com/theMILOlab/SPTCR-Seq-Pipeline/blob/main/Images/TCR_Expression.png">

##### Sample 275 Spatial TRD/TRG Infiltration
<img src="https://github.com/theMILOlab/SPTCR-Seq-Pipeline/blob/main/Images/GBM275_TRD_TRG_Infiltration.png" width="300">

## Installation

We recommend using the Conda C++ Drop-In package manager mamba (https://github.com/mamba-org/mamba) to resolve all the Dependencies faster. 
Execute the following Installation steps from conda base environment to automatically install the package manager. 
If you want to do installation by hand, just install as followed and subsequently, simply change your installation commands from 'conda install ...' to 'mamba install ...'.

```
conda install mamba -n base -c conda-forge
```

**For Installation:**
1. Clone this Repository with
   
        git clone https://github.com/theMILOlab/SPTCR-Seq-Pipeline.git

2. For parallel processing in correction step do:
   
        sudo apt-get install parallel   

3. Execute setup.sh to setup an Environment "SPTCR_ENV" with all Dependencies and build Tools from Source:

        cd ./SPTCR-Seq-Pipeline
        ./setup.sh

> !! If you have problems compiling RATTLE (especially contained spoa) from source see ./TOOLS/change_c++ versions.txt for some notes on Installation. RATTLE needs GCC/ G++ 9 to compile, see guide on how to maintain multiple compiler and c++ versions on your computer and to compile RATTLE. Also check the issue section of RATTLE (https://github.com/comprna/RATTLE) !!
     

4. Activate Environment to run the Pipeline
        conda activate -n SPTCR_ENV

5. For minimal user intervention Pipeline see Exemplary Pipeline: example.sh 

---
## Running Pipeline

### Overview: 
The Computational Workflow is compromised of two major steps: 

        1. Preprocessing Reads (2_Preprocess_Reads.sh)
        2. Cluster & Correct Reads (3_Cluster_Correct.sh)

When Calling SPTCR_Full_Pipeline.sh with the indicated args, these scripts will be called sequentially with standard args. This is the recommended form to call SPTCR_Pipeline.
If you want more control of the intermediate steps or reuse already calculated preprocessed reads, then you have to call the scripts sequentially.

### 1. Combined Demultiplexing & Preprocessing of Reads
***    
Performs preprocessing of the Reads for Correction & matching the Barcodes to the raw sequencing Result as well as generate a table of demultiplexed, annotated T-Cell Receptor Sequences and their adjoining UMI Region you can do with.:

#### Usage

        usage: 2_Preprocess_Reads.sh [-h] [-n NAME] -i INPUT_FASTQ [-o OUTFOLDER] [-t THREADS]
                                [-mem MEMORY] [-rep REPOSITORY] [-pri PRIMER] [-conf CONFIGURATION]
                                [-chop PYCHOPPER] [-trim ADAPTER_TRIM] [-igb IGBLAST]
                                [-demux DEMULTIPLEX]

        Pipeline to preprocess, demultiplex and extract UMI regions of Nanopore Reads for Libraries
        prepared for 10X Genomics

**Arguments**

        -h, --help            show this help message and exit
        -n NAME, --NAME NAME  Chosen Name for the Pipeline. Will Be Used as Name for the Pipelines
                                Outfolder.
        -i INPUT_FASTQ, --INPUT_FASTQ INPUT_FASTQ
                                Specify the Path to the raw Input Fastq File
        -o OUTFOLDER, --OUTFOLDER OUTFOLDER
                                Specify the Directory for the Output Folder. If not specified, will use
                                current working Directory.
        -t THREADS, --THREADS THREADS
                                Number of Threads to use.
        -mem MEMORY, --MEMORY MEMORY
                                Gigabytes of RAM to use for Demultiplexing
        -rep REPOSITORY, --REPOSITORY REPOSITORY
                                Specify the Location of the Github Repository Folder for SPTCR Seq
        -pri PRIMER, --PRIMER PRIMER
                                Specify Custom Primers if not having used either 10X Visium or Single Cell
                                for the reconstruction of Full Reads by Pychopper.
        -conf CONFIGURATION, --CONFIGURATION CONFIGURATION
                                Specify the possible Configurations of the Custom Primers for Pychopper if
                                not having used either 10X Visium or Single Cell for the reconstruction of
                                Full Reads by Pychopper. See PyChoppers Documentation
                                (https://github.com/epi2me-labs/pychopper) for explanation
        -chop PYCHOPPER, --PYCHOPPER PYCHOPPER
                                Specify if Pychopper should be performed on Input. Will use Input Fastq as
                                Pychopped File.
        -trim ADAPTER_TRIM, --ADAPTER_TRIM ADAPTER_TRIM
                                Specify if Reads should be trimmed from Adapters. If set to False will use
                                PyChopper Output
        -igb IGBLAST, --IGBLAST IGBLAST
                                If True, the preprocessed Fastq is aligned with IgBLAST Following
                                Processing. If already done, use the Path to the IgBlast Output and skip
        -demux DEMULTIPLEX, --DEMULTIPLEX DEMULTIPLEX
                                If set to True, extracts Barcode and UMI Region of the Reads and updates
                                the IgBlast Table. Form is default for downstream purposes, it is
                                recommended to leave as default if you intend to correct the SPTCR-seq
                                reads as well.



#### Example
        NAME="TEST"
        INPUT_FASTQ="PATH/TO/INPUT/FASTQ"
        THREADS=12 ## Number of Threads Given
        MEMORY=16 ## Gigabytes Given for Demultiplexing Step
        REPOSITORY= PATH/TO/GITHUB/REPO/SPTCR-Seq-Pipeline

        bash ./2_Preprocess_Reads.sh \
                -n ${NAME} \
                -i "${INPUT_FASTQ}" \
                -t ${THREADS} \
                -mem ${MEMORY} \
                -rep "${REPOSITORY}"

#### Example Output
see Example/PreProcessing for exemplary output of the PreProcessing Pipeline.

**Explanation of Output**

***./SAMPLENAME/PreProcessing/***

*/LOGS/*

Holds the Log File for the PreProcessing and demultiplexing pipeline.

*./SAMPLENAME_Cutadapt_trimmed_sana.fastq & ./SAMPLENAME_Cutadapt_trimmed_sana.fastq.fxi*

PreProcessed Fastq & its Fastq adjoining index

*./SAMPLENAME_preprocessed_IGB.tsv*

Raw IgBlast Call of PreProcessed Reads

*/Demultiplexing_SAMPLENAME/*

Folder Created by 1_Demultiplex_UMI_Extraction.sh. Holds the demultiplexed IgBlast Output.

*./SAMPLENAME_all_barcode_matches.csv*

Holds all relevant Barcode Matches in given edit distance for read found by scTagger. By Default we simply choose the First given Match as Barcode. 

*./SAMPLENAME_barcode_umi.csv*

Serves as a Demultiplexing, Deduplication Table. Holds Columns: Spatial Barcode,ReadID,UMI for given Input Fastq.

*./SAMPLENAME_vdj_umi_barcode_uncorrected_df.csv*

Overview Table of the demultiplexed IgBlast File. Holds Columns: ReadID,Locus,V,D,J,CDR3,CDR3_aa,Spatial Barcode,UMI


### 1.1 Demultiplexing Reads
***
Demultiplexing Pipeline that matches the Barcodes to the long Reads. The script extracts the UMI Region from the Long Read by Substracting the Strings Adapter-seq+16bp - Adapter-seq+28bp. Deduplication happens after Correction with SCRIPTS/umi_correction_from_summary.py. Part of 2_Preprocess_Reads.sh but can be called externally. 

#### Usage
1_Demultiplex_UMI_Extraction.sh [-h] [-n NAME] -i INPUT_FASTQ [-igb INPUT_IGB]
                                       [-o OUTFOLDER] [-t THREADS] [-mem MEMORY] [-rep REPOSITORY]
                                       [-a ADAPTER]

Pipeline to barcode and extract UMI regions of ONT Reads for Libraries prepared for 10X Genomics

**Arguments**

        -h, --help            show this help message and exit
        -n NAME, --NAME NAME  Name of Output Folder
        -i INPUT_FASTQ, --INPUT_FASTQ INPUT_FASTQ
                        Specify the Path to the Raw unmodified Input Fastq File
        -igb INPUT_IGB, --INPUT_IGB INPUT_IGB
                        Specify the Path to IgBlast File to be demultiplexed with
                        demultiplex_summarize.py. If not specified will only output table for later
                        demultiplexing.
        -o OUTFOLDER, --OUTFOLDER OUTFOLDER
                        Specify the Directory for the Outputfolder
        -t THREADS, --THREADS THREADS
                        Number of Threads
        -mem MEMORY, --MEMORY MEMORY
                        RAM to user
        -rep REPOSITORY, --REPOSITORY REPOSITORY
                        Specify the Location of the Repositroy Folder holding all References and
                        scripts for SPTCR Seq
        -a ADAPTER, --ADAPTER ADAPTER
                        Specify Illumina Read 1 Sequence. Adapter is matched as Anchor, to
                        demultiplex and extract the UMI Region..

#### Example
        NAME="TEST"
        INPUT_FASTQ="PATH/TO/INPUT/FASTQ"
        THREADS=12 ## Number of Threads Given
        MEMORY=16 ## Gigabytes Given for Demultiplexing Step
        REPOSITORY= PATH/TO/GITHUB/REPO/SPTCR-Seq-Pipeline
        IGBLAST=/OUTFOLDER/PREPROCESSING/SAMPLENAME_preprocessed_IGB.tsv

        bash ./1_Demultiplex_UMI_Extraction.sh \
                -i "${INPUT_FASTQ}" \
                -igb "${IGBLAST}" \
                -n "${SAMPLE_NAME}" \
                -t ${THREADS} \
                -mem ${MEMORY} \
                -rep "${REPOSITORY}"

#### Example Output
see Example/PreProcessing for exemplary output of the PreProcessing Pipeline.

**Explanation of Output**
***./OUTFOLDER/Demultiplexing_SAMPLENAME/***

Folder Created by 1_Demultiplex_UMI_Extraction.sh. Holds the demultiplexed IgBlast Output.

*/LOGS/*

Holds the Log File for the PreProcessing and demultiplexing pipeline.

*./SAMPLENAME_all_barcode_matches.csv*

Holds all relevant Barcode Matches in given edit distance for read found by scTagger. By Default we simply choose the First given Match as Barcode. 

*./SAMPLENAME_barcode_umi.csv*

Serves as a Demultiplexing, Deduplication Table. Holds Columns: Spatial Barcode,ReadID,UMI for given Input Fastq.

*./SAMPLENAME_vdj_umi_barcode_uncorrected_df.csv*

Overview Table of the demultiplexed IgBlast File. Holds Columns: ReadID,Locus,V,D,J,CDR3,CDR3_aa,Spatial Barcode,UMI



### 3. Cluster and Correct Reads
***
        Clusters Reads based on their Locus and Variable Gene-Family Annotation generated by IgBlast Alignment, next Read-groups are parallel corrected using Rattle Algorithm and finally annotated by IGBlast.
        
#### Usage
3_Cluster_Correct.sh [-h] [-n NAME] -i INPUT_FASTQ -b INPUT_IGB [-o OUTFOLDER] [-t THREADS]
                            [-lm LOWMEM] [-rep REPOSITORY] [-igb IGBLAST] [-cln CLEANUP]
                            [-bc BCUMI] [-grp GROUPER]

Pipeline to group & Correct T-Cell Receptor Reads generated by Oxford Nanopore Reads of Libraries
prepared for 10X Genomics


**Arguments**

        -h, --help            show this help message and exit
        -n NAME, --NAME NAME  Name of Output Folder
        -i INPUT_FASTQ, --INPUT_FASTQ INPUT_FASTQ
                                Specify the Path to (preprocessed) Fastq File
        -b INPUT_IGB, --INPUT_IGB INPUT_IGB
                                Path to IgBlast.tsv. Use Full IgBlast Output from PreProcessed Reads for
                                the grouping of TCR Reads.
        -o OUTFOLDER, --OUTFOLDER OUTFOLDER
                                Specify the Directory for the Outputfolder, default = PWD
        -t THREADS, --THREADS THREADS
                                Number of Threads to use, default=2
        -lm LOWMEM, --LOWMEM LOWMEM
                                Set to True if Memory & Compute Intensive Parallel-Correction Step should
                                be done sequentially to reduce System Pressure. default=False
        -rep REPOSITORY, --REPOSITORY REPOSITORY
                                Specify the Location of the Repositroy Folder holding all References and
                                scripts for SPTCR Seq,default,default=./
        -igb IGBLAST, --IGBLAST IGBLAST
                                If True, the corrected Fastq is aligned with IgBLAST Following Correction.
                                default=True
        -cln CLEANUP, --CLEANUP CLEANUP
                                If True, created intermediate Files and Folders will be deleted. For
                                Debugging you can set this to False. default=False
        -bc BCUMI, --BCUMI BCUMI
                                Specify the path to the SAMPLENAME_barcode_umi.csv generated by the
                                preprocessing pipeline. Defaults to
                                OUTFOLDER/PreProcessing/SAMPLENAME_barcode_umi.csv
        -grp GROUPER, --GROUPER GROUPER
                                Grouper to be used by TCR_GROUPING_IGB.py for grouping the TCR Reads into
                                Fastqs. Pick one or a combination of: locus,v_family, d_family, j_family,
                                default= locus, v_family


#### Example

        NAME="TEST"
        THREADS=12 ## Number of Threads Given
        MEMORY=16 ## Gigabytes Given for Demultiplexing Step
        REPOSITORY= PATH/TO/GITHUB/REPO/SPTCR-Seq-Pipeline

        bash "./SPTCR-seq Pipeline Scripts/3_Cluster_Correct.sh" \
                -i "./PreProcessing/${NAME}_Cutadapt_trimmed_sana.fastq" \
                -b "./PreProcessing/${NAME}_preprocessed_IGB.tsv" \
                -n ${NAME} \
                -t ${THREADS} \
                --GROUPER "locus,v_family" \
                -rep "${REPOSITORY}" 

#### Example Output
see Example/ClusterCorrect for exemplary output of the PreProcessing Pipeline.

**Explanation of Output**

see Example/Example/ClusterCorrect for exemplary output of the Clustering and Correction Pipeline. Output for CLEANUP=False:

*Example_corr_igb_overview_igb.csv*

Overview Table generated for Convenience. Holds VDJ Information as well as Barcode, Umi per ReadID. Columns: ReadID,Locus,V,D,J,CDR3,CDR3_aa,Spatial Barcode,UMI

*Example_corrected_IGB.tsv*

Raw IgBlast Output of corrected Fastq.

*Example_corrected_merged.fastq*

Corrected Fastq.

*Example_igb_corr_umi_corrected.csv*

Overview Table of the corrected IGB Output Count summarized and UMI COrrected. Hold Columns: Spatial Barcode,Uncorrected Count,UMI Corrected,Locus,V,D,J,CDR3_aa

***ClusterCorrect/LOGS/***

Holds the Logs for each Step of the Pipeline

***ClusterCorrect/IGB_CLUSTERS/***

Holds the by their Arrangement split Fastqs of the Input Fastq.

***ClusterCorrect/CORRECTION/***

*RATTLE_CLUSTERS*

Holds all Folders with Output of RATTLE Clustering of each grouped Fastq.

*RATTLE_CORRECT*

Holds the Folders for each corrected Fastq Group.


### Exemplary Pipeline
***
        An Exemplary minimal Usage Pipeline to demultiplex, preprocess and correct an Input Fastq.

        #! /usr/bin/env bash 

        full="../SPTCR-seq Pipeline Scripts/SPTCR_FULL_Pipeline.sh"

        sample="./Example_TCR.fastq"

        SAMPLE_NAME="$(basename "${sample}")"
        SAMPLE_NAME="$(cut -d'_' -f1 <<<${SAMPLE_NAME})"

        OUT="../Example/${SAMPLE_NAME}"
        echo "$SAMPLE_NAME"

        ############################
        "${full}" \
        -n "$SAMPLE_NAME" \
        -i"${sample}" \
        -t 8 \
        -mem 16 \
        -rep ".." \
        -cln False 


### Citations

Ravi, V.M., Will, P., Kueckelhaus, J., Sun, N., Joseph, K., Salié, H., Vollmer, L., Kuliesiute, U., von Ehr, J., Benotmane, J.K., Neidert, N., Follo, M., Scherer, F., Goeldner, J.M., Behringer, S.P., Franco, P., Khiat, M., Zhang, J., Hofmann, U.G., Fung, C., Ricklefs, F.L., Lamszus, K., Boerries, M., Ku, M., Beck, J., Sankowski, R., Schwabenland, M., Prinz, M., Schüller, U., Killmer, S., Bengsch, B., Walch, A.K., Delev, D., Schnell, O., Heiland, D.H., 2022. Spatially resolved multi-omics deciphers bidirectional tumor-host interdependence in glioblastoma. Cancer Cell 40, 639-655.e13. doi:10.1016/j.ccell.2022.05.009

Martin, M., 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet j. 17, 10. doi:10.14806/ej.17.1.200

de la Rubia, I., Indi, J.A., Carbonell, S., Lagarde, J., Albà, M.M., Eyras, E., 2020. Reference-free reconstruction and quantification of transcriptomes from long-read sequencing. BioRxiv. doi:10.1101/2020.02.08.939942

Martin, M., 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet j. 17, 10. doi:10.14806/ej.17.1.200

tange_2022_6570228,
        author       = {Tange, Ole},
        title        = {GNU Parallel 20220522 ('NATO')},
        month        = May,
        year         = 2022,
        note         = {{GNU Parallel is a general parallelizer to run
                        multiple serial command line programs in parallel
                        without changing them.}},
        publisher    = {Zenodo},
        doi          = {10.5281/zenodo.6570228},
        url          = {https://doi.org/10.5281/zenodo.6570228
