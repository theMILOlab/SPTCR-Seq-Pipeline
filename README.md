# SPTCR-seq
Pipeline to demultiplex, correct , annotate and UMI correct T-Cell Receptor Long Reads.

![image](https://user-images.githubusercontent.com/70334482/175873404-d5ef14b1-5be4-4789-8ae9-5214b924a89e.png)


## Installation

We recommend using the Conda C++ Drop-In package manager mamba to resolve all the Dependencies faster. Just install as followed and subsequently, simply change installation commands from 'conda install ...' to 'mamba install ...'.

```
conda install mamba -n base -c conda-forge
```

For Installation:
1. Clone this Repository with
        ```
        git clone https://github.com/jkbenotmane/SPTCR-seq.git
        ```

2. Install Requirements with
        ```
        cd ./SPTCR-seq
        conda env create -f SPTCR_ENV.yml
        ```
3. To do Demultiplexing do:
        ```
        cd ./Tools
        git clone https://github.com/vpc-ccg/scTagger.git
        ```
4. For Error Correction, do:
        ```
        cd./Tools
        git clone --recurse-submodules https://github.com/comprna/RATTLE
        cd RATTLE
        ./build.sh
        ```

5. Activate Environment to run the Pipeline
        ```
        conda activate -n SPTCR_ENV
        ```

## Running Pipeline

### 1. Demultiplexing Reads
***
        Demultiplexing Pipeline that matches the Barcodes to the long Reads. The script als extracts the UMI Region from the Long Read by Substracting the Strings Adapter-seq+16bp - Adapter-seq+28bp.

        ./1_Demultiplex_UMI_Extraction.sh
        usage: 1_Demultiplex_UMI_Extraction.sh [-h] [-n NAME] -i INPUT_FASTQ [-o OUTFOLDER]
                                                    [-t THREADS] [-mem MEMORY] [-rep REPOSITORY]
                                                    [-a ADAPTER]
        -i INPUT_FASTQ          Path to the Input Fastq
        -rep REPOSITORY         Path to the cloned Github Repository, default: ../
        -n NAME                 Sample Name, if not specified it will default to the basename of the fastq_%d_%Y
        -a ADAPTER              Sequence of the R1 Adapter used for Library Preparation default: 'CTACACGACGCTCTTCCGATCT'
        -t THREADS              Number of Threads to use
        -mem MEMORY             Maximum Memory to use for Barcode Matching

        #### Example:

        After Running the Pipeline you will get in the specified Outfolder following Table:


                ./OUTFOLDER/YOUR_SAMPLE_NAME_barcode_umi.csv

                Spatial Barcode,ReadID,UMI
                TGGGAGATACTTGTTT,0174f20e-b476-4ec3-a51f-a6895fcb846e,TTTGCAAGATTT
                TGCTAGTAGATTGTTT,4224b8ac-ba13-49df-b862-96746132b926,TCACACCCGTTC
                TGCTAGTAGATTGTTT,ec527ee2-49e7-43c1-bd07-b9a4425e22c5,TTTCGGGACGAT
                TGCTAGTAGATTGTTT,753a2c95-b76d-46c6-b080-04fa591b7eab,GATTGTTTGTTT
                TGCTAGTAGATTGTTT,37114111-fc84-4472-b06c-5bc67a85e615,GCGTGTCGCCAT

         As Barcodes are matched by String Distance, some reads have multiple equal matches. For the final Output Table above we just choose the first one. To see all look at:

                ./OUTFOLDER/YOUR_SAMPLE_NAME_DEMUX/YOUR_SAMPLE_NAME_barcode_umi.csv

                ReadID,Barcode,0,1,2
                1e3bc9e7-6e55-4242-806d-d6eae46b621f,GGTGGGATTAGGTCCC,GGTGGGATTAGGTCCC,,
                f8e40593-e913-4259-9a42-9266a7c14f09,TCCTTACGACGGTCCG,TCCTTACGACGGTCCG,,
                1c3e69d7-8858-43df-a329-f3f551fbf1e4,AACCAGACTCTTCGGT,AACCAGACTCTTCGGT,,
                b2e47859-226a-4446-9c86-3c1999fc90f7,AGTGATGACGCGTATA,AGTGATGACGCGTATA,,
                bba4f060-27ae-4eba-9566-cdc7a92a4d17,CGTAGCAGTAATGGAC,CGTAGCAGTAATGGAC,,


### 2. Preprocess Reads
***
        Preprocess the Reads for the Correction Pipeline. This script splits fusioned Reads, trims the Adapters from 10X & aligns the Reads for the TCR Annotations. If the Barcode UMI csv from Step 1 will be used to demultiplex the annotated TCR sequences.

        ./2_Preprocess_Reads.sh
        usage: 2_Preprocess_Reads.sh [-h] [-n NAME] -i INPUT_FASTQ [-o OUTFOLDER] [-t THREADS]
                                [-mem MEMORY] [-rep REPOSITORY] [-pri PRIMER]
                                [-conf CONFIGURATION] [-chop PYCHOPPER] [-trim ADAPTER_TRIM]
                                [-igb IGBLAST]

        -i INPUT_FASTQ          Path to the Input Fastq
        -o OUTFOLDER            Directory for the Outfolder, default: PWD
        -rep REPOSITORY         Path to the cloned Github Repository, default: ../
        -conf CONFIGURATION     Path to the Primer Configuration File for PyChopper, default: ./REFERENCE/Primer/Pychopper/                     
                                10XPrimers_pychopper_configuration.txt
        -pri PRIMER             Path to the Primer .fa File for PyChopper, default: ./REFERENCE/Primer/Pychopper/                     
                                10XPrimers_pychopper.fa
        -n NAME                 Sample Name, if not specified it will default to the basename of the fastq_%d_%Y
        -a ADAPTER              Sequence of the R1 Adapter used for Library Preparation default: 'CTACACGACGCTCTTCCGATCT'
        -t THREADS              Number of Threads to use
        -chop PYCHOPPER         Specify if Reads should be made full length by Pychopper, default: True
        -igb IGBLAST            If True, the preprocessed Fastq is aligned with IgBLAST Following Processing., default: True
        -trim ADAPTER_TRIM      Specify if Reads should be from Adapters

        #### Examplary Output:
                                                        Locus	V	D	J	CDR3	CDR3_aa	Spatial Barcode	UMI
        ReadID								
        87f02d5d-5e89-4fc7-971f-2e5b2d77f196	TRB	TRBV3-1*01	NaN	TRBJ1-2*01	NaN	NaN	CGGTTCAAGTAGGTGT	NaN
        c5e710cb-eefb-42b0-a02c-8f3a863bfac9	TRA	TRBV30*01	NaN	NaN	NaN	NaN	TTCGCACTGTACGACA	NaN
        7cf5db8f-1794-4bdf-b3f0-bf9f9e8d9f72	TRB	TRBV6-8*01	TRBD1*01	TRBJ1-4*01	GCCAGTGTGTCCGTGCTGGGGACGCCTCTGTCCAGGTGGTGAGTGT...	ASVSVLGTPLSRW*VWNGVCVCLCACYQRKT	AGGCCTATCAGGTACG	AATCGTTGTGTT
        c24fd351-a4ec-4156-85dd-69957d5b3d93	TRG	TRBV10-3*01	NaN	NaN	NaN	NaN	ACGGCCAACATGGACT	NaN

### 3. Cluster and Correct Reads
***
        Clusters Reads basd on their VJ Annotation generated by IgBlast Alignment, next Reads are corrected using Rattle Algorithm.
        
        ./3_Cluster_Correct.sh"
        usage: 3_Cluster_Correct.sh [-h] [-n NAME] -i INPUT_FASTQ -b INPUT_IGB [-o OUTFOLDER]
                                [-t THREADS] [-mem MEMORY] [-rep REPOSITORY] [-pri PRIMER]
                                [-conf CONFIGURATION] [-igb IGBLAST]

        -i INPUT_FASTQ          Path to the Input Fastq
        -i INPUT_IGB            Path to the Input IgBlast Result generated by 2_Preprocess_Reads.sh
        -n NAME                 Sample Name, if not specified it will default to the basename of the fastq_%d_%Y
        -o OUTFOLDER            Directory for the Outfolder, default: PWD
        -t THREADS              Number of Threads to use
        -igb IGBLAST            If True, the preprocessed Fastq is aligned with IgBLAST Following Processing., default: True
        -rep REPOSITORY         Path to the cloned Github Repository, default: ../

### Exemplary Pipeline
***
        An Exemplary minimal Usage Pipeline to demultiplex, preprocess and correct an Input Fastq.

        ```
        #! /usr/bin/env bash
        
        ## Command into Github Repository

        demux=./1_Demultiplex_UMI_Extraction.sh
        preproc=./2_Preprocess_Reads.sh
        clusco=./3_Cluster_Correct.sh

        ##############################

        ${demux}     \
        -i ./test_demux.fastq 

        ${preproc} \
        -i ./test_demux.fastq 

        ${clusco} \
        -i ./test_demux.fastq \
        -b ./PreProcessing/IGB_Trimmed/test_demux_preprocessed_IGB.tsv

        ```
