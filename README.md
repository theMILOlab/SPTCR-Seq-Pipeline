# SPTCR-seq
Pipeline to demultiplex, correct , annotate and UMI correct T-Cell Receptor Long Reads.


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

        Demultiplexing Pipeline that matches the Barcodes to the long Reads. The script als extracts the UMI Region from the Long Read by Substracting the Strings Adapter-seq+16bp - Adapter-seq+28bp.

            ./1_Demultiplex_UMI_Extraction.sh
                usage: 1_Demultiplex_UMI_Extraction.sh [-h] [-n NAME] -i INPUT_FASTQ [-o OUTFOLDER]
                                                    [-t THREADS] [-mem MEMORY] [-rep REPOSITORY]
                                                    [-a ADAPTER]
                -i INPUT_FASTQ  Path to the Input Fastq
                -rep REPOSITORY Path to the cloned Github Repository, default: ../
                -n NAME         Sample Name, if not specified it will default to the basename of the fastq_%d_%Y
                -a ADAPTER      Sequence of the R1 Adapter used for Library Preparation default: 'CTACACGACGCTCTTCCGATCT'
                -t THREADS      Number of Threads to use
                -mem MEMORY     Maximum Memory to use for Barcode Matching

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
        Preprocess the Reads for the Correction Pipeline. This script splits fusioned Reads, trims the Adapters from 10X & the aligns it for the TCR-sequences.

        ```
        ./2_Preprocess_Reads.sh
        usage: 2_Preprocess_Reads.sh [-h] [-n NAME] -i INPUT_FASTQ [-o OUTFOLDER] [-t THREADS]
                                [-mem MEMORY] [-rep REPOSITORY] [-pri PRIMER]
                                [-conf CONFIGURATION] [-chop PYCHOPPER] [-trim ADAPTER_TRIM]
                                [-igb IGBLAST]
                -rep REPOSITORY         Path to the cloned Github Repository, default: ../
                -conf CONFIGURATION     Path to the Primer Configuration File for PyChopper, default: ./REFERENCE/Primer/Pychopper/                     
                                        10XPrimers_pychopper_configuration.txt
                -pri PRIMER             Path to the Primer .fa File for PyChopper, default: ./REFERENCE/Primer/Pychopper/                     
                                        10XPrimers_pychopper.fa
                -n NAME                 Sample Name, if not specified it will default to the basename of the fastq_%d_%Y
                -a ADAPTER              Sequence of the R1 Adapter used for Library Preparation default: 'CTACACGACGCTCTTCCGATCT'
                -t THREADS              Number of Threads to use
                -mem MEMORY             Maximum Memory to use for Barcode Matching
                -chop PYCHOPPER         Specify if Reads should be made full length by Pychopper, default: True
                -igb IGBLAST            If True, the preprocessed Fastq is aligned with IgBLAST Following Processing., default: True
                -trim ADAPTER_TRIM      Specify if Reads should be from Adapters