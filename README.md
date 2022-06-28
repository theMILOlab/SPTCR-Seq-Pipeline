# ONLY FOR INTERNAL PURPOSES (YET) REPOSITORY CONTAINS in ./TOOLS THIRD PARTY CODE



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

3. For parallel processing do:
        ```
        sudo apt-get install parallel
        ```
        
4. To do Demultiplexing do:
        ```
        cd ./Tools
        git clone https://github.com/vpc-ccg/scTagger.git
        ```

5. For Error Correction, do:
        ```
        cd./Tools
        git clone --recurse-submodules https://github.com/comprna/RATTLE
        ### Compile from Source
        cd RATTLE
        ./build.sh
        ```
        !! If you have problems compiling contained spoa from source see ./TOOLS/change_c++ versions.txt for some notes on how to maintain multiple compiler and c++ versions on your computer for compilation.!!

6. Activate Environment to run the Pipeline
        ```
        conda activate -n SPTCR_ENV
        ```

7. For minimal Pipeline see Exemplary Pipeline


## Running Pipeline

### 1. Combined Demultiplexing & Preprocessing of Reads
***    
        Preprocessing the Reads for Correction & matching the Barcodes to the raw sequencing Result as well as generate an table of demultiplexed, annotated T-Cell Receptor Sequences and their adjoining UMI Region you can do with.:

        2_Preprocess_Reads.sh \
                -n test_demux \
                --INPUT_FASTQ path/to/input \
                --DEMULTIPLEX True
                -rep "/Path/to/Repository"

        #### Example Output
        
        ./Outfolder/Demultiplexing_sample/sample_vdj_umi_barcode_uncorrected_df.csv
        Important Output Table used for downstream scripts or raw. Shows VDj Arrangement, spatial Barcode and UMI Region.

        *ReadID,Locus,V,D,J,CDR3,CDR3_aa,Spatial Barcode,UMI*
        c2a37d19-be55-4976-bb8a-eb23d4056d54,TRA,TRBV11-3*01,,,,,GGCGGGTCTAGCCAGG,TGCCTACTTATG
        a26264a7-94e0-40a0-b39e-e4c0cc3795ec,TRA,TRBV4-1*01,,,,,TTGTATCCCATGGTCT,CGGAACCCTATT
        cdc8f537-7e95-41bf-befc-076810e148a1,TRG,TRBV21/OR9-2*01,,,,,GTGTCCGGAGAGCAGC,GACCGGCTGATC
        cd8ab80f-fad8-417e-a8b7-052ebf783f34,TRG,TRBV9*01,,,,,GAAACCACGCGACCAT,GCACCGAATACT
        43953ab9-6de4-4a95-bd4e-5b3bba50312b,TRB,TRBV6-2*01,TRBD1*01,TRBJ2-3*01,,,TATCCTATCTTAACTA,AGATTTAAGGCT
        e9a448e2-4ceb-42c8-b390-a493ad5f67ad,TRB,TRBV28*01,TRBD2*01,TRBJ1-1*01,,,CGCAGGCGATCCAAAC,AACGCACTCACT
        a0895494-d9e0-4e33-ab46-4bb16f63361f,,,,,,,GATTTCACGGCGGCTT,CCGACATCATTA
        48bc4795-3660-4fd2-b5d2-bea672ca5dde,TRA,TRBV17*01,,,,,GTCAGACGCTGACGGG,CAGTGATCCGCT
        50c59299-7df7-4b7b-a1b4-b96b048c7048,TRB,TRBV5-8*01,,TRBJ2-7*01,,,CGGTTACGACCATCCC,CTCATCGGCCCA
        92082556-e90f-490c-9d04-fcef48c67412,TRB,TRBV6-5*01,,TRBJ2-4*01,,,GTAGCCCTGTACTTAG,GGCCACTCTGTA
        5b0fc427-247c-4460-a7c8-b1640207d3f9,TRB,TRBV7-3*01,,TRBJ2-1*01,GCCAGCAGCTTAATGGCAGAAACAATGAGAGCAGTTC,PAA*WQKQ*EQF,GACGAACGGTAGATCC,ATCTAGCGCGGT
        ...
        
        
        ./Outfolder/Demultiplexing_sample/sample_demux_barcode_umi.csv
        Holds a Table with the spatial Barcode for a Read, its ReadID as well as extracted UMI.

        *Spatial Barcode,ReadID,UMI*
        AAACAAGTATCTCCCA,2ce0b730-1620-44ec-aa34-baf3d73409a9,TCAGAGTAATAG
        AAACAAGTATCTCCCA,a4d0c4bf-87a2-488a-bcf6-12dfcc94aef2,ATGCTGCTATAG
        AAACAAGTATCTCCCA,0f53b8d7-cfcf-440f-a164-7f16627a40ba,ACGGTGCATATG
        AAACAATCTACTAGCA,12a1e67b-6fc7-4c8a-ab35-d8656f141d04,ACGACCGAACAC
        AAACAATCTACTAGCA,7e3656a5-8ca5-4daa-8f80-92e2ca478384,GTCAGTTTCAGT
        AAACAATCTACTAGCA,8f04ac89-d118-4907-bcf7-2b81eeea2850,ACCTCCGGTTCC
        AAACAATCTACTAGCA,dcc8192d-57e4-4da9-95d7-0fec363de285,AGTCCAGGATTG
        AAACAATCTACTAGCA,3175222b-7325-4b41-9398-8cabb79e6d7b,TGAAACTGCTTT
        AAACAATCTACTAGCA,0a620e60-0375-4dc1-8126-146c6688ebb0,AACCGTTGGTTC
        AAACAATCTACTAGCA,ffbe151e-e3c5-4d7d-94e6-8c36bcb25676,ATATCTAAGAGC
        AAACAATCTACTAGCA,c35e7676-5166-4fcf-b604-65d83cef5f27,ATTTTTGTAGTG


        ./Outfolder/Demultiplexing_sample/sample_all_barcode_matches.csv
        Table with all Barcode Matches per ReadID.

        *ReadID,Barcode,0,1,2*
        1e3bc9e7-6e55-4242-806d-d6eae46b621f,GGTGGGATTAGGTCCC,GGTGGGATTAGGTCCC,,
        f8e40593-e913-4259-9a42-9266a7c14f09,TCCTTACGACGGTCCG,TCCTTACGACGGTCCG,,
        1c3e69d7-8858-43df-a329-f3f551fbf1e4,AACCAGACTCTTCGGT,AACCAGACTCTTCGGT,,
        b2e47859-226a-4446-9c86-3c1999fc90f7,AGTGATGACGCGTATA,AGTGATGACGCGTATA,,
        bba4f060-27ae-4eba-9566-cdc7a92a4d17,CGTAGCAGTAATGGAC,CGTAGCAGTAATGGAC,,

### 1.1 Demultiplexing Reads
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
        -igb, --INPUT_IGB       If Input IGB is given, it will be demultiplexed and the UMI added.

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



### 1.2 Preprocess Reads
***
        Preprocess the Reads for the Correction Pipeline. This script splits fusioned Reads, trims the Adapters from 10X & aligns the Reads for the TCR Annotations. If the Barcode UMI csv from Step 1 will be used to demultiplex the annotated TCR sequences.

        2_Preprocess_Reads.sh [-h] [-n NAME] -i INPUT_FASTQ [-o OUTFOLDER]
                                [-t THREADS] [-mem MEMORY] [-rep REPOSITORY]
                                [-pri PRIMER] [-conf CONFIGURATION]
                                [-chop PYCHOPPER] [-trim ADAPTER_TRIM]
                                [-igb IGBLAST] [-demux DEMULTIPLEX]

        -i INPUT_FASTQ          Path to the Input Fastq
        -o OUTFOLDER            Directory for the Outfolder, default: PWD
        -rep REPOSITORY         Path to the cloned Github Repository, default: ../
        -conf CONFIGURATION     Path to the Primer Configuration File for PyChopper, default: ./Reference/Primer/Pychopper/                     
                                10XPrimers_pychopper_configuration.txt
        -pri PRIMER             Path to the Primer .fa File for PyChopper, default: ./Reference/Primer/Pychopper/                     
                                10XPrimers_pychopper.fa
        -n NAME                 Sample Name, if not specified it will default to the basename of the fastq_%d_%Y
        -a ADAPTER              Sequence of the R1 Adapter used for Library Preparation default: 'CTACACGACGCTCTTCCGATCT'
        -t THREADS              Number of Threads to use
        -chop PYCHOPPER         Specify if Reads should be made full length by Pychopper, default: True
        -igb IGBLAST            If True, the preprocessed Fastq is aligned with IgBLAST Following Processing., default: True
        -trim ADAPTER_TRIM      Specify if Reads should be from Adapters
        -demux DEMULTIPLEX      If set to True, extracts Barcode and UMI Region of the Reads and updates the IgBlast Table. Form is default for downstream purposes. default=True


        #### Examplary Output:
        Outputs a modified fastq that was reoriented and adapters trimmed. If you want to use it be aware of the modified readnames.

        *./Outfolder/sample_Cutadapt_trimmed.fastq*
                @24:817|b2e47859-226a-4446-9c86-3c1999fc90f7 runid=bad19e6433e2af91b3d16445abff934b109b2e55 sampleid=SPTCR12 read=183614 ch=63 start_time=2021-10-02T01:08:59Z strand=+
                TTTTTTTTTTTTGAAGTGGTTGTGCGTTCTTTTGTGTGATCAAAACTTCACACAATTGGAAAATAAATGTTTCTTCGAAAATAGAATAATCAAACAAAATTATCCAGGACCTTATAGGGTTTTCAGTATGTACCAAGAGCGTACATCTTAGAAGACCAGGACCTTGTTATCACTGGGATCATTAGGTGGCTTTGAATTGTTTCTTGGGTAGCTTCGTCAGCTTCTCCTTAAACTTGTCAAAGGAACCAAAGTCACGTTGATGGCTTCGTAACTCCCATCTGGAAGTTCTCCACCACCGTTAGGGCTAGGTTGTAAATGCTAAGTTGATACCGTACCACCATTCGGAACTTCGTTGGAGGCTGAAGCTATCTGGGCGTGTAACATCTCCCTTGGCCAACGCCTCCTGGTACTTCTCCTCGGTGACGTTCAGGGAGTTGTGTTCACGTAGGCCGCATTGGTCTTGCTGGTCGTACATGATCTGCGCGTTAGTGTAGAGTTCCAGGGCGCCGTAGTCGTAGGGCAGGTCGGGGGAGGCTGTTGCTTCTGCCTGGTAGGCCCCAGATACCCCAAAGCCGAACCGGCTACCTGCAAGTGCCGCACGCTGCCAGCTCAACATGCTGCTAGTGCTAGTGCCACCGCTGATGCTGTAGTCTGCTGAAGCCGCTGCCCGACCACCGATACGAGTC
                +
                889::;<;84/))'&&'())/.'%%%%%&&-,('('''''())+))**&&++.01126522=;?ABA?@BBC;;;3&&&&&&&&&&&()+/3311//--+))*1.,,.11ABB:9(('''55.2*)))*.221.(''()++(%&**-9:100056665/+''$###$%%&*&&'&&''(10/.+('$%%%%&&&&&'))%%$%$%%(&&''()*--,,++,-//495)))/9<CB@A@><101.,('&''),0((((*-,6(((():&%%$$(,-2''%%%%%$$$'-65/(''''()04-)))-./..()'&&%%%%'$$$$&*++-+&%%$$%%&+0)((&&&'*46.*(&&&&%%&&%%$$$&&(*/10.--((&%&))8=>?95'')4==;821111642--32,,,*)*-.59?AABBBA@@A0***),..&&&&&.,,,56,,47//*)((()+(((**))+*('&'''(235423.+)*+%&&&&&''%$$&%&')22@B?>>44444700/--)((+-('')*)''''((%$$%%%&&&%$$$#####'((*+-0344200((()()'&%%$$&%$$$&''%$%%%'63311(('&&'&%&%&'((%%%((''((',++,***+,(''(+2>?>>.-('&%%$%%)888999;8))'&'&%%&&''*,&$$%&(*,/1
                @36:263|bba4f060-27ae-4eba-9566-cdc7a92a4d17 runid=bad19e6433e2af91b3d16445abff934b109b2e55 sampleid=SPTCR12 read=99091 ch=121 start_time=2021-10-02T01:11:00Z strand=-
                TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGTGTGATTGAGGGTGGAGTAGATTAGGCGCAGGGGTAGAAGTAGAGGTTAAGGAGGTGATGGCTATGATG
                +
                >>>>>>>>>>>>=====>>??????@AAABAACCCCBBB@>542346677;;87-++-(*.,,**+966//((((((()8779=><::==<81000010.../-&&&&'55:20.-
                @487:1018|92c04cea-5811-4b3a-87f8-0448f1ddd126 runid=bad19e6433e2af91b3d16445abff934b109b2e55 sampleid=SPTCR12 read=87480 ch=147 start_time=2021-10-02T01:10:21Z strand=+
                TGTCACCCAGGCTGGAGTGCAGTGGCATGATCTCAGGCTCACTGCAACCTCCGCTTCCCGGGTTCAAGCAATTCTCCTTCCTCAGCCTCCCAAGTAACTAGGATTACAGGCACGAGCCCCCACGCCCAGCTAATCTTTGTATTTTTAGTAGAGATGGGGTTTCGCCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCCCAGGTGATCTGCCTGCCTCGGCCTCCCAAAGTGCTGGGATTAGAGGCATGAGCCACCACGCCCGGCCAGTTTTCTGAGTTACTTACTGTCACTTAAGAGTTGTCCTCTTACCTTTTAGCCTGCTTTTTTTTTTTTTTCCATGCTGCCCTCCAGAGAGACGTCCTTTTTTTTTTTTTTTGAGATAGGTCAAACAGATTTTAATGTGAGGATAGTATGAAAGTTCAA
                +
                )))).0104=????>>@?<>A??@@@=:999<A<,,++,788535512333621113760////8:<<:/)((')*)))*(((4222230+*----)((&&(((4111/.-(%%%&(),2331**+,-@?>B?C2111*2<=@?CDMEBAEDDF=:8899@D>55555@DEDA?@>>>=??>>==<:9::;@<:::;60////0:@@AA>?=;;<;<>55546B?=<100019=<0.*(()*556688889>>AA=32222;:;;8.---.07.+*('''()+,/.(''&&%&&%##()++*++*012'&&&)(''&&%%%&&&,147899:985.,''&%&&%&11085,,+)(,))&'('(46<9=@@AAA>8+'&&%%&&&&&'''%$%&')11..,++++./))))***1/1.-./01,(

        *./Outfolder/sample_preprocessed_IGB.tsv*
        Holds the output table generated by IGBlast Alignment.

                sequence_id	sequence	locus	stop_codon	vj_in_frame	v_frameshift	productive	rev_comp	complete_vdj	v_call	d_call	j_call	sequence_alignment	germline_alignment	sequence_alignment_aa	germline_alignment_aa	v_alignment_start	v_alignment_end	d_alignment_start	d_alignment_end	j_alignment_start	j_alignment_end	v_sequence_alignment	v_sequence_alignment_aa	v_germline_alignment	v_germline_alignment_aa	d_sequence_alignment	d_sequence_alignment_aa	d_germline_alignment	d_germline_alignment_aa	j_sequence_alignment	j_sequence_alignment_aa	j_germline_alignment	j_germline_alignment_aa	fwr1	fwr1_aa	cdr1	cdr1_aa	fwr2	fwr2_aa	cdr2	cdr2_aa	fwr3	fwr3_aa	fwr4	fwr4_aa	cdr3	cdr3_aa	junction	junction_length	junction_aa	junction_aa_length	v_score	d_score	j_score	v_cigar	d_cigar	j_cigar	v_support	d_support	j_support	v_identity	d_identity	j_identity	v_sequence_start	v_sequence_end	v_germline_start	v_germline_end	d_sequence_start	d_sequence_end	d_germline_start	d_germline_end	j_sequence_start	j_sequence_end	j_germline_start	j_germline_end	fwr1_start	fwr1_end	cdr1_start	cdr1_end	fwr2_start	fwr2_end	cdr2_start	cdr2_end	fwr3_start	fwr3_end	fwr4_start	fwr4_end	cdr3_start	cdr3_end	np1	np1_length	np2	np2_length	v_family	d_family	j_family	cdr3_aa_length
                41:559|d9deb513-ed7a-4494-bb16-df2ab5adb397runid=bad19e6433e2af91b3d16445abff934b109b2e55sampleid=SPTCR12read=109140ch=417start_time=2021-10-02T01:07:29Zstrand=+	TGGGGCCAAACCCGTCACCCAGATCGTCAGCGCCGAGGCCCGGGGTAGAGCAGACTGTGGTTTTACCTCGGTGTCCTACCAGCAAGGGGTCCTGTCTGCCACCATCCTCTATGAGATCCTGCTAGGCTGCCACCCTGCATGTCAGCTGGCCAGCGCCCTTGTGTTGATGGCCATGGCCGCAGAAAGGATTTCTGAAGGCAGCCCTGAAGTGGAGTTAGGAGCTCTAACCCGTCATGGTTCTACACACATTCTTCTTTTGCCAGCGCTTCTGAAGAGCTGCTCTCACCTCTCTGCATCCCAACAGATATCCCCCCATGTGCATGCACACCTGCACACTCACGGCCGAAATCTCCCTAACCCAGGGGACCTTAGCATGCCTAAGTGACTAAACCAATAAAAATGAAACT	TRA	F				T	F	TRBV5-1*01			ACACATTCTTCTTTTGCCAGCGCTTCTGAAGAGCTGCTCTCACCTCTCTG	GATCAAAACGAGAGGACAGCAAGTGACACTGAGCTGCTCCCCTATCTCTG	HILLLPALLKSCSHLS	IKTRGQQVTLSCSPIS	1	50					ACACATTCTTCTTTTGCCAGCGCTTCTGAAGAGCTGCTCTCACCTCTCTG	HILLLPALLKSCSHLS	GATCAAAACGAGAGGACAGCAAGTGACACTGAGCTGCTCCCCTATCTCTG	IKTRGQQVTLSCSPIS																											20.534			243S32N50M114S204N			9.489e+00			40.000			244	293	33	82																											TRBV5-1			0



### 3. Cluster and Correct Reads
***
        Clusters Reads based on their VJ-Family Annotation generated by IgBlast Alignment, next Read groups are parallel corrected using Rattle Algorithm and finally annotated by IGBlast.
        
        ./3_Cluster_Correct.sh"
        usage: 3_Cluster_Correct.sh [-h] [-n NAME] -i INPUT_FASTQ -b INPUT_IGB [-o OUTFOLDER]
                                [-t THREADS] [-rep REPOSITORY] [-igb IGBLAST] [-cln CLEANUP]

        -i INPUT_FASTQ          Path to the Input Fastq
        -i INPUT_IGB            Path to the Input IgBlast Result generated by 2_Preprocess_Reads.sh
        -n NAME                 Sample Name, if not specified it will default to the basename of the fastq_%d_%Y
        -o OUTFOLDER            Directory for the Outfolder, default: PWD
        -t THREADS              Number of Threads to use
        -igb IGBLAST            If True, the corrected Fastq is aligned with IgBLAST following correction., default: True
        -rep REPOSITORY         Path to the cloned Github Repository, default: ./
        -cln CLEANUP            If True, created intermediate Files and Folders will be deleted. For Debugging you can set 
                                this to False., default: True

        Following Correction you will find the corrected TCR Reads .fastq in the defined Outfolder:

                OUTFOLDER/YOUR_SAMPLE_NAME_corrected_merged.fastq

        The Final Fully Annotated IgBlast File can be found in:

                OUTFOLDER/YOUR_SAMPLE_NAME_corrected_IGB.tsv
        
        For a VDj-Annotation Summary:

                OUTFOLDER/YOUR_SAMPLE_NAME_corrected_IGB_summary.csv
                
                ReadID,Locus,V,D,J,CDR3,CDR3_aa,Spatial Barcode,UMI
                7017f836-e631-407b-a236-e0e9f2f4fef3,TRB,TRBV14*01,TRBD2*01,TRBJ1-6*01,GGGGTTCTGCCAGAAGGTGGCCGAGACCCTCAGGCGGCTGCTCAGGCAGTATCTGGAGTCATTGAGGGCGGGCTGCTCCTTGAGGGGCTGCGGGTCTGTGCTGACCCCACTGTGCACCTCCTTCCCATTCACCC,GVLPEGGRDPQAAAQAVSGVIEGGLLLEGLRVCADPTVHLLPIHP,ACACCTGACACTAGGG,AAGATCAGTTAA
                311a0300-de74-41e9-83d8-91841700a79f,TRB,TRBV14*01,TRBD2*01,TRBJ1-6*01,GGGGTTCTGCCAGAAGGTGGCCGAGACCCTCAGGCGGCTGCTCAGGCAGTATCTGGAGTCATTGAGGGCGGGCTGCTCCTTGAGGGGCTGCGGGTCTGTGCTGACCCCACTGTGCACCTCCTTCCCATTCACCC,GVLPEGGRDPQAAAQAVSGVIEGGLLLEGLRVCADPTVHLLPIHP,ATTGATGAGGAGCGCC,ATGACACTGAGT
                f1ab263f-76ce-48a6-9ec0-51b7b29cef78,TRB,TRBV14*01,TRBD2*01,TRBJ1-6*01,


### Exemplary Pipeline
***
        An Exemplary minimal Usage Pipeline to demultiplex, preprocess and correct an Input Fastq.

        #! /usr/bin/env bash
        
        ## Command into Github Repository

        preproc=./2_Preprocess_Reads.sh
        cluscorr=./3_Cluster_Correct.sh
        input=PATH/TO/INPUT/FASTQ
        MEMORY=32
        THREADS=12
        SAMPLE_NAME=SAMPLE_XY
        OUT=PATH/TO/OUTFOLDER

        ##############################

        bash "${preproc}" \
                -n ${SAMPLE_NAME} \
                -i "${input}" \
                -t ${THREADS} \
                -mem ${MEMORY} \
                -o ${OUT} \

        bash "${cluscorr}" \
                -i "${input}" \
                -b "${OUT}/PreProcessing/${SAMPLE_NAME}_preprocessed_IGB.tsv" \
                -n ${SAMPLE_NAME} \
                -t ${THREADS} \



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