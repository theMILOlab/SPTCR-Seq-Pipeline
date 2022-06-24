#!/usr/bin/env python
import argparse
from email.policy import default

parser = argparse.ArgumentParser()
parser.add_argument('-I','--INPUT_FASTQ',help="Path to Input Fastq File", required=True )
parser.add_argument('-B','--BED',help="Path to Input Bed", required=False, default="" )
parser.add_argument('-O','--OUT',help= "Path to Outfolder", required=True )
parser.add_argument('-L','--LOCI',help="Split by Loci",default=False,type=bool)
parser.add_argument('-V','--VDJ',help="Split by VDJ Arrangement",default=False,type=bool)
parser.add_argument('-igb','--IGB',help="Path to IgBlast Tsv",default="")
parser.add_argument('-VW','--VDJ_WRITE',help="Read in Arrangements & from csv & split FASTQ",default=False,type=bool)
parser.add_argument('-VWT','--VDJ_WRITE_CSVPATH',help="Path to CSVs to split FASTQ into Arrangements")

#parser.add_argument('-BAM','--BAM',help="Path to Input BAM", required=True )

args = parser.parse_args()

arg_vars = vars(args)

############### Import Python Module ###############
import pyfastx
import pandas as pd
import dask.dataframe as dd
from joblib import Parallel, delayed
from collections import defaultdict
import os  
import glob
from datetime import datetime

#######################################################

OUT=str(arg_vars["OUT"])
LOCI=arg_vars["LOCI"]
VDJ=arg_vars["VDJ"]
VDJ_WRITE=arg_vars["VDJ_WRITE"]
ARRANGEMENT_CSVs=arg_vars["VDJ_WRITE_CSVPATH"]

############# Read in Fastq ############# 
FASTQ=pyfastx.Fastq(arg_vars["INPUT_FASTQ"])

############# Read in Alignment BED  #############
if arg_vars["BED"] not in "":
    BED_ALN=dd.read_csv(arg_vars["BED"],sep="\t",low_memory=False,names=["Aligned Transcript",'Start_a', 'End_a','ReadID','Score','Strand'])
    BED_ALN=BED_ALN.compute()

#VJC_ALN=pysam.AlignmentFile(str(arg_vars["BAM"]),'rb')

############# Splitting by IgBlast VDJ Arrangement ############# 

if arg_vars["IGB"] not in "":
    print("Reading in IgBlast Result")
    ### Read in IgBlast Result
    
    IGB=dd.read_csv(arg_vars["IGB"],sep="\t",low_memory=True,usecols=['sequence_id','locus','cdr3', 'v_call','d_call','j_call'], dtype={'fwr4': 'object','fwr4_aa': 'object','cdr3': 'object','v_alignment_end': 'float64','v_alignment_start': 'float64','v_germline_end': 'float64','v_germline_start': 'float64','v_sequence_end': 'float64','v_sequence_start': 'float64'},lineterminator='\n').set_index('sequence_id')

    IGB=IGB.rename(columns={'sequence':'locus'})
    IGB=IGB.dropna(subset=["locus"])
    #,axis=0
    cdr3_assigned=IGB[(~IGB["v_call"].isna()) & (~IGB["j_call"].isna())]
    #IGB=IGB.compute()
    
    
    
    ### Ignoring Exon
    
    #(~SPTCR16["cdr3"].isna()) & 
    #cdr3_assigned=cdr3_assigned[cdr3_assigned['cdr3'].str.len() > 5]

    cdr3_assigned["v_call Exon"]=cdr3_assigned['v_call'].str.split(pat="*",expand=True, n=1)[1]
    cdr3_assigned["v_call"]=cdr3_assigned['v_call'].str.split(pat="*",expand=True, n=1)[0]

    cdr3_assigned["j_call Exon"]=cdr3_assigned['j_call'].str.split(pat="*",expand=True, n=1)[1]
    cdr3_assigned["j_call"]=cdr3_assigned['j_call'].str.split(pat="*",expand=True, n=1)[0]

    ##### Split Readname to deal with modified Readname
    cdr3_assigned=cdr3_assigned.reset_index(drop=False)
    
    print(cdr3_assigned['sequence_id'].compute())
    
    index_cols=cdr3_assigned["sequence_id"].str.split(pat="runid=",expand=True,n=1)
    cdr3_assigned["sequence_id"]=index_cols[0]
    cdr3_assigned["sequence_id"]=cdr3_assigned["sequence_id"].str.split(pat="|",expand=False,n=1).str[1]
    #index_cols=cdr3_assigned["index"].str.split(pat="|",expand=True,n=1)
    #cdr3_assigned["index"]=index_cols[1]
    cdr3_assigned=cdr3_assigned.set_index("sequence_id")
    #cdr3_assigned=cdr3_assigned.compute()
    print(cdr3_assigned.head(10))

    

############# Splitting by Locus in Lists Block #############

            
if LOCI == True:
    now = datetime.now().time()
    def SPLIT_BY_LOCI(GROUPS, LOCUS, DICTIONARY):
        for READGROUP, read_aln in GROUPS:
            for row in read_aln["Aligned Transcript"]:
                if LOCUS in row.split("|")[1]:
                    #for readid in read_aln["ReadID"]:
                    DICTIONARY[LOCUS].add(READGROUP)
        txt_name="{0}.txt".format(LOCUS)
        with open(txt_name,"a") as f:
            for readid in DICTIONARY[LOCUS]:
                f.write(readid)
                f.write('\n')
                
        return DICTIONARY
    
    print("Splitting Input Fastq by Alignment in Constant Loci",now)
    loci = defaultdict(set)
    grouper=BED_ALN.groupby("ReadID")
    
    #### Parallel Process TCR Loci
    
    Parallel(verbose=5,n_jobs=-1,backend='multiprocessing')(delayed(SPLIT_BY_LOCI)(GROUPS=grouper, LOCUS=Locus, DICTIONARY=loci) for Locus in ["TRA","TRB","TRD","TRG"])  
    
    ### read Files into dictionary
    
    for Locus in ["TRA","TRB","TRD","TRG"]:
        filename="{0}.txt".format(Locus)
        with open(filename) as f:
            for line in f:
                line=line.strip()
                loci[Locus].add(line) 

############# Splitting by Minimap VDJ Arrangement Block #############
elif VDJ == True:
    print(arg_vars["BED"])
    if "VDJ" in os.path.basename(arg_vars["BED"]):
        print("Writing VDJ Clonal Arrangement DF")
        only_full_len=BED_ALN[(BED_ALN["End_a"]-BED_ALN["Start_a"]>200)]
        grouper=only_full_len.groupby("ReadID")
        #VDJ_arrangements=pd.DataFrame(columns=["V","D","J","C"])
        tuple_list=[]
        for (readname, Exons) in grouper:
            for exon in Exons["Aligned Transcript"]:
                constructed_tuple=(readname,str(exon).split("|")[1],str(exon).split("|+++|")[2].split("|")[1],str(exon).split("|+++|")[1].split("|")[1],str(os.path.basename(arg_vars["BED"])).split("_")[0])
                tuple_list.append(constructed_tuple)
                #VDJ_arrangements.loc[readname,"V"]=str(exon).split("|")[1]
                #VDJ_arrangements.loc[readname,"J"]=str(exon).split("|+++|")[2].split("|")[1]
                #VDJ_arrangements.loc[readname,"D"]=str(exon).split("|+++|")[1].split("|")[1]
                #VDJ_arrangements.loc[readname,"C"]=str(os.path.basename(arg_vars["BED"])).split("_")[0]
        VDJ_arrangements=pd.DataFrame.from_records(tuple_list, columns =["ReadID","V","D","J","C"])
        name=OUT+"/"+str(os.path.basename(arg_vars["BED"])).split("_")[0]+"_full_len_VDJ_arrangements.csv"
        VDJ_arrangements.to_csv(name)
        print(VDJ_arrangements)
        
    if "VJ" in os.path.basename(arg_vars["BED"]):
        print("Writing VJ Clonal Arrangement DF")
        only_full_len=BED_ALN[(BED_ALN["End_a"]-BED_ALN["Start_a"]>200)]
        grouper=only_full_len.groupby("ReadID")
        #VJ_arrangements=pd.DataFrame(columns=["V","J","C"])
        tuple_list=[]
        for (readname, Exons) in grouper:
            for exon in Exons["Aligned Transcript"]:
                constructed_tuple=(readname,str(exon).split("|")[1],str(exon).split("|+++|")[1].split("|")[1],str(os.path.basename(arg_vars["BED"])).split("_")[0])
                tuple_list.append(constructed_tuple)
                #VJ_arrangements.loc[readname,"V"]=str(exon).split("|")[1]
                #VJ_arrangements.loc[readname,"J"]=str(exon).split("|+++|")[1].split("|")[1]
                #VJ_arrangements.loc[readname,"C"]=str(os.path.basename(arg_vars["BED"])).split("_")[0]
        VJ_arrangements=pd.DataFrame.from_records(tuple_list, columns =["ReadID","V","J","C"])
        name=OUT+"/"+str(os.path.basename(arg_vars["BED"])).split("_")[0]+"_full_len_VJ_arrangements.csv"
        VJ_arrangements.to_csv(name)
        print(VJ_arrangements)

if VDJ_WRITE == True:
    arrangements = defaultdict(dict)
    pattern=str(ARRANGEMENT_CSVs)+"/*.csv"
    for file in glob.iglob(pattern, recursive=False):
        vdj_pattern=pd.read_csv(file,index_col=0).groupby(["V","C"])
        for (group, reads) in vdj_pattern: 
            arrangements[group[1]][group[0]] = set(reads.index)
            #arrangements[group[1]] = defaultdict(dict)
            #arrangements[group[1]][group[0]] = defaultdict(set)
            #arrangements[group[1]][group[0]].add(set(reads.index))
    
############# Define parallelized Writing Functions #############  

def write_Cluster_Fastq_from_list(READ_LIST,INPUT_FASTQ,LOCUS,V_SEG ,OUTPATH):
    #print(LOCUS)
    V_SEG=V_SEG.replace("/","_")
    CLUSTER=LOCUS+"_"+V_SEG
    
    path=os.path.join(OUTPATH,CLUSTER)
    #try:
    #    path=os.path.join(OUTPATH,CLUSTER)
    #    os.mkdir(path)
    #except FileExistsError:
    #    pass
    FILENAME=OUTPATH+'/{0}.fastq'.format(CLUSTER)
    #FILENAME=FILENAME.replace("/","_")
    #write_path=os.path.join(path,FILENAME)
    with open(FILENAME,'a') as f:
        for read in READ_LIST:
            f.write(INPUT_FASTQ[read].raw)
    return None

def write_Locus_Fastq_from_list(READ_LIST,INPUT_FASTQ,LOCUS ,OUTPATH):
    #print(LOCUS)
    try:
        path=os.path.join(OUTPATH,LOCUS)
        os.mkdir(path)
    except FileExistsError:
        pass
    FILENAME='{0}.fastq'.format(LOCUS)
    write_path=os.path.join(path,FILENAME)
    with open(write_path,'a') as f:
        for read in READ_LIST:
            f.write(INPUT_FASTQ[read].raw)
    return None


def write_Arrang_Fastq_from_IGB(GROUP,GROUPNAME,INPUT_FASTQ, OUTPATH):
    ### Group By VJ Arrangements
    
    #for (groupname, group) in grouper:
    #GROUP=GROUP.compute()
    LOCUS=GROUPNAME[0]
    name="_".join(GROUPNAME)
    name=name.replace("/","__")
    name=name.replace("*","+")
    READ_LIST=list(GROUP.index)

    #try:
    #    path=os.path.join(OUTPATH,LOCUS)
    #    os.mkdir(path)
    #except FileExistsError:
    #    pass
    #FILENAME='{0}.fastq'.format(groupname)
    #write_path=os.path.join(path,FILENAME)
    FILENAME=OUTPATH+'/{0}.fastq'.format(name)
    with open(FILENAME,'a') as f:
        for read in READ_LIST:
            try:
                f.write(INPUT_FASTQ[read].raw)
            except:
                with open(OUTPATH+'/Failed_reads_.txt','a') as p:
                    p.write(read)
                    p.write("\n")
    return None



#### Execution Block ####

now = datetime.now().time()

if LOCI == True:
    now = datetime.now().time()
    print('Process based parallel writing of Locus', now)
    #prefer="processes",n_jobs=-1,backend='multiprocessing'

    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Locus_Fastq_from_list)(READ_LIST=readlist, INPUT_FASTQ=FASTQ, LOCUS=locus, OUTPATH=OUT) for locus, readlist in loci.items())

    finished=datetime.now().time()
    print('Finished parallel writing of Clusters', finished)
    
if VDJ_WRITE == True:
    now = datetime.now().time()
    print('Process based parallel writing of Clusters', now)

    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Cluster_Fastq_from_list)(READ_LIST=reads, INPUT_FASTQ=FASTQ, LOCUS=locus,V_SEG=v_seg, OUTPATH=OUT) for locus,v_dict in arrangements.items() for v_seg, reads in v_dict.items())

    finished=datetime.now().time()
    print('Finished parallel writing of Clusters', finished)
    

if arg_vars["IGB"] not in "":
    now = datetime.now().time()
    
    TR=cdr3_assigned[cdr3_assigned["locus"]=="TRB"]
    TR=TR.compute()
    print(TR.head())
    grouper=TR.groupby(by=["locus","v_call","j_call"])
    #,"d_call"
    print('Process based parallel writing of TRB VDJ Clusters', now)
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)
    finished=datetime.now().time()
    print('Finished parallel writing of TRB Clusters', finished)
    
    now = datetime.now().time()
    
    TR=cdr3_assigned[cdr3_assigned["locus"]=="TRD"]
    TR=TR.compute()
    print(TR.head())
    grouper=TR.groupby(by=["locus","v_call","j_call"])
    #,"d_call"
    print('Process based parallel writing of TRD VDJ Clusters', now)
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)
    finished=datetime.now().time()
    print('Finished parallel writing of TRD Clusters', finished)
    
    now = datetime.now().time()
    print('Process based parallel writing of VJ Clusters', now)
    
    TR=cdr3_assigned[cdr3_assigned["locus"]=="TRA"]
    TR=TR.compute()
    print(TR.head())
    grouper=TR.groupby(by=["locus","v_call","j_call"])
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)
    finished=datetime.now().time()
    print('Finished parallel writing of TRA Clusters', finished)
    now = datetime.now().time()
    print('Process based parallel writing of VJ Clusters', now)
    
    TR=cdr3_assigned[cdr3_assigned["locus"]=="TRG"]
    TR=TR.compute()
    print(TR.head())
    grouper=TR.groupby(by=["locus","v_call","j_call"])
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)

    finished=datetime.now().time()
    print('Finished parallel writing of TRG Clusters', finished)
 