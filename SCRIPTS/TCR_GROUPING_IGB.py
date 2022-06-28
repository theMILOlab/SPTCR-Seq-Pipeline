#!/usr/bin/env python
import argparse
from email.policy import default

parser = argparse.ArgumentParser()
parser.add_argument('-I','--INPUT_FASTQ',help="Path to Input Fastq File", required=True )
parser.add_argument('-O','--OUT',help= "Path to Outfolder", required=True )
parser.add_argument('-igb','--IGB',help="Path to IgBlast Tsv",default="")

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
IGB_PATH=arg_vars["IGB"]

############# Read in Fastq ############# 
FASTQ=pyfastx.Fastq(arg_vars["INPUT_FASTQ"])


############# Splitting by IgBlast VDJ Arrangement ############# 

if arg_vars["IGB"] not in "":
    print("Reading in IgBlast Result")
    ### Read in IgBlast Result
    
    IGB=pd.read_csv(IGB_PATH,low_memory=False,index_col=0)
    cdr3_assigned=IGB[(~IGB["V"].isna()) & (~IGB["J"].isna())]

    ### Ignoring Exon
    cdr3_assigned["V Exon"]=cdr3_assigned['V'].str.split(pat="*",expand=False, n=1).str[1]
    cdr3_assigned["V"]=cdr3_assigned['V'].str.split(pat="*",expand=False, n=1).str[0]

    cdr3_assigned["J Exon"]=cdr3_assigned['J'].str.split(pat="*",expand=False, n=1).str[1]
    cdr3_assigned["J"]=cdr3_assigned['J'].str.split(pat="*",expand=False, n=1).str[0]

    print(cdr3_assigned.head(10))


    
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


if arg_vars["IGB"] not in "":
    now = datetime.now().time()
    
    TR=cdr3_assigned[cdr3_assigned["Locus"]=="TRB"]
    #print(TR.head())
    grouper=TR.groupby(by=["Locus","V","J"])
    print('Process based parallel writing of TRB VDJ Clusters', now)
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)
    finished=datetime.now().time()
    print('Finished parallel writing of TRB Clusters', finished)
    
    now = datetime.now().time()
    
    TR=cdr3_assigned[cdr3_assigned["Locus"]=="TRD"]
    #print(TR.head())
    grouper=TR.groupby(by=["Locus","V","J"])
    print('Process based parallel writing of TRD VDJ Clusters', now)
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)
    finished=datetime.now().time()
    print('Finished parallel writing of TRD Clusters', finished)
    
    now = datetime.now().time()
    
    print('Process based parallel writing of VJ Clusters', now)
    
    TR=cdr3_assigned[cdr3_assigned["Locus"]=="TRA"]
    #print(TR.head())
    grouper=TR.groupby(by=["Locus","V","J"])
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)
    finished=datetime.now().time()
    print('Finished parallel writing of TRA Clusters', finished)
    now = datetime.now().time()
    
    print('Process based parallel writing of VJ Clusters', now)
    
    TR=cdr3_assigned[cdr3_assigned["Locus"]=="TRG"]
    #print(TR.head())
    grouper=TR.groupby(by=["Locus","V","J"])
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)

    finished=datetime.now().time()
    print('Finished parallel writing of TRG Clusters', finished)