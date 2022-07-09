#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-I','--INPUT_FASTQ',help="Path to Input Fastq File", required=True )
parser.add_argument('-O','--OUT',help= "Path to Outfolder", required=True )
parser.add_argument('-igb','--IGB',help="Path to IgBlast Tsv",default="")
parser.add_argument('-grp','--GROUPER',help="Grouper to use for grouping the TCR Reads into Fastqs",default=["Locus","V","J"])

args = parser.parse_args()

arg_vars = vars(args)

############### Import Python Module ###############
import pyfastx
import pandas as pd
from joblib import Parallel, delayed
from collections import defaultdict
import os
from datetime import datetime

############### Define Input Variables #####################
GROUPER=str(arg_vars["GROUPER"]).split(',')
print(GROUPER)
#GROUPER=ast.literal_eval(GROUPER)

OUT=str(arg_vars["OUT"])
IGB_PATH=arg_vars["IGB"]

print("Grouping TCRs by:",GROUPER)
############# Read in Fastq ############# 
FASTQ=pyfastx.Fastq(arg_vars["INPUT_FASTQ"])

############# Splitting by IgBlast VJ-Family Arrangement ############# 

if arg_vars["IGB"] not in "":
    print("Reading in IgBlast Result")
    
    ### Read in IgBlast Result
    IGB=pd.read_csv(IGB_PATH,low_memory=False,index_col=0)
    cdr3_assigned=IGB[(~IGB["V"].isna()) & (~IGB["J"].isna())]

    ### Only VJ Family
    cdr3_assigned["V Exon"]=cdr3_assigned['V'].str.split(pat="*",expand=False, n=1).str[1]
    cdr3_assigned["V"]=cdr3_assigned['V'].str.split(pat="*",expand=False, n=1).str[0]

    cdr3_assigned["J Exon"]=cdr3_assigned['J'].str.split(pat="*",expand=False, n=1).str[1]
    cdr3_assigned["J"]=cdr3_assigned['J'].str.split(pat="*",expand=False, n=1).str[0]

    print(cdr3_assigned.head(10))

############# Define parallelized Writing Functions #############  
mod_read_ids=list(FASTQ.keys())
true_read_ids=[str(id.name).split('|')[1] for id in FASTQ]
FASTQ_ROSETTA=dict(zip(true_read_ids,mod_read_ids))
#print(FASTQ_ROSETTA)

def write_Arrang_Fastq_from_IGB(GROUP,GROUPNAME,INPUT_FASTQ, OUTPATH,FASTQ_ROSETTA=FASTQ_ROSETTA):
    ### Group By VJ Arrangements
    LOCUS=GROUPNAME[0]
    name="_".join(GROUPNAME)
    name=name.replace("/","__")
    name=name.replace("*","+")
    
    ## Translate Read IDs modified by PyChopper
    READ_LIST=[]
    READ_LIST_true=list(GROUP.index)
    for read in READ_LIST_true:
        READ_LIST.append(FASTQ_ROSETTA[read])

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
    grouper=TR.groupby(by=GROUPER)
    print('Process based parallel writing of TRB VDJ Clusters', now)
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)
    finished=datetime.now().time()
    print('Finished parallel writing of TRB Clusters', finished)
    
    now = datetime.now().time()
    
    TR=cdr3_assigned[cdr3_assigned["Locus"]=="TRD"]
    #print(TR.head())
    grouper=TR.groupby(by=GROUPER)
    print('Process based parallel writing of TRD VDJ Clusters', now)
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)
    finished=datetime.now().time()
    print('Finished parallel writing of TRD Clusters', finished)
    
    now = datetime.now().time()
    
    print('Process based parallel writing of VJ Clusters', now)
    
    TR=cdr3_assigned[cdr3_assigned["Locus"]=="TRA"]
    #print(TR.head())
    grouper=TR.groupby(by=GROUPER)
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)
    finished=datetime.now().time()
    print('Finished parallel writing of TRA Clusters', finished)
    now = datetime.now().time()
    
    print('Process based parallel writing of VJ Clusters', now)
    
    TR=cdr3_assigned[cdr3_assigned["Locus"]=="TRG"]
    #print(TR.head())
    grouper=TR.groupby(by=GROUPER)
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)

    finished=datetime.now().time()
    print('Finished parallel writing of TRG Clusters', finished)