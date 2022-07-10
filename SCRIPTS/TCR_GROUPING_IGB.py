#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-I','--INPUT_FASTQ',help="Path to Input Fastq File", required=True )
parser.add_argument('-O','--OUT',help= "Path to Outfolder", required=True )
parser.add_argument('-igb','--IGB',help="Path to IgBlast Tsv. Use Full IgBlast Output from PreProcessed Reads",default="")
parser.add_argument('-grp','--GROUPER',help="Grouper to use for grouping the TCR Reads into Fastqs. Pick one or a combination of: locus,v_family, d_family, j_family",default=["locus,v_family"])
parser.add_argument('-ids','--READID_COL',help="Column for ReadIds",default="sequence_id")
parser.add_argument('-loc','--LOCUS_COL',help="Column for TCR Loci",default="locus")
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
GROUPER=[str(groupie).strip() for groupie in GROUPER]
READID_COL=str(arg_vars["READID_COL"])
LOCUS_COL=str(arg_vars["LOCUS_COL"])

OUT=str(arg_vars["OUT"])
IGB_PATH=arg_vars["IGB"]

print("Grouping TCRs by:",GROUPER, "Using ReadIds from Column:",READID_COL)

############# Read in Fastq ############# 
FASTQ=pyfastx.Fastq(arg_vars["INPUT_FASTQ"])

############# Splitting by IgBlast VJ-Family Arrangement ############# 
print("Reading in IgBlast Result and preparing it")

### Read in IgBlast Result
if IGB_PATH.endswith('.tsv'):
    ending='\t'
elif IGB_PATH.endswith('.csv'):
    ending=','

########
cols=GROUPER.copy()
cols.append(READID_COL)
print('Only Reading:',cols)
IGB=pd.read_csv(IGB_PATH,sep=ending,low_memory=False,usecols=cols,index_col=READID_COL)

## Modify Read ID
#IGB[READID_COL]=IGB[READID_COL].str.split('|').str[1]
IGB.index=IGB.index.str.split('runid').str[0]

### Drop Rows with Na in Grouper
cdr3_assigned=IGB[~IGB[GROUPER[0]].isna()]
for count,grouper in enumerate(GROUPER[1:]):
    count+=1
    cdr3_assigned=cdr3_assigned[~cdr3_assigned[GROUPER[count]].isna()]
print(cdr3_assigned)

### Only VJ Family
#cdr3_assigned["V Exon"]=cdr3_assigned['V'].str.split(pat="*",expand=False, n=1).str[1]
#cdr3_assigned["V"]=cdr3_assigned['V'].str.split(pat="*",expand=False, n=1).str[0]

#cdr3_assigned["J Exon"]=cdr3_assigned['J'].str.split(pat="*",expand=False, n=1).str[1]
#cdr3_assigned["J"]=cdr3_assigned['J'].str.split(pat="*",expand=False, n=1).str[0]

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
    #for read in READ_LIST_true:
    #    READ_LIST.append(FASTQ_ROSETTA[read])

    FILENAME=OUTPATH+'/{0}.fastq'.format(name)
    
    with open(FILENAME,'a') as f:
        for read in READ_LIST_true:
            #READ_LIST_true
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
    
    TR=cdr3_assigned[cdr3_assigned[LOCUS_COL]=="TRB"]
    #print(TR.head())
    grouper=TR.groupby(by=GROUPER)
    print('Process based parallel writing of TRB VDJ Clusters', now)
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)
    finished=datetime.now().time()
    print('Finished parallel writing of TRB Clusters', finished)
    
    now = datetime.now().time()
    
    TR=cdr3_assigned[cdr3_assigned[LOCUS_COL]=="TRD"]
    #print(TR.head())
    grouper=TR.groupby(by=GROUPER)
    print('Process based parallel writing of TRD VDJ Clusters', now)
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)
    finished=datetime.now().time()
    print('Finished parallel writing of TRD Clusters', finished)
    
    now = datetime.now().time()
    
    print('Process based parallel writing of VJ Clusters', now)
    
    TR=cdr3_assigned[cdr3_assigned[LOCUS_COL]=="TRA"]
    #print(TR.head())
    grouper=TR.groupby(by=GROUPER)
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)
    finished=datetime.now().time()
    print('Finished parallel writing of TRA Clusters', finished)
    now = datetime.now().time()
    
    print('Process based parallel writing of VJ Clusters', now)
    
    TR=cdr3_assigned[cdr3_assigned[LOCUS_COL]=="TRG"]
    #print(TR.head())
    grouper=TR.groupby(by=GROUPER)
    #dropna=False,
    Parallel(verbose=5,backend='multiprocessing')(delayed(write_Arrang_Fastq_from_IGB)(GROUP=group, GROUPNAME=groupname ,INPUT_FASTQ=FASTQ, OUTPATH=OUT) for (groupname, group) in grouper)

    finished=datetime.now().time()
    print('Finished parallel writing of TRG Clusters', finished)