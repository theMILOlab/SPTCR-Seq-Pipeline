#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-O','--OUT',help= "Path to Outfolder", required=True )
parser.add_argument('-igb','--IGB',help="Path to IgBlast csv output from SPTCR Cluster Correct Pipeline",default="")
parser.add_argument('-n','--NAME',help="Name of the Pipeline", required=True )
parser.add_argument('-u','--UNCORR',help="Path to uncorrected IGB", required=True )

##################### Import Modules ############################
import modin.pandas as pmd
from umi_tools import UMIClusterer
from collections import defaultdict
import sys
import os
import json
import pandas as pd
from itertools import chain
import numpy as np
import ast
import sklearn.preprocessing as sk
#######################################################
#######################################################

OUT=str(arg_vars["OUT"])
sample_name=arg_vars["NAME"]
read_dir=arg_vars["IGB"]
uncorr=arg_vars["UNCORR"]


#######################################################
################ Umi Correction Pipeline ##############
#######################################################
print('::::::::',sample_name,'::::::::')
print(read_dir)
sample_df=pmd.read_csv(read_dir)

##############################################
#sample_df=sample_df.sample(n=100)
##############################################

barcodes=sample_df['Visium Barcode'].unique()

####
clusterer = UMIClusterer(cluster_method="directional")
all_clustered_umis=defaultdict(list)
distance=2
####

for count, barcode in enumerate(barcodes):
    #if count < 10:
    print('Clustering UMIs:',sample_name,count)
    selected_BC=sample_df[sample_df['Visium Barcode']==barcode]
    umi_col=selected_BC['UMI'].str.encode(encoding='utf-8')
    umis=umi_col.groupby('UMI').size().to_dict()
    clustered_umis = clusterer(umis, threshold=distance)
    all_clustered_umis[barcode].append(clustered_umis)
    
construc_dict={}
for barcode in all_clustered_umis.keys():
    a=list(chain.from_iterable(all_clustered_umis[barcode]))
    umi_series=pd.Series(a)
    umi_series=umi_series.apply(lambda x: [thing.decode("utf-8") for thing in x])
    construc_dict[barcode]=umi_series
umis_per_BC=pd.DataFrame.from_dict(construc_dict)
print(umis_per_BC.head(5))
umis_per_BC.to_csv('{0}}/UMI_CLUSTERING/{1}_clustered_umis_edit_2.csv'.format(OUT,sample_name),index=False)

#### Modify UMI Counts

### Read Clustered UMIs as DF and modify to Default Dict
clustered_umis=pd.read_csv('{0}/UMI_CLUSTERING/{1}_clustered_umis_edit_2.csv'.format(OUT,sample_name),low_memory=False)

clustered_umis_dict=clustered_umis.to_dict(orient='list')

for key, value in clustered_umis_dict.items():
    clustered_umis_dict[key]= [str(x) for x in value]
    clustered_umis_dict[key]=[value for value in clustered_umis_dict[key] if value != 'nan']
    clustered_umis_dict[key]= [ ast.literal_eval(x) for x in clustered_umis_dict[key]]
    clustered_umis_dict[key]=[set(x) for x in clustered_umis_dict[key]]


### Read VDJ & Modify ######

####################################
#sample_name='SPTCR16'
#sample_name_extended='GBM275'
keys=list(clustered_umis_dict.keys())
######################################

SPTCR=pmd.read_csv('/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Spatial_Data/SPTCR_Demultiplexed/scTagger/UMI_Extraction/SPTCR_UMI_BARCODE_REVCOMP_CORRECTED/{0}_{1}_vdj_umi_barcode_revcomp_corrected.csv'.format(sample_name,sample_name_extended),index_col=0).dropna(subset=['locus','v_call','j_call','cdr3_aa'])
#pd.read_csv(DF_NAME, converters={'COLUMN_NAME': pd.eval})
################################
SPTCR=SPTCR[SPTCR['Visium Barcode'].isin(keys)]
########################################
print(SPTCR)

### Concat to TCR & get Relevant Columns
SPTCR[['locus','v_call','d_call','j_call','cdr3_aa']]=SPTCR[['locus','v_call','d_call','j_call','cdr3_aa']].astype(str)
SPTCR['TCR']=SPTCR[['locus','v_call','d_call','j_call','cdr3_aa']].agg('_'.join,axis=1)
SPTCR_umi_corrector=SPTCR[['TCR','Visium Barcode','UMI']].reset_index(drop=True)

#Get TCR UMI Table
SPTCR_umi_corrector=SPTCR_umi_corrector.groupby(['TCR','Visium Barcode']).agg({'UMI': lambda x: set(x)})
SPTCR_umi_corrector_transposed=SPTCR_umi_corrector.T

### Get TCR BC Count Table
UMI_uncorrected=SPTCR[['TCR','Visium Barcode','UMI']].reset_index(drop=True)
UMI_uncorrected_count_table=UMI_uncorrected.reset_index(drop=False).groupby(['TCR','Visium Barcode']).size().reset_index(name='UMI uncorrected Count').sort_values(by='UMI uncorrected Count',ascending=False).reset_index(drop=True)


TRA=UMI_uncorrected_count_table[UMI_uncorrected_count_table['TCR'].str.contains('TRA_')]
TRA=TRA.rename(columns={'TCR':'TRA'})

TRB=UMI_uncorrected_count_table[UMI_uncorrected_count_table['TCR'].str.contains('TRB_')]
TRB=TRB.rename(columns={'TCR':'TRB'})

TRA_count=TRA.sort_values(by='UMI uncorrected Count',ascending=False).reset_index(drop=True)

TRB_count=TRB.sort_values(by='UMI uncorrected Count',ascending=False).reset_index(drop=True)

## Scale TRA & TRB Count

TRB_count_tansformed=TRB_count
TRA_count_tansformed=TRA_count

if len(TRB_count['UMI uncorrected Count']) > len(TRA_count['UMI uncorrected Count']):
    #print('TRB')
    MAX=TRB_count['UMI uncorrected Count'].max()
    MIN=TRB_count['UMI uncorrected Count'].min()
else:
    #print('TRA')
    MAX=TRA_count['UMI uncorrected Count'].max()
    MIN=TRA_count['UMI uncorrected Count'].min()
#print(MAX,MIN)

### Min Max Scaler
scaler=sk.MinMaxScaler(feature_range=(MIN, MAX))


## Prepare Data and Scale

Tranformed_TCR_count=pmd.DataFrame()
Tranformed_TCR_count['TRA Unscaled UMI Uncorrected']=TRA_count['UMI uncorrected Count']
Tranformed_TCR_count['TRB Unscaled UMI Uncorrected']=TRB_count['UMI uncorrected Count']

Tranformed_TCR_count[['TRA Transformed','TRB Transformed']]=scaler.fit_transform(Tranformed_TCR_count[['TRA Unscaled UMI Uncorrected','TRB Unscaled UMI Uncorrected']])
Tranformed_TCR_count=Tranformed_TCR_count.round().astype(int)
#display('Tranformed_TCR_count',Tranformed_TCR_count)



### Add Scaled Data to old count DF and cutoff at decided value
TRB_count['Count Transformed']=Tranformed_TCR_count['TRB Transformed'].round().astype(int)
TRB_count['Count Transformed']=TRB_count['Count Transformed'].fillna(0)
TRB_count['Count Transformed']=TRB_count['Count Transformed'].round().astype(int)
#TRB_count=TRB_count[TRB_count['Count Transformed']>=cutoff]

TRA_count['Count Transformed']=Tranformed_TCR_count['TRA Transformed']
TRA_count['Count Transformed']=TRA_count['Count Transformed'].fillna(0)
TRA_count['Count Transformed']=TRA_count['Count Transformed'].round().astype(int)
#TRA_count=TRA_count[TRA_count['Count Transformed']>=cutoff]


TRA_count=TRA_count.rename(columns={'TRA':'TCR'})
TRB_count=TRB_count.rename(columns={'TRB':'TCR'})
UMI_uncorrected_count_table=pmd.concat([TRA_count,TRB_count])

###

### Subset DF only for TCRs with relevant Count
cutoff=5

print('len(UMI_uncorrected_count_table)',len(UMI_uncorrected_count_table))
UMI_uncorrected_count_table=UMI_uncorrected_count_table[UMI_uncorrected_count_table['UMI uncorrected Count']>cutoff]
print('len(UMI_uncorrected_count_table)',len(UMI_uncorrected_count_table))

### Add Column Umi Corrected to fill
UMI_uncorrected_count_table['UMI Corrected']=0


def umi_correct(TCR, BC,clustered_umis_dict=clustered_umis_dict):
    umi_corrected_sum=[]
    example_umi_clusters=clustered_umis_dict[BC]
    for umi_set in SPTCR_umi_corrector_transposed.loc['UMI',TCR]:
        for barcode_umi_set in example_umi_clusters:
            if bool(set(umi_set).intersection(barcode_umi_set)) == True:
                umi_corrected_sum.append(bool(umi_set.intersection(barcode_umi_set)))
            
    umi_sum=sum(umi_corrected_sum)
    return umi_sum

## Apply Function across count table and add UMI Corrected column
UMI_uncorrected_count_table['UMI Corrected']=UMI_uncorrected_count_table.apply(lambda x: umi_correct(TCR=x['TCR'],BC=x['Visium Barcode']),axis='columns')

### Write DF to disk
print(UMI_uncorrected_count_table)   
UMI_uncorrected_count_table.to_csv('/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/Spatial_Data/SPTCR_Demultiplexed/scTagger/UMI_Extraction/UMI_Correction/UMI_CORRECTED/{0}_UMI_corrected_count_table.csv'.format(sample_name),index=False)