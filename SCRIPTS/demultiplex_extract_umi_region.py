
import pandas as pd 
import modin.pandas as pmd
import os
from collections import defaultdict
from Bio.Seq import Seq

import argparse
from unicodedata import name

parser = argparse.ArgumentParser()
parser.add_argument('-i','--INPUT',help="Output Tsv from scTagger", required=True )
parser.add_argument('-u','--UMI_FOLDER',help="Input Folder for Umi Extraction files from_demux_umi_dfs.sh", required=True )
parser.add_argument('-o','--OUT',help="Path to Outfolder", required=True )
parser.add_argument('-n','--NAME',help="Name of the Pipeline", required=True )
parser.add_argument('-b','--BCs',help="Path to Reference Barcodes", required=False,default='../REFERENCE/Barcodes/visium_bc.tsv' )
#parser.add_argument('-igb','--IGB',help="Path to IGB Folder with Initial IG Blast TSV", required=True )

args = parser.parse_args()
arg_vars = vars(args)
arg_vars["NAME"]

INPUT=arg_vars["INPUT"]
output_directory=arg_vars["OUT"]
umi_folder=arg_vars['UMI_FOLDER']
sample=arg_vars['NAME']
BARCODES=arg_vars["BCs"]

#######################################################
print('::::{0} Demultiplexing::::'.format(sample))
UMI_TSVs_dict=defaultdict(list)
Demux_Dict=defaultdict(list)

### Demultiplexing
demux_frame=pmd.read_csv(INPUT,names=['ReadID','StringDist','Dist','Adapter','Barcode'],usecols=['ReadID','Barcode'],sep='\t')

demux_frame=pmd.concat([demux_frame,demux_frame['Barcode'].str.split(',',expand=True)],axis=1)
demux_frame_all=demux_frame
demux_frame['Barcode']=demux_frame[0] ### Choose randomly first Barcode -> can be optimized to choose spatially closest
demux_frame=demux_frame[['ReadID','Barcode']]

out_name='{0}_DEMUX/{0}_all_barcode_matches.csv'.format(sample)
write_path=os.path.join(output_directory,out_name)
demux_frame_all.to_csv(write_path,index=False)
print(demux_frame)


print('::::{0} Add UMI Column::::'.format(sample))
#### UMI TSV

#### For UMI Containing Part
BC_adap_umi=str(umi_folder)+'/'+'{0}_Extracted_Adapter_Barcode_UMI.tsv'.format(sample)

adapter_bc_Umi=pmd.read_csv(BC_adap_umi,sep='\t',names=['ReadID','StringDist','Start','Adapter Sequence'])
adapter_bc_Umi['Adapter Len']=adapter_bc_Umi['Adapter Sequence'].str.len()
adapter_bc_Umi=adapter_bc_Umi.dropna()
print('Adapter and Barcode and UMI',adapter_bc_Umi.head(5))

#### For not UMI Containing Part
BC_adap=str(umi_folder)+'/'+'{0}_Extracted_Adapter_Barcode.tsv'.format(sample)
adapter_bc=pmd.read_csv(BC_adap,sep='\t',names=['ReadID','StringDist','Start','Adapter Sequence'])
adapter_bc['Adapter Len']=adapter_bc['Adapter Sequence'].str.len()
adapter_bc=adapter_bc.dropna()
print('Adapter and Barcode',adapter_bc.head(5))


### Get UMI Column
#for sample, umi_bc_dfs in UMI_TSVs_dict.items():

print(':::: {0} ::::'.format(sample))

print ('--- Merging UMI Extraction Dfs ---')
UMI_ISOLATION=pmd.merge(left=adapter_bc_Umi,right=adapter_bc,on='ReadID',how='inner')
print(UMI_ISOLATION.head(5))

print ('--- Substracting Barcode-UMI Regions ---')
Substract_Adapters=UMI_ISOLATION[['ReadID','Adapter Sequence_x','Adapter Sequence_y']]
Substract_Adapters=Substract_Adapters.rename(columns={'Adapter Sequence_x':'Adapter with BC and UMI','Adapter Sequence_y':'Adapter with BC'})
Substract_Adapters['UMI']=Substract_Adapters.apply(lambda x: x["Adapter with BC and UMI"].replace(x["Adapter with BC"], "").strip(), axis=1)

### Merge UMI DF and Barcode DF
UMI_DF=Substract_Adapters[['ReadID','UMI']]

#demux_df=Demux_Dict[sample]

umi_barcode_df=pmd.merge(left=UMI_DF,right=demux_frame,on='ReadID')


#Demux_Dict[sample]=umi_barcode_df

print('Merged UMI Barcode DF {0}'.format(sample))
print(umi_barcode_df.head(5))

############# Revcomp Correct Barcodes#############


Visium_Bcs=pd.read_csv(BARCODES,sep='\t',names=['Spatial Barcode'])

### Get All forms of Visium BC and concat to new Demultiplexing Frame
Visium_Bcs['Forward Indicator']='Forward'
Visium_Bcs['Complement']=Visium_Bcs['Spatial Barcode'].apply(lambda x: ''.join(list(Seq(x).complement())))
Visium_Bcs['Complement Indicator']='Complement'
Visium_Bcs['Reverse']=Visium_Bcs['Spatial Barcode'].apply(lambda x: ''.join(list(Seq(x).reverse_complement())))
Visium_Bcs['Reverse']=Visium_Bcs['Reverse'].apply(lambda x: ''.join(list(Seq(x).complement())))
Visium_Bcs['Reverse Indicator']='Reverse'
Visium_Bcs['Reverse Complement']=Visium_Bcs['Spatial Barcode'].apply(lambda x: ''.join(list(Seq(x).reverse_complement())))
Visium_Bcs['Reverse Complement Indicator']='Reverse Complement'

concat1=pmd.DataFrame()
concat2=pmd.DataFrame()
concat3=pmd.DataFrame()
concat4=pmd.DataFrame()

concat1[['Spatial Barcode','Indicator']]=Visium_Bcs[['Spatial Barcode','Forward Indicator']]
concat2[['Spatial Barcode','Indicator']]=Visium_Bcs[['Complement','Complement Indicator']]
concat3[['Spatial Barcode','Indicator']]=Visium_Bcs[['Reverse Complement','Reverse Complement Indicator']]
concat4[['Spatial Barcode','Indicator']]=Visium_Bcs[['Reverse','Reverse Indicator']]

test_all_forms=pmd.concat([concat1,concat2,concat3,concat4],axis=0)

## Demultiplex by Merging Both Dataframes
demux=pmd.merge(left=test_all_forms,left_on='Spatial Barcode',right=umi_barcode_df, right_on=['Barcode'],indicator=True)

## Split into Forward and Reverse Barcode Matches
forward_barcodes=demux[demux['Indicator']=='Forward'][['Spatial Barcode','ReadID']]
reverse_barcodes=demux[demux['Indicator']=='Reverse Complement'][['Spatial Barcode','ReadID','UMI']]
print(reverse_barcodes)

## Reverse Complement the Reverse UMIs
reverse_barcodes['UMI']=reverse_barcodes['UMI'].astype(str)                                          
reverse_barcodes['UMI']=reverse_barcodes['UMI'].apply(lambda x: ''.join(list(Seq(x).reverse_complement())),axis=1)


## Merge Both Dataframes Again
demux_frame=pmd.concat([reverse_barcodes,forward_barcodes],axis=0)
print(demux_frame)

### Write Demultiplexed Dataframe
out_name='{0}_barcode_umi.csv'.format(sample)
write_path=os.path.join(output_directory,out_name)
demux_frame.to_csv(write_path, index=False)


