##############################################
## Read in VDJ Annotation File
vdj=pmd.read_csv('/media/jkbuntu/SAMSUNG2TB/Dropbox/KBJasim/Projects/Capture_Sequencing/MEASURE_CORRECTION/VDJ_Annotation_per_Read/UNCORRECTED/{0}/{0}_{0}_vdj_annotation.csv'.format(sample),index_col=0)
### Modify ReadID
vdj.index=vdj.index.str.split(pat='|',expand=False,n=1)
vdj=vdj.reset_index(drop=False)
print(vdj)
vdj['index']=vdj['index'].apply(lambda x: str(x[1]).split('runid')[0])
vdj=vdj.rename(columns={'index':'ReadID'})
print('VDJ Annotation File with',vdj.head(5))

## Demultiplex VDJ Annotation File with Barcode and UMI File
demux_umi_df=pmd.merge(left=vdj,right=umi_barcode_df,on='ReadID')
demux_umi_df.to_csv('{0}/{1}_vdj_umi_barcode_uncorrected_df.csv'.format(output_directory,sample