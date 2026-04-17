# Figure 1 Data Processing Code

Code to process omics data used for figure 1

## 00_Diffbind_dbacount_TFs.R

Counts RNA-seq data onto promoter regions.

Requires:
-   RNA_promoters.bed - promoter regions for reads to be counted on
-   samplesheet_TFs_nocontrol.csv - sample sheet of data files to be counted (example given in repository)

Creates:
-   promoter_counts_nov19.RData 
-   promoter_counts_histonemarks_nov19.RData -with different sample sheet

## 01_ProcessDiffbindChIP.R

Creates RPKM and read tables as well as control corrected versions for each timepoint. Whilst filtering out black listed regions

Requires:
-   promoter_counts_histonemarks_nov19.RData
-   promoter_counts_nov19.RData
-   BL_proms.bed - blacklisted regions

Creates:
-   CT0data_nc.RData
-   CT4data_nc.RData
-   CT8data_nc.RData
-   CT12data_nc.RData
-   CT16data_nc.RData
-   CT20data_nc.RData

## 02_CombineChIPandRNA.R

Creates RPKM and read tables as well as control corrected versions for each timepoint. Whilst filtering out black listed regions

Requires:
-   CT0data_nc.RData
-   CT4data_nc.RData
-   CT8data_nc.RData
-   CT12data_nc.RData
-   CT16data_nc.RData
-   CT20data_nc.RData
-   BL_proms.bed - blacklisted regions
-   RNAseqFC.RData - RNAseq data from Koike et al.

Creates:
-   ChIP_regions.bed
-   FC_RNA_edgeRnorm_forJTK.txt -table used for running JTK-cycle
-   RNA_promoters_edgeRnorm.bed - promoters for finding overlaps with ChIP regions
-   ChIP_RNAproms_intersect.bed - overlaps created using bedtools (command given in R script)
-   RNAproms_ChIP_intersect.bed - overlaps created using bedtools (command given in R script)
-   RNAChIPtogether.RData - main output

## 03_IntersectChIPandRNAwithGRChIP.R

Adds in GR ChIP-seq data

Requires:
-   RNAChIPtogether.RData
-   GR_limdata.RData - Lim et al., GR ChIP-seq data, not used in final publication, but was used early as a check.
-   GR_limdata_6am.RData - Lim et al., GR ChIP-seq data, not used in final publication, but was used early as a check.
-   GR_limdata_6pm.RData - Lim et al., GR ChIP-seq data, not used in final publication, but was used early as a check.
-   GR_limdata_pred.RData - Lim et al., GR ChIP-seq data, not used in final publication, but was used early as a check.
-   GR_dex.RData - Novel GR data, dexamethosoline
-   GR_WT.RData - Novel GR data, wild type
-   prom_overlapping_enhancers_sort.bed

Creates:
-   RNAChIPtogether_GRdist.RData - main output

Other: 
-   Code writes GR regions (from the individual GR data files) into bed files, which are then compared with prom_overlapping_enhancers_sort.bed using bedtools to find closest. This is done upstream and downstream and the resulting files are read in and combined with main data.
