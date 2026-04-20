# Figure 2 Data Processing Code

Code to process omics data used for figure 1

## 00_dbacount_h3k4me1.R

Merges h3k4me1 reads to define enhancer regions.

Requires:
-   samplesheetH3K4me1_nocontrol_temp.csv - sample sheet of h3k4me1 data files

Creates:
-   enhancer_peaklist_h3k4me1.RData - enhancer definitions, peaks of h3k4me1 +-500bp

## 00_dbacount_TFs.R

Counts ChIP-seq data onto enhancer regions.

Requires:
-   enhancer_peaklist_h3k4me1.RData - promoter regions for reads to be counted on
-   samplesheet_TFs_nocontrol.csv - sample sheet of data files to be counted (example given in fig1_dataprocessing)

Creates:
-   enhancer_counts_TFs_nc_h3k4me1.RData

## 01_ProcessDiffbindChIP.R

Creates RPKM and read tables for each timepoint. Whilst filtering out black listed regions

Requires:
-   enhancer_counts_TFs_nc_h3k4me1.RData
-   BL_enh.bed - blacklisted regions

Creates:
-   CT0data_nc.RData
-   CT4data_nc.RData
-   CT8data_nc.RData
-   CT12data_nc.RData
-   CT16data_nc.RData
-   CT20data_nc.RData

## 02_CombineChIP.R

combines individual timepoint data into one data frame (and runs JTK-cycle on the h3k27ac)

Requires:
-   CT0data_nc.RData
-   CT4data_nc.RData
-   CT8data_nc.RData
-   CT12data_nc.RData
-   CT16data_nc.RData
-   CT20data_nc.RData
-   BL_enh.bed - blacklisted regions

Creates:
-   ChIP_regions.bed
-   h3k27ac_data_detrend_forJTK.txt  -table used for running JTK-cycle on h3k27ac data
-   ChIPtogether.RData - main output

## 03_IntersectChIPwithGRChIP.R

Adds in GR ChIP-seq data

Requires:
-   ChIPtogether.RData
-   GR_limdata.RData - Lim et al., GR ChIP-seq data, not used in final publication, but was used early as a check.
-   GR_limdata_6am.RData - Lim et al., GR ChIP-seq data, not used in final publication, but was used early as a check.
-   GR_limdata_6pm.RData - Lim et al., GR ChIP-seq data, not used in final publication, but was used early as a check.
-   GR_limdata_pred.RData - Lim et al., GR ChIP-seq data, not used in final publication, but was used early as a check.
-   GR_dex.RData - Novel GR data, dexamethosoline
-   GR_WT.RData - Novel GR data, wild type
-   prom_overlapping_enhancers_sort.bed

Creates:
-   enhancers.bed -> is then sorted to form enhancers_sort.bed which is used for intersections
-   ChIPtogether_GRdist.RData - main output

Other: 
-   Code writes GR regions (from the individual GR data files) into bed files, which are then compared with prom_overlapping_enhancers_sort.bed using bedtools to find closest. This is done upstream and downstream and the resulting files are read in and combined with main data.

## 04_Filteroutpromoters.R