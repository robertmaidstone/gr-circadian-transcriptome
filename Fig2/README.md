# Figure 2

Code to recreate figure 2.
Fig 2A is a cartoon made using powerpoint, all other plots are created using the following code.


## fig2DE_rhythmicity_analysis.R

Requires:
-   enhancer_withRNA.RData
-   pascho.csv - RNA-seq data from Yang et al., 2016
-   conv_pascho.text - geneID to gene symbol conversion file

Creates:
-   RNA_n_forJTK.txt - needed for Metacycle to read, not needed after this script is run
-   paschodata_enhancer.RData - processed RNA data with JTK results
-   fig 2D - phase plots with stats and venn diagram
-   fig 2E - boxplot comparing gene expression across groups
-   GRboundRhythmicenhancers_forHomer.bed - used for fig 2F
-   sup fig 2A - logistic modelling of rhythmicity adjusting for log10 expression and GR binding
-   sup fig 2B - median gene expression against proportion of rhythmic genes
-   sup fig 2C - logistic modelling of rhythmicity adjusting for GR binding, stratified by deciles of gene expression, values and forest plot


## rhythmicity_analysis_fitzgerald_data.R

Figure 2D and E and bed file creation for use by Homer to generate figure 2F. Modelling of rhythmicity for supplemental figure 2A,B and C

## genomiclocationplots.R

fig 2B

## compareRhythms_fitzgerald_data.R

fig2D using compareRhythms instead of JTK-cycle - this might go into supplemental

## supfig1C_dnaseanalysis_multtime_rev.R

Requires:
-   enhancer_data_final_rev.RData
-   GSM1479701_WT_DNAse_ZT02.bw
-   GSM1479702_WT_DNAse_ZT06.bw
-   GSM1479703_WT_DNAse_ZT10.bw
-   GSM1479704_WT_DNAse_ZT14.bw
-   GSM1479705_WT_DNAse_ZT18.bw
-   GSM1479706_WT_DNAse_ZT22.bw