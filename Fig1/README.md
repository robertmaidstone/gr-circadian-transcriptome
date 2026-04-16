# Figure 1

Code to recreate figure 1

Fig 1A is a cartoon made using powerpoint, all other plots are created using the following code.

## fig1CE_rhythmicity_analysis.R

Requires:
-   RNAChIPtogether_GRdist.RData
-   pascho.csv - RNA-seq data from Yang et al., 2016
-   conv_pascho.text - geneID to gene symbol conversion file

Creates:
-   RNA_n_forJTK.txt - needed for Metacycle to read, not needed after this script is run
-   pachodata.RData - processed RNA data with JTK results
-   fig 1C - phase plots with stats and venn diagram
-   fig 1E - logistic modelling of rhythmicity adjusting for log10 expression and GR binding
-   sup fig 1A - median gene expression against proportion of rhythmic genes
-   sup fig 1B - logistic modelling of rhythmicity adjusting for GR binding, stratified by deciles of gene expression, values and forest plot

## fig1B_TF_RNA_densityplots.R

Requires:
-   RNAChIPtogether_GRdist.RData
-   promoter_counts_reverb.RData
-   paschodata.RData

Creates:
-   promoter_data_final.RData - ChIP and RNA from Koike et al., and Reverb ChIP together
-   fig1B

## fig1D_clustering_dynamic_profiles.R

Requires:
-   RNAChIPtogether_GRdist.RData"
-   paschodata.RData

Creates:
-   fig 1D

## compareRhythms_fitzgerald_data.R

fig1C using compareRhythms instead of JTK-cycle - this might go into supplemental

Yang G, Chen L, Grant GR, Paschos G, Song WL, Musiek ES, Lee V, McLoughlin SC, Grosser T, Cotsarelis G, FitzGerald GA. Timing of expression of the core clock gene Bmal1 influences its effects on aging and survival. Sci Transl Med. 2016 Feb 3;8(324):324ra16. doi: 10.1126/scitranslmed.aad3305. Epub 2016 Feb 3. PMID: 26843191; PMCID: PMC4870001.
