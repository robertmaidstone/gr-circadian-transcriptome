#source("https://bioconductor.org/biocLite.R")
#biocLite("DiffBind")


library(DiffBind)
setwd("~/rds-loudon-ray-rattray/rds-shared-loudon-ray-rattray/Takahashi_data/")

load(file="Enhancer_Counting/enhancer_peaklist_h3k4me1.RData")


ddmm<-dba(sampleSheet="samplesheet_TFs_nocontrol.csv")
ddmm_count<-dba.count(ddmm,peaks = peaks,minOverlap=1)
save(ddmm_count,file="Data/enhancer_counts_TFs_nc_h3k4me1.RData")