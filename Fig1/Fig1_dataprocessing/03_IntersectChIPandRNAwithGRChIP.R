library(tidyverse)

setwd("~/GR_Nov2019/Data")

load(file="RNAChIPtogether.RData")

data_master %>% dplyr::select(ENHANCERID)  %>%
  separate(ENHANCERID, c("seqnames", "starts", "ends"),sep="[.]") %>% unique -> prom_overlaping_enhancers

write.table(prom_overlaping_enhancers, file="prom_overlapping_enhancers.bed", quote=F, sep="\t", row.names=F, col.names=F)

#### Lim GR

load("GRdata/GR_limdata.RData")
GR_regions <- GR_count$peaks[[1]] %>% dplyr::select(Chr,Start,End)
GR_regions %>% dplyr::rename(seqnames=Chr,starts=Start,ends=End) -> GR_regions

write.table(GR_regions, file="GR_lim_regions.bed", quote=F, sep="\t", row.names=F, col.names=F)

############# find closest using bedtools
############# bedtools sort -i prom_overlapping_enhancers.bed > prom_overlapping_enhancers_sort.bed
############# bedtools closest -iu -D a -t first -a prom_overlapping_enhancers_sort.bed -b GR_lim_regions.bed > GR_lim_intersect_ds.bed #######################
############# bedtools closest -id -D a -t first -a prom_overlapping_enhancers_sort.bed -b GR_lim_regions.bed > GR_lim_intersect_us.bed #######################

rhythmic_GR_int <- read.csv("GR_lim_intersect_us.bed",sep="\t",header = F) 
rhythmic_GR_int %>% mutate(ENHANCERID=paste(V1,V2,V3,sep=".")) %>% dplyr::rename(GRlim_us=V7) %>%
  dplyr::select(ENHANCERID,GRlim_us) %>% unique() -> rhythmic_GR_int_us

rhythmic_GR_int <- read.csv("GR_lim_intersect_ds.bed",sep="\t",header = F) 
rhythmic_GR_int %>% mutate(ENHANCERID=paste(V1,V2,V3,sep=".")) %>% dplyr::rename(GRlim_ds=V7) %>%
  dplyr::select(ENHANCERID,GRlim_ds) %>% unique() -> rhythmic_GR_int_ds

merge(data_master,rhythmic_GR_int_us,by=c("ENHANCERID")) -> tempK
merge(tempK,rhythmic_GR_int_ds,by=c("ENHANCERID")) -> tempK

#### Lim GR (6am)

load("GRdata/GR_limdata_6am.RData")
GR_regions <- GR_count$peaks[[1]] %>% dplyr::select(Chr,Start,End)
GR_regions %>% dplyr::rename(seqnames=Chr,starts=Start,ends=End) -> GR_regions

write.table(GR_regions, file="GR_lim_regions_6am.bed", quote=F, sep="\t", row.names=F, col.names=F)

############# find closest using bedtools
############# bedtools sort -i prom_overlapping_enhancers.bed > prom_overlapping_enhancers_sort.bed
############# bedtools closest -iu -D a -t first -a prom_overlapping_enhancers_sort.bed -b GR_lim_regions_6am.bed > GR_lim_intersect_6am_ds.bed #######################
############# bedtools closest -id -D a -t first -a prom_overlapping_enhancers_sort.bed -b GR_lim_regions_6am.bed > GR_lim_intersect_6am_us.bed #######################

rhythmic_GR_int <- read.csv("GR_lim_intersect_6am_us.bed",sep="\t",header = F) 
rhythmic_GR_int %>% mutate(ENHANCERID=paste(V1,V2,V3,sep=".")) %>% dplyr::rename(GRlim_6am_us=V7) %>%
  dplyr::select(ENHANCERID,GRlim_6am_us) %>% unique() -> rhythmic_GR_int_us

rhythmic_GR_int <- read.csv("GR_lim_intersect_6am_ds.bed",sep="\t",header = F) 
rhythmic_GR_int %>% mutate(ENHANCERID=paste(V1,V2,V3,sep=".")) %>% dplyr::rename(GRlim_6am_ds=V7) %>%
  dplyr::select(ENHANCERID,GRlim_6am_ds) %>% unique() -> rhythmic_GR_int_ds

merge(tempK,rhythmic_GR_int_us,by=c("ENHANCERID")) -> tempK
merge(tempK,rhythmic_GR_int_ds,by=c("ENHANCERID")) -> tempK

#### Lim GR (6pm)

load("GRdata/GR_limdata_6pm.RData")
GR_regions <- GR_count$peaks[[1]] %>% dplyr::select(Chr,Start,End)
GR_regions %>% dplyr::rename(seqnames=Chr,starts=Start,ends=End) -> GR_regions

write.table(GR_regions, file="GR_lim_regions_6pm.bed", quote=F, sep="\t", row.names=F, col.names=F)

############# find closest using bedtools
############# bedtools sort -i prom_overlapping_enhancers.bed > prom_overlapping_enhancers_sort.bed
############# bedtools closest -iu -D a -t first -a prom_overlapping_enhancers_sort.bed -b GR_lim_regions_6pm.bed > GR_lim_intersect_6pm_ds.bed #######################
############# bedtools closest -id -D a -t first -a prom_overlapping_enhancers_sort.bed -b GR_lim_regions_6pm.bed > GR_lim_intersect_6pm_us.bed #######################

rhythmic_GR_int <- read.csv("GR_lim_intersect_6pm_us.bed",sep="\t",header = F) 
rhythmic_GR_int %>% mutate(ENHANCERID=paste(V1,V2,V3,sep=".")) %>% dplyr::rename(GRlim_6pm_us=V7) %>%
  dplyr::select(ENHANCERID,GRlim_6pm_us) %>% unique() -> rhythmic_GR_int_us

rhythmic_GR_int <- read.csv("GR_lim_intersect_6pm_ds.bed",sep="\t",header = F) 
rhythmic_GR_int %>% mutate(ENHANCERID=paste(V1,V2,V3,sep=".")) %>% dplyr::rename(GRlim_6pm_ds=V7) %>%
  dplyr::select(ENHANCERID,GRlim_6pm_ds) %>% unique() -> rhythmic_GR_int_ds

merge(tempK,rhythmic_GR_int_us,by=c("ENHANCERID")) -> tempK
merge(tempK,rhythmic_GR_int_ds,by=c("ENHANCERID")) -> tempK

#### Lim GR (pred)

load("GRdata/GR_limdata_pred.RData")
GR_regions <- GR_count$peaks[[1]] %>% dplyr::select(Chr,Start,End)
GR_regions %>% dplyr::rename(seqnames=Chr,starts=Start,ends=End) -> GR_regions

write.table(GR_regions, file="GR_lim_pred_regions.bed", quote=F, sep="\t", row.names=F, col.names=F)

############# find closest using bedtools
############# bedtools sort -i prom_overlapping_enhancers.bed > prom_overlapping_enhancers_sort.bed
############# bedtools closest -iu -D a -t first -a prom_overlapping_enhancers_sort.bed -b GR_lim_pred_regions.bed > GR_lim_pred_intersect_ds.bed #######################
############# bedtools closest -id -D a -t first -a prom_overlapping_enhancers_sort.bed -b GR_lim_pred_regions.bed > GR_lim_pred_intersect_us.bed #######################

rhythmic_GR_int <- read.csv("GR_lim_pred_intersect_us.bed",sep="\t",header = F) 
rhythmic_GR_int %>% mutate(ENHANCERID=paste(V1,V2,V3,sep=".")) %>% dplyr::rename(GRlim_pred_us=V7) %>%
  dplyr::select(ENHANCERID,GRlim_pred_us) %>% unique() -> rhythmic_GR_int_us

rhythmic_GR_int <- read.csv("GR_lim_pred_intersect_ds.bed",sep="\t",header = F) 
rhythmic_GR_int %>% mutate(ENHANCERID=paste(V1,V2,V3,sep=".")) %>% dplyr::rename(GRlim_pred_ds=V7) %>%
  dplyr::select(ENHANCERID,GRlim_pred_ds) %>% unique() -> rhythmic_GR_int_ds

merge(tempK,rhythmic_GR_int_us,by=c("ENHANCERID")) -> tempK
merge(tempK,rhythmic_GR_int_ds,by=c("ENHANCERID")) -> tempK

# #### Louise GR dex
# 
load("GRdata/GR_dex.RData")
GR_regions <- GR_count$peaks[[1]] %>% dplyr::select(Chr,Start,End)
GR_regions %>% dplyr::rename(seqnames=Chr,starts=Start,ends=End) -> GR_regions

write.table(GR_regions, file="GR_lou_dex_regions.bed", quote=F, sep="\t", row.names=F, col.names=F)
# 
# ############# intersect using bedtools
############# bedtools closest -iu -D a -t first -a prom_overlapping_enhancers_sort.bed -b GR_lou_dex_regions.bed > GR_lou_dex_intersect_ds.bed #######################
############# bedtools closest -id -D a -t first -a prom_overlapping_enhancers_sort.bed -b GR_lou_dex_regions.bed > GR_lou_dex_intersect_us.bed #######################

rhythmic_GR_int <- read.csv("GR_lou_dex_intersect_us.bed",sep="\t",header = F) 
rhythmic_GR_int %>% mutate(ENHANCERID=paste(V1,V2,V3,sep=".")) %>% dplyr::rename(GRlou_dex_us=V7) %>%
  dplyr::select(ENHANCERID,GRlou_dex_us) %>% unique() -> rhythmic_GR_int_us

rhythmic_GR_int <- read.csv("GR_lou_dex_intersect_ds.bed",sep="\t",header = F) 
rhythmic_GR_int %>% mutate(ENHANCERID=paste(V1,V2,V3,sep=".")) %>% dplyr::rename(GRlou_dex_ds=V7) %>%
  dplyr::select(ENHANCERID,GRlou_dex_ds) %>% unique() -> rhythmic_GR_int_ds

merge(tempK,rhythmic_GR_int_us,by=c("ENHANCERID")) -> tempK
merge(tempK,rhythmic_GR_int_ds,by=c("ENHANCERID")) -> tempK
# 
# 
# #### Louise GR WT
# 
load("GRdata/GR_WT.RData")
GR_regions <- GR_count$peaks[[1]] %>% dplyr::select(Chr,Start,End)
GR_regions %>% dplyr::rename(seqnames=Chr,starts=Start,ends=End) -> GR_regions

write.table(GR_regions, file="GR_lou_WT_regions.bed", quote=F, sep="\t", row.names=F, col.names=F)
# 
# ############# intersect using bedtools
############# bedtools closest -iu -D a -t first -a prom_overlapping_enhancers_sort.bed -b GR_lou_WT_regions.bed > GR_lou_WT_intersect_ds.bed #######################
############# bedtools closest -id -D a -t first -a prom_overlapping_enhancers_sort.bed -b GR_lou_WT_regions.bed > GR_lou_WT_intersect_us.bed #######################

rhythmic_GR_int <- read.csv("GR_lou_WT_intersect_us.bed",sep="\t",header = F) 
rhythmic_GR_int %>% mutate(ENHANCERID=paste(V1,V2,V3,sep=".")) %>% dplyr::rename(GRlou_WT_us=V7) %>%
  dplyr::select(ENHANCERID,GRlou_WT_us) %>% unique() -> rhythmic_GR_int_us

rhythmic_GR_int <- read.csv("GR_lou_WT_intersect_ds.bed",sep="\t",header = F) 
rhythmic_GR_int %>% mutate(ENHANCERID=paste(V1,V2,V3,sep=".")) %>% dplyr::rename(GRlou_WT_ds=V7) %>%
  dplyr::select(ENHANCERID,GRlou_WT_ds) %>% unique() -> rhythmic_GR_int_ds

merge(tempK,rhythmic_GR_int_us,by=c("ENHANCERID")) -> tempK
merge(tempK,rhythmic_GR_int_ds,by=c("ENHANCERID")) -> tempK


data_master<-tempK

save(data_master, file="RNAChIPtogether_GRdist.RData")


