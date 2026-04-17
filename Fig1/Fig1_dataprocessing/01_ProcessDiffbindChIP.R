#source("https://bioconductor.org/biocLite.R")
#biocLite("DiffBind")

library(DiffBind)
library(tidyverse)
library(MASS)

setwd("Data")

load("promoter_counts_histonemarks_nov19.RData")
histonemarks_count <- ddmm_count

#load("promoter_counts_TFs_nocontrol.RData")
load("promoter_counts_nov19.RData")

peak_locs<-ddmm_count$peaks[[1]][,1:3]

for(CT in c(0,4,8,12,16,20)){ 
RPKM<-as_tibble(peak_locs)
RPKM$ID<-1:(dim(peak_locs)[1])
RPKM$BMAL1<-ddmm_count$peaks[[(CT/4)+1]][,5]
RPKM$CLOCK<-ddmm_count$peaks[[(CT/4)+7]][,5]
RPKM$CRY1<-ddmm_count$peaks[[(CT/4)+13]][,5]
RPKM$CRY2<-ddmm_count$peaks[[(CT/4)+19]][,5]
RPKM$PER1<-ddmm_count$peaks[[(CT/4)+25]][,5]
RPKM$PER2<-ddmm_count$peaks[[(CT/4)+31]][,5]
RPKM$NPAS2<-ddmm_count$peaks[[(CT/4)+37]][,5]

RPKM$H3K27ac<-histonemarks_count$peaks[[CT/4+1]][,5]
RPKM$H3K4me3<-histonemarks_count$peaks[[CT/4+7]][,5]

RPKM %>% 
  gather(Target,RPKM,-ID,-Chr,-Start,-End) %>%
  dplyr::select(ID,Chr,Start,End,Target,RPKM)-> RPKM

RPKM %>%
  spread(Target,RPKM) -> RPKM_spread

# prom_sites <- GRanges(seqnames = RPKM_spread$Chr,
#                     ranges = IRanges(RPKM_spread$Start, end = RPKM_spread$End))
#
# df <- data.frame(seqnames=seqnames(prom_sites),
#                  starts=start(prom_sites)-1,
#                  ends=end(prom_sites),
#                  names=c(rep(".", length(prom_sites))),
#                  scores=c(rep(".", length(prom_sites))),
#                  strands=strand(prom_sites))
# 
# write.table(df, file="~/Desktop/Sites.bed", quote=F, sep="\t", row.names=F, col.names=F)

############# intersect using bedtools
############# bedtools intersect -bed -wa -a Sites.bed -b mm10.blacklist.bed > BL_proms.bed ####################
BlackRegions<- read.csv("~/GRwholepromoter/BL_proms.bed",sep="\t",header = F)
head(BlackRegions)

RPKM_spread %>%
  filter(!(paste(Chr,Start-1) %in% paste(BlackRegions$V1,BlackRegions$V2))) -> RPKM_filt

###### cRPKM  ########
cRPKM<-as_tibble(peak_locs)
cRPKM$ID<-1:(dim(peak_locs)[1])
cRPKM$BMAL1<-ddmm_count$peaks[[(CT/4)+1]][,7]
cRPKM$CLOCK<-ddmm_count$peaks[[(CT/4)+7]][,7]
cRPKM$CRY1<-ddmm_count$peaks[[(CT/4)+13]][,7]
cRPKM$CRY2<-ddmm_count$peaks[[(CT/4)+19]][,7]
cRPKM$PER1<-ddmm_count$peaks[[(CT/4)+25]][,7]
cRPKM$PER2<-ddmm_count$peaks[[(CT/4)+31]][,7]
cRPKM$NPAS2<-ddmm_count$peaks[[(CT/4)+37]][,7]
cRPKM$H3K27ac<-histonemarks_count$peaks[[CT/4+1]][,7]
cRPKM$H3K4me3<-histonemarks_count$peaks[[CT/4+7]][,7]

cRPKM %>% 
  gather(Target,cRPKM,-ID,-Chr,-Start,-End) %>%
  dplyr::select(ID,Chr,Start,End,Target,cRPKM)-> cRPKM

cRPKM %>%
  spread(Target,cRPKM) -> cRPKM_spread

BlackRegions<- read.csv("~/GRwholepromoter/BL_proms.bed",sep="\t",header = F)
head(BlackRegions)

cRPKM_spread %>%
  filter(!(paste(Chr,Start-1) %in% paste(BlackRegions$V1,BlackRegions$V2))) -> cRPKM_filt

###### READS  ########
READS<-as_tibble(peak_locs)
READS$ID<-1:(dim(peak_locs)[1])
READS$BMAL1<-ddmm_count$peaks[[(CT/4)+1]][,6]
READS$CLOCK<-ddmm_count$peaks[[(CT/4)+7]][,6]
READS$CRY1<-ddmm_count$peaks[[(CT/4)+13]][,6]
READS$CRY2<-ddmm_count$peaks[[(CT/4)+19]][,6]
READS$PER1<-ddmm_count$peaks[[(CT/4)+25]][,6]
READS$PER2<-ddmm_count$peaks[[(CT/4)+31]][,6]
READS$NPAS2<-ddmm_count$peaks[[(CT/4)+37]][,6]

READS$H3K27ac<-histonemarks_count$peaks[[CT/4+1]][,6]
READS$H3K4me3<-histonemarks_count$peaks[[CT/4+7]][,6]

READS %>% 
  gather(Target,Reads,-ID,-Chr,-Start,-End) %>%
  dplyr::select(ID,Chr,Start,End,Target,Reads)-> READS

READS %>%
  spread(Target,Reads) -> READS_spread

BlackRegions<- read.csv("~/GRwholepromoter/BL_proms.bed",sep="\t",header = F)
head(BlackRegions)

READS_spread %>%
  filter(!(paste(Chr,Start-1) %in% paste(BlackRegions$V1,BlackRegions$V2))) -> READS_filt

###### cREADS  ########
cREADS<-as_tibble(peak_locs)
cREADS$ID<-1:(dim(peak_locs)[1])
cREADS$BMAL1<-ddmm_count$peaks[[(CT/4)+1]][,8]
cREADS$CLOCK<-ddmm_count$peaks[[(CT/4)+7]][,8]
cREADS$CRY1<-ddmm_count$peaks[[(CT/4)+13]][,8]
cREADS$CRY2<-ddmm_count$peaks[[(CT/4)+19]][,8]
cREADS$PER1<-ddmm_count$peaks[[(CT/4)+25]][,8]
cREADS$PER2<-ddmm_count$peaks[[(CT/4)+31]][,8]
cREADS$NPAS2<-ddmm_count$peaks[[(CT/4)+37]][,8]

cREADS$H3K27ac<-histonemarks_count$peaks[[CT/4+1]][,8]
cREADS$H3K4me3<-histonemarks_count$peaks[[CT/4+7]][,8]

cREADS %>% 
  gather(Target,cReads,-ID,-Chr,-Start,-End) %>%
  dplyr::select(ID,Chr,Start,End,Target,cReads)-> cREADS

cREADS %>%
  spread(Target,cReads) -> cREADS_spread

BlackRegions<- read.csv("BL_proms.bed",sep="\t",header = F)
head(BlackRegions)

cREADS_spread %>%
  filter(!(paste(Chr,Start-1) %in% paste(BlackRegions$V1,BlackRegions$V2))) -> cREADS_filt

###### saving #######

assign(paste0("RPKM_CT",CT,"_filt"),RPKM_filt)
assign(paste0("cRPKM_CT",CT,"_filt"),cRPKM_filt)
assign(paste0("READS_CT",CT,"_filt"),READS_filt)
assign(paste0("cREADS_CT",CT,"_filt"),cREADS_filt)

save(list=c(paste0("RPKM_CT",CT,"_filt"),paste0("cRPKM_CT",CT,"_filt"),paste0("READS_CT",CT,"_filt"),paste0("cREADS_CT",CT,"_filt")),file = paste0("CT",CT,"data_nc.RData"))
}
