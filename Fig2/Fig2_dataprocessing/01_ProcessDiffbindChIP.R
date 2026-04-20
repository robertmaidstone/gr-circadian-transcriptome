#source("https://bioconductor.org/biocLite.R")
#biocLite("DiffBind")

library(DiffBind)
library(tidyverse)
library(MASS)
setwd("Data")

load("enhancer_counts_TFs_nc_h3k4me1.RData")
peak_locs<-ddmm_count$peaks[[1]][,1:3]
####################

prom_sites <- GRanges(seqnames = peak_locs$Chr,
                    ranges = IRanges(peak_locs$Start, end = peak_locs$End))

df <- data.frame(seqnames=seqnames(prom_sites),
                 starts=start(prom_sites)-1,
                 ends=end(prom_sites),
                 names=c(rep(".", length(prom_sites))),
                 scores=c(rep(".", length(prom_sites))),
                 strands=strand(prom_sites))

write.table(df, file="Sites.bed", quote=F, sep="\t", row.names=F, col.names=F)

############# intersect using bedtools
############# bedtools intersect -bed -wa -a Sites.bed -b mm10.blacklist.bed > BL_enh.bed ####################


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
RPKM$H3K27ac<-ddmm_count$peaks[[(CT/4)+43]][,5]
RPKM$H3K4me1<-ddmm_count$peaks[[(CT/4)+49]][,5]

RPKM %>% 
  gather(Target,RPKM,-ID,-Chr,-Start,-End) %>%
  dplyr::select(ID,Chr,Start,End,Target,RPKM)-> RPKM

RPKM %>%
  spread(Target,RPKM) -> RPKM_spread


BlackRegions<- read.csv("BL_enh.bed",sep="\t",header = F)

RPKM_spread %>%
  filter(!(paste(Chr,Start-1) %in% paste(BlackRegions$V1,BlackRegions$V2))) -> RPKM_filt

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
READS$H3K27ac<-ddmm_count$peaks[[(CT/4)+43]][,6]
READS$H3K4me1<-ddmm_count$peaks[[(CT/4)+49]][,6]

READS %>% 
  gather(Target,Reads,-ID,-Chr,-Start,-End) %>%
  dplyr::select(ID,Chr,Start,End,Target,Reads)-> READS

READS %>%
  spread(Target,Reads) -> READS_spread

READS_spread %>%
  filter(!(paste(Chr,Start-1) %in% paste(BlackRegions$V1,BlackRegions$V2))) -> READS_filt

###### saving #######

assign(paste0("RPKM_CT",CT,"_filt"),RPKM_filt)
assign(paste0("READS_CT",CT,"_filt"),READS_filt)

save(list=c(paste0("RPKM_CT",CT,"_filt"),paste0("READS_CT",CT,"_filt")),file = paste0("CT",CT,"data_nc.RData"))
}

