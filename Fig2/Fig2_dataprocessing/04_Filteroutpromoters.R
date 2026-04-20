library(tidyverse)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

setwd("Data")

#### Get Promoters

load("~/GR_promoters_fig1/Data/RNAChIPtogether_GRdist.RData")
data_master %>% dplyr::select(ENHANCERID)  %>%
  separate(ENHANCERID, c("seqnames", "starts", "ends"),sep="[.]") %>% unique -> merged_promoters

write.table(merged_promoters, file="merged_promoters.bed", quote=F, sep="\t", row.names=F, col.names=F)
################
load("ChIPtogether_GRdist.RData")

data_master %>% dplyr::select(ENHANCERID)  %>%
  separate(ENHANCERID, c("seqnames", "starts", "ends"),sep="[.]") %>% unique -> enhancers

write.table(enhancers, file="enhancers.bed", quote=F, sep="\t", row.names=F, col.names=F)
############# find closest using bedtools
############# bedtools sort -i enhancers.bed > enhancers_sort.bed
############# bedtools sort -i merged_promoters.bed > merged_promoters_sort.bed
############# bedtools closest -iu -D a -t first -a enhancers_sort.bed -b merged_promoters_sort.bed > enh_prom_intersect_ds.bed #######################
############# bedtools closest -id -D a -t first -a enhancers_sort.bed -b merged_promoters_sort.bed > enh_prom_intersect_us.bed #######################

enh_prom_intersect <- read.csv("enh_prom_intersect_ds.bed",sep="\t",header = F) 
enh_prom_intersect %>% mutate(ENHANCERID=paste(V1,V2,V3,sep=".")) %>% dplyr::rename(Prom_dist_ds=V7) %>%
  dplyr::select(ENHANCERID,Prom_dist_ds) %>% unique() -> enh_prom_intersect_ds

enh_prom_intersect <- read.csv("enh_prom_intersect_us.bed",sep="\t",header = F) 
enh_prom_intersect %>% mutate(ENHANCERID=paste(V1,V2,V3,sep=".")) %>% dplyr::rename(Prom_dist_us=V7) %>%
  dplyr::select(ENHANCERID,Prom_dist_us) %>% unique() -> enh_prom_intersect_us

merge(data_master,enh_prom_intersect_us,by=c("ENHANCERID")) -> tempK
merge(tempK,enh_prom_intersect_ds,by=c("ENHANCERID")) -> tempK


data_master<-tempK


mm10 = TxDb.Mmusculus.UCSC.mm10.knownGene 

makeGRangesFromDataFrame(enhancers,seqnames.field = "seqnames",start.field = "starts",end.field = "ends") -> enh_granges

peakAnno <- annotatePeak(enh_granges, tssRegion=c(-2000, 400),
                         TxDb=mm10, annoDb="org.Mm.eg.db")

###############

peakAnno %>% as_tibble %>% mutate(ENHANCERID=paste(seqnames,start,end,sep=".")) %>%
  mutate(location=(str_split(annotation,pattern = " ") %>% map(.f = function(x){x[1]}) %>% unlist)) %>%
  dplyr::select(ENHANCERID,distanceToTSS,location) -> peakAnno_t

data_master %>% mutate(Prom=ifelse((Prom_dist_ds<=400)|(Prom_dist_us>=2000),TRUE,FALSE)) %>%
  mutate(Prom_dist=ifelse((Prom_dist_ds<abs(Prom_dist_us))&(Prom_dist_ds!=-1),Prom_dist_ds,abs(Prom_dist_ds))) %>%
  merge(peakAnno_t,by="ENHANCERID") -> data_master_m


data_master_m  -> data_master

save(data_master, file="ChIPtogether_GRdist_locdata.RData")


