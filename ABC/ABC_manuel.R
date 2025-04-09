library(GenomicRanges)
library(PopSV)
library(tidyverse)

load("~/GRanalysis_master/Fig2/data/enhancer_master.RData")

data_master %>% dplyr::select(ENHANCERID) %>% separate(ENHANCERID,into = c("chr","start","end")) %>% unique %>%
  mutate(start=as.numeric(start)) %>% mutate(end=as.numeric(end))-> enh_regions

setwd("~/GRanalysis_master/ABC/")
DNase_counts<-bin.bw("~/GRanalysis_master/ABC/GSM1479701_WT_DNAse_ZT02.bw", enh_regions,appendIndex.outfile = FALSE)

DNase_counts %>% head

enh_regions_gr <- GenomicRanges::makeGRangesFromDataFrame(enh_regions)

proms<-read.table(file="~/GRanalysis_master/Fig2/data/RNA_promoters_V2.bed")
names(proms)<-c("chr","start","end")

proms_gr <- GenomicRanges::makeGRangesFromDataFrame(proms)

## chr 1

DNase_counts$bc %>% filter(chr=="chr1") -> dnase_chr
enh_regions_gr <- GenomicRanges::makeGRangesFromDataFrame(dnase_chr,keep.extra.columns = TRUE)

proms_gr <- GenomicRanges::makeGRangesFromDataFrame(proms %>% filter(chr=="chr1"))

denom <- c()
for(i in 1:length(proms_gr)){
  GenomicRanges::distance(proms_gr[i],enh_regions_gr,select="all") -> ll
  temp<-as.data.frame(enh_regions_gr) %>% mutate(dist=ll) #%>%GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  temp %>% filter(dist<5000000) %>% mutate(mult=bc*(dist)^(-1)) %>%dplyr::select(mult) %>% sum -> denom[i]
}

ABC_df <- data.frame(seqnames=NA,   start=NA,     end=NA, width=NA, strand=NA,   dist=NA,        denom =NA,      ABC=NA,enh_i=NA)
for(i in 1:length(enh_regions_gr)){
  dnase<-enh_regions_gr[i]$bc
  GenomicRanges::distance(enh_regions_gr[i],proms_gr,select="all") -> ll
  temp<-as.data.frame(proms_gr) %>% mutate(dist=ll) %>% mutate(denom=denom) %>% mutate(ABC=(dist)^(-1)*dnase/denom) #%>%GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  temp %>% filter(ABC==max(temp$ABC,na.rm = T)) %>% mutate(enh_i=i) %>% rbind(ABC_df,.) -> ABC_df
}

ABC_df %>% filter(ABC!=0) %>% dplyr::select(ABC) %>% unlist %>% hist

ABC_df %>% dplyr::select(enh_i,ABC) %>% unique %>% dplyr::select(ABC) %>% unlist %>% hist

## chr 1 - with HiC

load("~/GRanalysis_master/ABC/hic_data.RData")
ordered_granges_enh <- enh_regions_gr[order(seqnames(enh_regions_gr), start(enh_regions_gr), end(enh_regions_gr))]
ordered_granges_prom <- proms_gr[order(seqnames(proms_gr), start(proms_gr), end(proms_gr))]

DNase_counts$bc %>% filter(chr=="chr1") -> dnase_chr


enh_regions_gr <- GenomicRanges::makeGRangesFromDataFrame(dnase_chr,keep.extra.columns = TRUE)
mcols(enh_regions_gr)$hic_region <- ordered_granges_enh[seqnames(ordered_granges_enh) == "chr1"]$hic_region

proms_gr <- GenomicRanges::makeGRangesFromDataFrame(proms %>% filter(chr=="chr1"))
mcols(proms_gr)$hic_region <- ordered_granges_prom[seqnames(ordered_granges_prom) == "chr1"]$hic_region

denom <- c()
for(i in 156:length(proms_gr)){
  GenomicRanges::distance(proms_gr[i],enh_regions_gr,select="all") -> ll
  temp<-as.data.frame(enh_regions_gr) %>% mutate(dist=ll) %>% filter(dist<5000000) #%>%GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  if(!is.na(proms_gr[i]$hic_region)){
  data_f %>% filter(data_f$V1%in%temp$hic_region,data_f$V2==proms_gr[i]$hic_region) -> df_hic
  if(dim(df_hic)[1]==0){denom[i]<-NA}else{
  temp %>% filter(hic_region %in% df_hic$V1) %>%
    merge(.,df_hic,by.x="hic_region",by.y="V1") %>%
    mutate(mult=bc*V3) %>%
    dplyr::select(mult) %>% sum -> denom[i]
  }
  }
}

save(denom, file = "~/GRanalysis_master/ABC/denom.RData")

ABC_df <- data.frame(seqnames=NA,   start=NA,     end=NA, width=NA, strand=NA,   dist=NA,        denom =NA,      ABC=NA,enh_i=NA)
for(i in 1:length(enh_regions_gr)){
  dnase<-enh_regions_gr[i]$bc
  GenomicRanges::distance(enh_regions_gr[i],proms_gr,select="all") -> ll
  temp<-as.data.frame(proms_gr) %>% mutate(dist=ll) %>% mutate(denom=denom)
  
  if(!is.na(enh_regions_gr[i]$hic_region)){
   
     data_f %>% filter(data_f$V1%in%enh_regions_gr[i]$hic_region,data_f$V2==temp$hic_region) -> df_hic
    if(dim(df_hic)[1]==0){denom[i]<-NA}else{
      temp %>% filter(hic_region %in% df_hic$V2) %>%
        merge(.,df_hic,by.x="hic_region",by.y="V2") %>%
        mutate(ABC=ABC=V3*dnase/denom) -> temp
      
      temp %>% filter(ABC==max(temp$ABC,na.rm = T)) %>% mutate(enh_i=i) %>% rbind(ABC_df,.) -> ABC_df
  }
}

ABC_df %>% filter(ABC!=0) %>% dplyr::select(ABC) %>% unlist %>% hist

ABC_df %>% dplyr::select(enh_i,ABC) %>% unique %>% dplyr::select(ABC) %>% unlist %>% hist




