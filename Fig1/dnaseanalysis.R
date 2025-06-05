
library(GenomicRanges)
library(PopSV)
library(tidyverse)

load("~/GRanalysis_master/Fig2/data/enhancer_master.RData")

data_master %>% dplyr::select(ENHANCERID) %>% separate(ENHANCERID,into = c("chr","start","end")) %>% unique %>%
  mutate(start=as.numeric(start)) %>% mutate(end=as.numeric(end))-> enh_regions

setwd("~/GRanalysis_master/ABC/")
DNase_counts<-bin.bw("~/GRanalysis_master/ABC/GSM1479701_WT_DNAse_ZT02.bw", enh_regions,appendIndex.outfile = FALSE)

DNase_counts %>% head

DNase_counts$bc %>% mutate(ENHANCERID=paste(chr,start,end,sep=".")) %>% merge(data_master,by="ENHANCERID") -> ll_enh

ll_enh %>% mutate(GR=(GRlou_dex_us==0)|(GRlou_WT_us==0)) %>% dplyr::select(ENHANCERID,bc,GR) %>% unique %>% as_tibble %>%
  ggplot(aes(x=bc,fill=GR)) + geom_density(alpha=0.4) + scale_x_log10(lim=c(0.08,10)) + scale_fill_manual(values = c("black","red")) +
  theme_bw() + theme(legend.position = "none") + ylab("Density") + xlab("Normalised Counts") -> p_enhancers

ggsave(plot = p_enhancers,filename = "~/GRanalysis_master/GR_CircadianLiverTranscriptome/Fig2/Plots/dnase_enhancers.png",width=4,height=3)

### promoters

load("~/GRanalysis_master/GR_CircadianLiverTranscriptome/Data/RNAChIPtogether_GRdist.RData")

data_master %>% dplyr::select(PROMID) %>% separate(PROMID,into = c("chr","start","end")) %>% unique %>%
  mutate(start=as.numeric(start)) %>% mutate(end=as.numeric(end))-> prom_regions

prom_regions %>% filter(!is.na(start)) -> prom_regions

DNase_counts<-bin.bw("~/GRanalysis_master/ABC/GSM1479701_WT_DNAse_ZT02.bw", prom_regions,appendIndex.outfile = FALSE)

DNase_counts %>% head

DNase_counts$bc %>% mutate(PROMID=paste(chr,start,end,sep=".")) %>% merge(data_master,by="PROMID") -> ll_prom

ll_prom %>% mutate(GR=(GRlou_dex_us==0)|(GRlou_WT_us==0)) %>% dplyr::select(PROMID,bc,GR) %>% unique %>% as_tibble %>%
  ggplot(aes(x=bc,fill=GR)) + geom_density(alpha=0.4) + scale_x_log10(lim=c(0.08,10)) + scale_fill_manual(values = c("black","red")) + 
  theme_bw() + theme(legend.position = "none") + ylab("Density") + xlab("Normalised Counts") -> p_promoters

ggsave(plot = p_promoters,filename = "~/GRanalysis_master/GR_CircadianLiverTranscriptome/Fig1/Plots/dnase_promoters.png",width=4,height=3)

