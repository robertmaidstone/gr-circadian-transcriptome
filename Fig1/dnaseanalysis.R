
library(GenomicRanges)
library(PopSV)
library(tidyverse)
library(effectsize)

load("~/GRanalysis_master/Fig2/data/enhancer_master.RData")

data_master %>% dplyr::select(ENHANCERID) %>% separate(ENHANCERID,into = c("chr","start","end"),remove = FALSE) %>% unique %>%
  mutate(start=as.numeric(start)) %>% mutate(end=as.numeric(end))-> enh_regions

enh_regions <- makeGRangesFromDataFrame(enh_regions,
                                        keep.extra.columns = TRUE,
                                        seqnames.field = "chr",
                                        start.field = "start",
                                        end.field = "end")

enh_regions_mm9<-as.data.frame(easyLift::easyLiftOver(enh_regions,map = "mm10_mm9"))
enh_regions_mm9 %>% dplyr::select(chr=seqnames,start,end,ENHANCERID) -> enh_regions_mm9

setwd("~/GRanalysis_master/ABC/")
DNase_counts<-bin.bw("~/GRanalysis_master/ABC/GSM1479701_WT_DNAse_ZT02.bw", enh_regions_mm9,appendIndex.outfile = FALSE)



enh_regions_mm9 %>% arrange(chr,start,end) %>% cbind(bc=DNase_counts$bc$bc) %>% merge(data_master,by="ENHANCERID") -> ll_enh

ll_enh %>% mutate(GR=(GRlou_dex_us==0)|(GRlou_WT_us==0)) %>% dplyr::select(ENHANCERID,bc,GR) %>% unique %>% as_tibble %>%
  ggplot(aes(x=bc+.1,fill=GR)) + geom_density(alpha=0.4) + scale_x_log10(lim=c(0.08,100)) + scale_fill_manual(values = c("black","red")) +
  theme_bw() + theme(legend.position = "none") + ylab("Density") + xlab("Normalised Counts") -> p_enhancers

ggsave(plot = p_enhancers,filename = "~/GRanalysis_master/GR_CircadianLiverTranscriptome/Fig2/Plots/dnase_enhancers.png",width=4,height=3)

### promoters

load("~/GRanalysis_master/GR_CircadianLiverTranscriptome/Data/RNAChIPtogether_GRdist.RData")

data_master %>% dplyr::select(PROMID) %>% separate(PROMID,into = c("chr","start","end"),remove = FALSE) %>% unique %>%
  mutate(start=as.numeric(start)) %>% mutate(end=as.numeric(end))-> prom_regions

prom_regions %>% filter(!is.na(start)) -> prom_regions

prom_regions <- makeGRangesFromDataFrame(prom_regions,
                                        keep.extra.columns = TRUE,
                                        seqnames.field = "chr",
                                        start.field = "start",
                                        end.field = "end")

prom_regions_mm9<-as.data.frame(easyLift::easyLiftOver(prom_regions,map = "mm10_mm9"))
prom_regions_mm9 %>% dplyr::select(chr=seqnames,start,end,PROMID) -> prom_regions_mm9

DNase_counts<-bin.bw("~/GRanalysis_master/ABC/GSM1479701_WT_DNAse_ZT02.bw", prom_regions_mm9,appendIndex.outfile = FALSE)

DNase_counts %>% head

prom_regions_mm9 %>% arrange(chr,start,end) %>% cbind(bc=DNase_counts$bc$bc) %>% merge(data_master,by="PROMID") -> ll_prom

DNase_counts$bc %>% mutate(PROMID=paste(chr,start,end,sep=".")) %>% merge(data_master,by="PROMID") -> ll_prom

ll_prom %>% mutate(GR=(GRlou_dex_us==0)|(GRlou_WT_us==0)) %>% dplyr::select(PROMID,bc,GR) %>% unique %>% as_tibble %>%
  ggplot(aes(x=bc+.1,fill=GR)) + geom_density(alpha=0.4) + scale_x_log10(lim=c(0.08,100)) + scale_fill_manual(values = c("black","red")) + 
  theme_bw() + theme(legend.position = "none") + ylab("Density") + xlab("Normalised Counts") -> p_promoters

ggsave(plot = p_promoters,filename = "~/GRanalysis_master/GR_CircadianLiverTranscriptome/Fig1/Plots/dnase_promoters.png",width=4,height=3)




ll_enh %>% mutate(GR=(GRlou_dex_us==0)|(GRlou_WT_us==0)) %>% 
  group_by(ENHANCERID) %>%
  mutate(BMAL1=mean(BMAL1,na.rm=T)) %>% 
  mutate(CLOCK=mean(CLOCK,na.rm=T)) %>% 
  mutate(CRY1=mean(CRY1,na.rm=T)) %>% 
  mutate(CRY2=mean(CRY2,na.rm=T)) %>% 
  mutate(NPAS2=mean(NPAS2,na.rm=T)) %>% 
  mutate(PER1=mean(PER1,na.rm=T)) %>% 
  mutate(PER2=mean(PER2,na.rm=T)) %>% 
  #mutate(=mean(CLOCK,na.rm=T)) %>%  #need reverb
  dplyr::select(ENHANCERID,bc,GR,BMAL1,CLOCK,CRY1,CRY2,NPAS2,PER1,PER2) %>% unique %>% as_tibble -> mod_data

lm(BMAL1 ~ bc + GR,data = mod_data) -> model_1
standardize(model_1, method = "refit")  # standardizes predictors and outcome
effectsize::standardize_parameters(model_1, method = "refit")

lm(CLOCK ~ bc + GR,data = mod_data) -> model_1
standardize(model_1, method = "refit")  # standardizes predictors and outcome
effectsize::standardize_parameters(model_1, method = "refit")

lm(CRY1 ~ bc + GR,data = mod_data) -> model_1
standardize(model_1, method = "refit")  # standardizes predictors and outcome
effectsize::standardize_parameters(model_1, method = "refit")

lm(CRY2 ~ bc + GR,data = mod_data) -> model_1
standardize(model_1, method = "refit")  # standardizes predictors and outcome
effectsize::standardize_parameters(model_1, method = "refit")

lm(NPAS2 ~ bc + GR,data = mod_data) -> model_1
standardize(model_1, method = "refit")  # standardizes predictors and outcome
effectsize::standardize_parameters(model_1, method = "refit")

lm(PER1 ~ bc + GR,data = mod_data) -> model_1
standardize(model_1, method = "refit")  # standardizes predictors and outcome
effectsize::standardize_parameters(model_1, method = "refit")

lm(PER2 ~ bc + GR,data = mod_data) -> model_1
standardize(model_1, method = "refit")  # standardizes predictors and outcome
effectsize::standardize_parameters(model_1, method = "refit")

ll_prom %>% mutate(GR=(GRlou_dex_us==0)|(GRlou_WT_us==0)) %>% 
  group_by(ENHANCERID) %>%
  mutate(BMAL1=mean(BMAL1,na.rm=T)) %>% 
  mutate(CLOCK=mean(CLOCK,na.rm=T)) %>% 
  mutate(CRY1=mean(CRY1,na.rm=T)) %>% 
  mutate(CRY2=mean(CRY2,na.rm=T)) %>% 
  mutate(NPAS2=mean(NPAS2,na.rm=T)) %>% 
  mutate(PER1=mean(PER1,na.rm=T)) %>% 
  mutate(PER2=mean(PER2,na.rm=T)) %>% 
  #mutate(=mean(CLOCK,na.rm=T)) %>%  #need reverb
  dplyr::select(ENHANCERID,bc,GR,BMAL1,CLOCK,CRY1,CRY2,NPAS2,PER1,PER2) %>% unique %>% as_tibble -> mod_data

lm(BMAL1 ~ bc + GR,data = mod_data) -> model_1
standardize(model_1, method = "refit")  # standardizes predictors and outcome
effectsize::standardize_parameters(model_1, method = "refit")

lm(CLOCK ~ bc + GR,data = mod_data) -> model_1
standardize(model_1, method = "refit")  # standardizes predictors and outcome
effectsize::standardize_parameters(model_1, method = "refit")

lm(CRY1 ~ bc + GR,data = mod_data) -> model_1
standardize(model_1, method = "refit")  # standardizes predictors and outcome
effectsize::standardize_parameters(model_1, method = "refit")

lm(CRY2 ~ bc + GR,data = mod_data) -> model_1
standardize(model_1, method = "refit")  # standardizes predictors and outcome
effectsize::standardize_parameters(model_1, method = "refit")

lm(NPAS2 ~ bc + GR,data = mod_data) -> model_1
standardize(model_1, method = "refit")  # standardizes predictors and outcome
effectsize::standardize_parameters(model_1, method = "refit")

lm(PER1 ~ bc + GR,data = mod_data) -> model_1
standardize(model_1, method = "refit")  # standardizes predictors and outcome
effectsize::standardize_parameters(model_1, method = "refit")

lm(PER2 ~ bc + GR,data = mod_data) -> model_1
standardize(model_1, method = "refit")  # standardizes predictors and outcome
effectsize::standardize_parameters(model_1, method = "refit")

