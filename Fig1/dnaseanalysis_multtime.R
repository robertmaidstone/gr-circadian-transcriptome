
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

setwd("~/GRanalysis_master/ABC/dnase_data/")
DNase_counts_2<-bin.bw("GSM1479701_WT_DNAse_ZT02.bw", enh_regions_mm9,appendIndex.outfile = FALSE)
DNase_counts_6<-bin.bw("GSM1479702_WT_DNAse_ZT06.bw", enh_regions_mm9,appendIndex.outfile = FALSE)
DNase_counts_10<-bin.bw("GSM1479703_WT_DNAse_ZT10.bw", enh_regions_mm9,appendIndex.outfile = FALSE)
DNase_counts_14<-bin.bw("GSM1479704_WT_DNAse_ZT14.bw", enh_regions_mm9,appendIndex.outfile = FALSE)
DNase_counts_18<-bin.bw("GSM1479705_WT_DNAse_ZT18.bw", enh_regions_mm9,appendIndex.outfile = FALSE)
DNase_counts_22<-bin.bw("GSM1479706_WT_DNAse_ZT22.bw", enh_regions_mm9,appendIndex.outfile = FALSE)



enh_regions_mm9 %>% arrange(chr,start,end) %>% cbind(.,bc2=DNase_counts_2$bc$bc,
                                                     bc6=DNase_counts_6$bc$bc,
                                                     bc10=DNase_counts_10$bc$bc,
                                                     bc14=DNase_counts_14$bc$bc,
                                                     bc18=DNase_counts_18$bc$bc,
                                                     bc22=DNase_counts_22$bc$bc) %>% merge(data_master,by="ENHANCERID") %>%
  rowwise() %>%
  mutate(mean_bc = mean(c_across(bc2:bc22), na.rm = TRUE)) %>%
  ungroup() -> ll_enh

ll_enh %>% mutate(GR=(GRlou_dex_us==0)|(GRlou_WT_us==0)) %>% dplyr::select(ENHANCERID,mean_bc,GR) %>% unique %>% as_tibble %>%
  ggplot(aes(x=mean_bc+.1,fill=GR)) + geom_density(alpha=0.4) + scale_x_log10(lim=c(0.08,100)) + scale_fill_manual(values = c("black","red")) +
  theme_bw() + theme(legend.position = "none") + ylab("Density") + xlab("Normalised Counts") -> p_enhancers

ggsave(plot = p_enhancers,filename = "~/GRanalysis_master/GR_CircadianLiverTranscriptome/Fig2/Plots/dnase_mult_enhancers.png",width=4,height=3)

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

setwd("~/GRanalysis_master/ABC/dnase_data/")
DNase_counts_2<-bin.bw("GSM1479701_WT_DNAse_ZT02.bw", prom_regions_mm9,appendIndex.outfile = FALSE)
DNase_counts_6<-bin.bw("GSM1479702_WT_DNAse_ZT06.bw", prom_regions_mm9,appendIndex.outfile = FALSE)
DNase_counts_10<-bin.bw("GSM1479703_WT_DNAse_ZT10.bw", prom_regions_mm9,appendIndex.outfile = FALSE)
DNase_counts_14<-bin.bw("GSM1479704_WT_DNAse_ZT14.bw", prom_regions_mm9,appendIndex.outfile = FALSE)
DNase_counts_18<-bin.bw("GSM1479705_WT_DNAse_ZT18.bw", prom_regions_mm9,appendIndex.outfile = FALSE)
DNase_counts_22<-bin.bw("GSM1479706_WT_DNAse_ZT22.bw", prom_regions_mm9,appendIndex.outfile = FALSE)

DNase_counts %>% head

prom_regions_mm9 %>% arrange(chr,start,end) %>% cbind(.,bc2=DNase_counts_2$bc$bc,
                                                      bc6=DNase_counts_6$bc$bc,
                                                      bc10=DNase_counts_10$bc$bc,
                                                      bc14=DNase_counts_14$bc$bc,
                                                      bc18=DNase_counts_18$bc$bc,
                                                      bc22=DNase_counts_22$bc$bc)  %>% merge(data_master,by="PROMID")   %>%
  rowwise() %>%
  mutate(mean_bc = mean(c_across(bc2:bc22), na.rm = TRUE)) %>%
  ungroup() -> ll_prom

##
merge(merge(merge((DNase_counts_2$bc %>% mutate(PROMID=paste(chr,start,end,sep="."))),(DNase_counts_6$bc %>% mutate(PROMID=paste(chr,start,end,sep="."))),by="PROMID"),
(DNase_counts_10$bc %>% mutate(PROMID=paste(chr,start,end,sep="."))),by="PROMID"),
(DNase_counts_14$bc %>% mutate(PROMID=paste(chr,start,end,sep="."))),by="PROMID"),
(DNase_counts_18$bc %>% mutate(PROMID=paste(chr,start,end,sep="."))),by="PROMID"),
(DNase_counts_22$bc %>% mutate(PROMID=paste(chr,start,end,sep="."))),by="PROMID") %>% merge(data_master,by="PROMID") -> ll_prom


##

ll_prom %>% mutate(GR=(GRlou_dex_us==0)|(GRlou_WT_us==0)) %>% dplyr::select(PROMID,mean_bc,GR) %>% unique %>% as_tibble %>%
  ggplot(aes(x=mean_bc+.1,fill=GR)) + geom_density(alpha=0.4) + scale_x_log10(lim=c(0.08,100)) + scale_fill_manual(values = c("black","red")) + 
  theme_bw() + theme(legend.position = "none") + ylab("Density") + xlab("Normalised Counts") -> p_promoters

ggsave(plot = p_promoters,filename = "~/GRanalysis_master/GR_CircadianLiverTranscriptome/Fig1/Plots/dnase_mult_promoters.png",width=4,height=3)


save(ll_enh,ll_prom,file = "~/GRanalysis_master/GR_CircadianLiverTranscriptome/Data/dnasedata.RData")

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
  dplyr::select(ENHANCERID,bc=mean_bc,GR,BMAL1,CLOCK,CRY1,CRY2,NPAS2,PER1,PER2) %>% unique %>% as_tibble -> mod_data

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
  dplyr::select(ENHANCERID,bc=mean_bc,GR,BMAL1,CLOCK,CRY1,CRY2,NPAS2,PER1,PER2) %>% unique %>% as_tibble -> mod_data

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

