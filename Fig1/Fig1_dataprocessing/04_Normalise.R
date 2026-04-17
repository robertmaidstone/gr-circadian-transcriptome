
library(tidyverse)
library(mclust)
setwd("Data")

load(file="RNAChIPtogether_GRdist.RData")
#load(file="C:/Users/mqbssrmj/Dropbox/Gene_Regulatory_Networks_Data/regression_data_master_all.RData")

data_master <- data_master %>% unique() #only look at unique entries
data_master %>% mutate(RNA=RNA*1000000) -> data_master
###### group by gene/enhancer/promoter trips and normalise

data_master %>% group_by(GENEID,ENHANCERID,PROMID) %>% mutate(m.RNA=mean(RNA)) %>%
  mutate(v.RNA=var(RNA)) %>%
  mutate(n.RNA=ifelse(v.RNA==0,0,(RNA-m.RNA)/sqrt(v.RNA))) %>%
  dplyr::select(-m.RNA,-v.RNA)-> data_norm

data_norm %>% group_by(GENEID,ENHANCERID,PROMID) %>% mutate(m.TF=mean(BMAL1)) %>%
  mutate(v.TF=var(BMAL1)) %>%
  mutate(n.BMAL1=ifelse(v.TF==0,0,(BMAL1-m.TF)/sqrt(v.TF))) %>%
  dplyr::select(-m.TF,-v.TF) %>% group_by(GENEID,ENHANCERID,PROMID) %>% mutate(m.TF=mean(CLOCK)) %>%
  mutate(v.TF=var(CLOCK)) %>%
  mutate(n.CLOCK=ifelse(v.TF==0,0,(CLOCK-m.TF)/sqrt(v.TF))) %>%
  dplyr::select(-m.TF,-v.TF) %>% group_by(GENEID,ENHANCERID,PROMID) %>% mutate(m.TF=mean(CRY1)) %>%
  mutate(v.TF=var(CRY1)) %>%
  mutate(n.CRY1=ifelse(v.TF==0,0,(CRY1-m.TF)/sqrt(v.TF))) %>%
  dplyr::select(-m.TF,-v.TF) %>% group_by(GENEID,ENHANCERID,PROMID) %>% mutate(m.TF=mean(CRY2)) %>%
  mutate(v.TF=var(CRY2)) %>%
  mutate(n.CRY2=ifelse(v.TF==0,0,(CRY2-m.TF)/sqrt(v.TF))) %>%
  dplyr::select(-m.TF,-v.TF) %>% group_by(GENEID,ENHANCERID,PROMID) %>% mutate(m.TF=mean(PER1)) %>%
  mutate(v.TF=var(PER1)) %>%
  mutate(n.PER1=ifelse(v.TF==0,0,(PER1-m.TF)/sqrt(v.TF))) %>%
  dplyr::select(-m.TF,-v.TF)   %>% group_by(GENEID,ENHANCERID,PROMID) %>% mutate(m.TF=mean(PER2)) %>%
  mutate(v.TF=var(PER2)) %>%
  mutate(n.PER2=ifelse(v.TF==0,0,(PER2-m.TF)/sqrt(v.TF))) %>%
  dplyr::select(-m.TF,-v.TF)  %>% group_by(GENEID,ENHANCERID,PROMID) %>% mutate(m.TF=mean(H3K27ac)) %>%
  mutate(v.TF=var(H3K27ac)) %>%
  mutate(n.H3K27ac=ifelse(v.TF==0,0,(H3K27ac-m.TF)/sqrt(v.TF))) %>%
  dplyr::select(-m.TF,-v.TF) -> data_norm

save(data_norm,file="RNAChIPtogether_GRdist_norm.RData")
