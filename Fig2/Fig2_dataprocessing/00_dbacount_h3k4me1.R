
library(DiffBind)
library(tidyverse)
setwd("/run/user/1003/gvfs/sftp:host=csf3.itservices.manchester.ac.uk,user=mqbssrmj/mnt/mr01-home01/mqbssrmj/rds-loudon-ray-rattray/rds-shared-loudon-ray-rattray/Takahashi_data/")

ddmm<-dba(sampleSheet="samplesheetH3K4me1_nocontrol_temp.csv")
consensusObj<-dba.peakset(ddmm,consensus=DBA_TREATMENT,minOverlap=0)
data.peakset1 <- dba.peakset(consensusObj, bRetrieve=TRUE)

ddmm_count<-dba.count(ddmm,data.peakset1,summits = T)

rbind(ddmm_count$peaks[[1]] %>% dplyr::select(Chr,Start,End,Summits) %>% mutate(Sample=1),
      ddmm_count$peaks[[2]] %>% dplyr::select(Chr,Start,End,Summits) %>% mutate(Sample=2),
      ddmm_count$peaks[[3]] %>% dplyr::select(Chr,Start,End,Summits) %>% mutate(Sample=3),
      ddmm_count$peaks[[4]] %>% dplyr::select(Chr,Start,End,Summits) %>% mutate(Sample=4),
      ddmm_count$peaks[[5]] %>% dplyr::select(Chr,Start,End,Summits) %>% mutate(Sample=5),
      ddmm_count$peaks[[6]] %>% dplyr::select(Chr,Start,End,Summits) %>% mutate(Sample=6)) %>%
  as_tibble()-> peaks

peaks %>% group_by(Chr,Start,End) %>% mutate(m.summit=ceiling(median(Summits))) %>%
  ungroup %>%
  mutate(n.Start=m.summit-500) %>%
  mutate(n.End=m.summit+500) %>%
  dplyr::select(seqnames=Chr,start=n.Start,end=n.End) %>%
  unique -> peaks
save(peaks,file="Enhancer_Counting/enhancer_peaklist_h3k4me1.RData")
