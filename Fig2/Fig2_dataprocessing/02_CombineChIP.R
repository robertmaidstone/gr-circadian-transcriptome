library(tidyverse)
library(edgeR)
library(pracma)
setwd("Data")


RSS_fit<-c()
RSS_test<-c()
RSS_fit_int<-c()
RSS_test_int<-c()
coeffs_sep<-c()

for(timepoint in c(0,4,8,12,16,20)){
  load(paste0("CT",timepoint,"data_nc.RData"))
  
  READS_filt<-eval(parse(text=paste0("READS_CT",timepoint,"_filt")) )
  
  ############## filter out   blacklisted #############
  
  BlackRegions<- read.csv("BL_enh.bed",sep="\t",header = F)
  
  !(paste(READS_filt$Chr,READS_filt$Start,READS_filt$End) %in% paste(BlackRegions$V1,BlackRegions$V2,BlackRegions$V3)) -> bl_logic
  READS_filt %>%
    filter(bl_logic) -> READS_filt

  ####################################################################
  ################# Convert histone marks to CPM ##################### from macs files (total tags)
  c(4744198,6970058,4492392,5671641,4828132,11796353) -> h3k27ac_filt_tags
  #c(24299934,15842142,17336065,19634619,23683261,18270194) -> h3k4me3_filt_tags
  c(19727506,20780588,18628639,18220262,21368998,20959840) -> h3k4me1_filt_tags
  
  READS_filt$H3K27ac<-(READS_filt$H3K27ac*1000000)/h3k27ac_filt_tags[which(timepoint== c(0,4,8,12,16,20))]
  READS_filt$H3K4me1<-(READS_filt$H3K4me1*1000000)/h3k4me1_filt_tags[which(timepoint== c(0,4,8,12,16,20))]
  ####################################################################
  
 # cbind(READS_filt %>% dplyr::select(ID,Chr,Start,End,maxBF) , log(READS_filt %>% dplyr::select(-ID,-Chr,-Start,-End,-maxBF)) ) %>%
 #   as_tibble() ->READS_filt
  
  assign(paste0("READS_CT",timepoint,"_filt"),READS_filt)
  
}



READS_CT0_filt
READS_CT0_filt %>% 
  dplyr::rename(seqnames=Chr,starts=Start,ends=End) %>%
  dplyr::select(c(2,3,4)) %>%
  write.table(file="ChIP_regions.bed", quote=F, sep="\t", row.names=F, col.names=F)

##########################################################

rbind(cbind(READS_CT0_filt,Time=0),
      cbind(READS_CT4_filt,Time=4),
      cbind(READS_CT8_filt,Time=8),
      cbind(READS_CT12_filt,Time=12),
      cbind(READS_CT16_filt,Time=16),
      cbind(READS_CT20_filt,Time=20),
      cbind(READS_CT0_filt,Time=24),
      cbind(READS_CT4_filt,Time=28),
      cbind(READS_CT8_filt,Time=32),
      cbind(READS_CT12_filt,Time=36),
      cbind(READS_CT16_filt,Time=40),
      cbind(READS_CT20_filt,Time=44)) %>% mutate(ENHANCERID=paste(Chr,Start,End,sep=".")) -> READS_temp

READS_temp  %>% dplyr::select(-ID,-Chr,-Start,-End) -> data_master


head(data_master)
dim(data_master)
dim(data_master)[1]/12

data_master <- data_master %>% unique
head(data_master)
dim(data_master)
dim(data_master)[1]/12


# h3k27ac JTK -------------------------------------------------------------

data_master %>% as_tibble %>% dplyr::select(ENHANCERID,Time,H3K27ac) %>%
  filter(Time < 24) %>%
  spread(Time,H3K27ac) -> h3k27ac_data

h3k27ac_data_detrend <- as_tibble(cbind(ENHANCERID=h3k27ac_data$ENHANCERID,t(detrend(t(h3k27ac_data[,-1]),tt="linear"))))

##### JTK Cycle analysis of RNA counts (after normalisation)

library(MetaCycle)
write.table(h3k27ac_data_detrend, file="h3k27ac_data_detrend_forJTK.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

Meta_Output <- meta2d(infile="h3k27ac_data_detrend_forJTK.txt", filestyle="txt",
                      outputFile = FALSE, timepoints=seq(0, 20, by=4),
                      cycMethod=c("JTK"), outIntegration="both",outRawData = F)

merge(data_master,Meta_Output$JTK,by.x="ENHANCERID",by.y="CycID")-> dm
dm %>% as_tibble -> data_master


save(data_master, file="ChIPtogether.RData")


