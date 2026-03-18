
library(tidyverse)
library(patchwork)

setwd("~/GRanalysis_master/GR_CircadianLiverTranscriptome/")
load(file="Data/RNAChIPtogether_GRdist.RData")
load("Data/promoter_counts_reverb.RData")

# add reverb data ---------------------------------------------------------

ddmm_count$peaks[[1]] %>% as_tibble %>%
  mutate(ENHANCERID=paste(Chr,Start,End,sep=".")) %>%
  dplyr::select(ENHANCERID,Sample1=Reads) -> rev_sam1

ddmm_count$peaks[[2]] %>% as_tibble %>%
  mutate(ENHANCERID=paste(Chr,Start,End,sep=".")) %>%
  dplyr::select(ENHANCERID,Sample2=Reads) -> rev_sam2

rev_data <- rev_sam1 %>% mutate(Sample2=rev_sam2$Sample2) %>%
  mutate(REVERBa=(Sample1+Sample2)/2) %>%
  dplyr::select(-Sample1,-Sample2)

merge(data_master,rev_data,by="ENHANCERID") %>% mutate(REVERBa=ifelse(Time==8,REVERBa,NA)) %>% as_tibble -> data_master_rev

save(data_master_rev,file = "Data/promoter_data_final.RData")

# Density plots -----------------------------------------------------------

#for transcripts (merged promoters)

plot_TF<-function(plot_data,TF,x.axis=T,y.axis=T,legend=T,y.title=F,x.title=F){
  temp<-TF
  (plot_data %>% filter(TF==temp) %>% ggplot(aes(x=m.val+1,fill=GR)) + geom_density(alpha=0.5) + ylim(c(0,3.5)) + theme_bw() + scale_fill_manual(values = c("red","black","green","blue")) +
     theme(legend.title = element_blank(),plot.title = element_text(size = 12)) +
      scale_x_log10(limits=c(1,10000))+ xlab("Reads") + ylab("Density")+ ggtitle(TF)) -> p1
  if(legend==FALSE){p1+guides(fill=FALSE) -> p1 }
  if(x.axis==FALSE){p1+ theme(axis.text.x=element_blank(),
                                  axis.ticks.x=element_blank()) -> p1 }
  if(y.axis==FALSE){p1+ theme(axis.text.y=element_blank(),
                              axis.ticks.y=element_blank()) -> p1 }
  if(y.title==FALSE){p1+ theme(axis.title.y=element_blank()) -> p1 }
  if(x.title==FALSE){p1+ theme(axis.title.x=element_blank()) -> p1 }
  return(p1)
}

plot_TF_bp<-function(plot_data,TF,x.axis=T,y.axis=T,legend=T,y.title=F,x.title=F){
  temp<-TF
  (plot_data %>% filter(TF==temp) %>% ggplot(aes(x=x_groups,y=m.val+1,fill=GR)) + geom_boxplot(alpha=0.5) + 
      theme_bw() + scale_fill_manual(values = c("red","black")) +
      theme(legend.title = element_blank(),plot.title = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      scale_y_log10(limits=c(5,1000)) + ylab("Reads") + xlab("Period")+ ggtitle(TF)) -> p1
  if(legend==FALSE){p1+guides(fill=FALSE) -> p1 }
  if(x.axis==FALSE){p1+ theme(axis.text.x=element_blank(),
                              axis.ticks.x=element_blank()) -> p1 }
  if(y.axis==FALSE){p1+ theme(axis.text.y=element_blank(),
                              axis.ticks.y=element_blank()) -> p1 }
  if(y.title==FALSE){p1+ theme(axis.title.y=element_blank()) -> p1 }
  if(x.title==FALSE){p1+ theme(axis.title.x=element_blank()) -> p1 }
  return(p1)
}

load("Data/paschodata.RData")

data_master_all_unorm %>%
  mutate(GENEID=To) %>%
  mutate(GR=(GRlou_dex_us==0)) %>%
  group_by(GENEID,Time) %>%
  dplyr::mutate(GR=any(GR)) %>%
  mutate(GR=c("No GR binding","GR binding")[GR+1]) %>%
  ungroup %>%
  dplyr::select(GENEID,GR,RNA1:RNA24) %>%
  unique %>%
  gather(RNA_sample,val,-GENEID,-GR) %>%
  group_by(GENEID,GR) %>%
  dplyr::mutate(m.val=mean(val)) %>%
  ungroup %>%
  dplyr::select(GENEID,GR,m.val) %>%
  unique -> plot_RNA_data_2

plot_RNA_data_2  %>% ggplot(aes(x=m.val+1,fill=GR)) + geom_density(alpha=0.5)  + theme_bw() + scale_fill_manual(values = c("red","black")) +
  theme(legend.title = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y = element_blank(),
        plot.title = element_text(size = 12)) +
  scale_x_log10(limits=c(1,10000))+ guides(fill="none") + ggtitle("RNA") +
  ylim(c(0,3.5)) + xlab("") -> p_RNA

data_master_rev %>%
  mutate(GR=(GRlou_dex_us==0)) %>%
  mutate(GR=c("No GR binding","GR binding")[GR+1]) %>%
  dplyr::select(ENHANCERID,GENEID,PROMID,Time,BMAL1,CLOCK,CRY1,CRY2,NPAS2,PER1,PER2,GR,RNA,`REVERBa*`=REVERBa) %>% 
  unique %>% 
  gather(TF,val,-ENHANCERID,-GENEID,-PROMID,-Time,-GR) %>%
  mutate(Time=as.numeric(as.character(Time))) %>%
  filter(Time<24) %>%
  group_by(ENHANCERID,GENEID,PROMID,GR,TF) %>%
  dplyr::mutate(m.val=mean(val,na.rm=T)) %>%
  ungroup %>%
  dplyr::select(ENHANCERID,GENEID,PROMID,GR,TF,m.val) %>% 
  unique -> plot_data

plot_TF(plot_data,"BMAL1",legend=F,x.axis = F)+
  plot_TF(plot_data,"CLOCK",legend=F,x.axis = F,y.axis = F)+
  plot_TF(plot_data,"CRY1",legend=T,x.axis = F,y.axis = F)+ theme(legend.position = "inside",legend.position.inside = c(.6, .75)) +
  plot_TF(plot_data,"CRY2",legend=F,x.axis = F,y.title = T)+
  plot_TF(plot_data,"NPAS2",legend=F,x.axis = F,y.axis = F)+
  plot_TF(plot_data,"PER1",legend=F,y.axis = F,x.axis=F)+
  plot_TF(plot_data,"PER2",legend=F) + 
  plot_TF(plot_data,"REVERBa*",legend=F,y.axis=F,x.title = "Reads")+ggtitle("REVERBa")  + 
  p_RNA -> p_TFs

ggsave(plot = p_TFs,"Fig1/Plots/TF_GRbinding_promoter_reverb_pascho2.png",width=7,height=6)
 ###
