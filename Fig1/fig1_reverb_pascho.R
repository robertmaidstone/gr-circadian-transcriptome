
library(tidyverse)
library(patchwork)

setwd("~/GRanalysis_master/Fig1/")
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
# venn  -------------------------------------------------------------------

# Gene 
data_master %>% 
  dplyr::select(PROMID,GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID,Rhythmic,Both_us) %>% 
  as_tibble %>% group_by(GENEID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(GENEID,Rhythmic.g,Both_us.g) %>%
  unique %>%
  dplyr::select(-GENEID) %>%
  table 

#Transcripts
data_master %>% 
  dplyr::select(PROMID,GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID,Rhythmic,Both_us) %>% 
  as_tibble %>% group_by(GENEID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(Rhythmic,Both_us) %>%
  table 

#Transcripts (merged promoters)
data_master %>% 
  dplyr::select(GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID,Rhythmic,Both_us) %>% 
  as_tibble %>% group_by(GENEID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(Rhythmic,Both_us) %>%
  table 

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

# data_master %>%
#   mutate(GR=(GRlou_dex_us==0)) %>%
#   group_by(GENEID,Time) %>%
#   dplyr::mutate(GR=any(GR)) %>%
#   mutate(GR=c("No GR binding","GR binding")[GR+1]) %>%
#   ungroup %>%
#   dplyr::select(GENEID,Time,GR,RNA,JTK_adjphase,JTK_pvalue) %>%
#   unique %>%
#   gather(TF,val,-GENEID,-Time,-GR,-JTK_adjphase,-JTK_pvalue) %>%
#   mutate(Time=as.numeric(as.character(Time))) %>%
#   filter(Time<24) %>%
#   group_by(GENEID,GR,TF) %>%
#   dplyr::mutate(m.val=mean(val)) %>%
#   ungroup %>%
#   dplyr::select(GENEID,GR,TF,m.val,JTK_adjphase,JTK_pvalue) %>%
#   unique -> plot_RNA_data

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

plot_RNA_data  %>% ggplot(aes(x=m.val+1,fill=GR)) + geom_density(alpha=0.5)  + theme_bw() + scale_fill_manual(values = c("red","black")) +
  theme(legend.title = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y = element_blank(),
        plot.title = element_text(size = 12)) +
  scale_x_log10(limits=c(1,10000))+ guides(fill=F) + ggtitle("RNA") +
  ylim(c(0,5)) + xlab("Reads") -> p_RNA


plot_RNA_data_2  %>% ggplot(aes(x=m.val+1,fill=GR)) + geom_density(alpha=0.5)  + theme_bw() + scale_fill_manual(values = c("red","black")) +
  theme(legend.title = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y = element_blank(),
        plot.title = element_text(size = 12)) +
  scale_x_log10(limits=c(1,10000))+ guides(fill="none") + ggtitle("RNA") +
  ylim(c(0,3.5)) + xlab("") -> p_RNA

data_master %>%
  mutate(GR=(GRlou_dex_us==0)) %>%
  mutate(GR=c("No GR binding","GR binding")[GR+1]) %>%
  dplyr::select(ENHANCERID,GENEID,PROMID,Time,BMAL1,CLOCK,CRY1,CRY2,NPAS2,PER1,PER2,GR,RNA) %>% 
  unique %>% 
  gather(TF,val,-ENHANCERID,-GENEID,-PROMID,-Time,-GR) %>%
  mutate(Time=as.numeric(as.character(Time))) %>%
  filter(Time<24) %>%
  group_by(ENHANCERID,GENEID,PROMID,GR,TF) %>%
  dplyr::mutate(m.val=mean(val)) %>%
  ungroup %>%
  dplyr::select(ENHANCERID,GENEID,PROMID,GR,TF,m.val) %>% 
  unique -> plot_data

# plot_TF(plot_data,"BMAL1",legend=F,x.axis = F)+
#   plot_TF(plot_data,"CLOCK",legend=F,x.axis = F,y.axis = F)+
#   plot_TF(plot_data,"CRY1",legend=F,x.axis = F,y.axis = F)+
#   plot_TF(plot_data,"CRY2",legend=F,x.axis = F,y.title = T)+
#   plot_TF(plot_data,"NPAS2",legend=F,x.axis = F,y.axis = F)+
#   plot_TF(plot_data,"PER1",legend=T,y.axis = F)+
#   plot_TF(plot_data,"PER2",legend=F) + p_RNA -> p_TFs
# 
# #ggsave(plot = p_TFs,"~/GR_Nov2019/TF_GRbinding_promoterV2.png",width=7,height=5)

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

ggsave(plot = p_TFs,"TF_GRbinding_promoter_reverb_pascho2.png",width=7,height=6)

data_master_rev %>%
  mutate(GR=(GRlou_dex_us==0)) %>%
  mutate(GR=c("No GR binding","GR binding")[GR+1]) %>%
  mutate(rhy=JTK_pvalue<0.05) %>%
  mutate(per_plot=ifelse(JTK_adjphase %in% c(22,0,2),"22-2",ifelse(JTK_adjphase %in% c(10,12,14),"10-14",NA))) %>%
  dplyr::select(ENHANCERID,GENEID,PROMID,Time,BMAL1,CLOCK,CRY1,CRY2,NPAS2,PER1,PER2,GR,RNA,`REVERBa*`=REVERBa,rhy,per_plot) %>% 
  unique %>% 
  gather(TF,val,-ENHANCERID,-GENEID,-PROMID,-Time,-GR,-rhy,-per_plot) %>%
  mutate(Time=as.numeric(as.character(Time))) %>%
  filter(Time<24) %>%
  group_by(ENHANCERID,GENEID,PROMID,GR,TF) %>%
  dplyr::mutate(m.val=mean(val,na.rm=T)) %>%
  ungroup %>%
  dplyr::select(ENHANCERID,GENEID,PROMID,GR,TF,m.val,rhy,per_plot) %>% 
  unique %>%
  filter(!is.na(per_plot)) %>%
  mutate(x_groups=paste(GR,per_plot)) %>%
  filter(rhy==T)-> plot_data

plot_RNA_data   %>% mutate(rhy=JTK_pvalue<0.05) %>%
  mutate(per_plot=ifelse(JTK_adjphase %in% c(22,0,2),"22-2",ifelse(JTK_adjphase %in% c(10,12,14),"10-14",NA))) %>%
  filter(!is.na(per_plot)) %>%
  mutate(x_groups=paste(GR,per_plot)) %>%
  filter(rhy==T) %>%
  ggplot(aes(x=x_groups,y=m.val+1,fill=GR)) + geom_boxplot(alpha=0.5) + 
  theme_bw() + scale_fill_manual(values = c("red","black")) +
  theme(legend.title = element_blank(),plot.title = element_text(size = 12),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_log10(limits=c(5,1000)) + ylab("") + xlab("")+ ggtitle("RNA") + guides(fill="none") -> p_RNA_bp

plot_TF_bp(plot_data,"BMAL1",legend=F,x.axis = F)+
  plot_TF_bp(plot_data,"CLOCK",legend=F,x.axis = F,y.axis = F)+
  plot_TF_bp(plot_data,"CRY1",legend=F,x.axis = F,y.axis = F)+
  plot_TF_bp(plot_data,"CRY2",legend=F,x.axis = F,y.title = T)+
  plot_TF_bp(plot_data,"NPAS2",legend=F,x.axis = F,y.axis = F)+
  plot_TF_bp(plot_data,"PER1",legend=T,x.axis = F,y.axis = F)+
  plot_TF_bp(plot_data,"PER2",legend=F) + plot_TF_bp(plot_data,"REVERBa*",legend=F,y.axis=F,x.title=T) + p_RNA_bp  -> p_TFs

ggsave(plot = p_TFs,"~/GRanalysis_master/Fig1/TFbinding_periodgroups.png",width=7,height=5)
# supp venn ---------------------------------------------------------------

data_master %>% 
  dplyr::select(PROMID,GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us) %>%
  unique %>%
  mutate(WT_us=(WT_us==0),Dex_us=(Dex_us==0)) %>% 
  dplyr::select(WT_us,Dex_us) %>% 
  table

data_master %>% 
  dplyr::select(PROMID,GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue) %>%
  unique %>%
  mutate(WT_us=(WT_us==0),Dex_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  group_by(WT_us,Dex_us) %>%
  mutate(pct.rhy=sum(Rhythmic)*100/length(Rhythmic)) %>%
  dplyr::select(pct.rhy) %>%
  unique

# Gene 
data_master %>% 
  dplyr::select(PROMID,GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue) %>%
  unique %>%
  mutate(WT_us=(WT_us==0),Dex_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID,WT_us,Dex_us) %>% 
  as_tibble %>% group_by(GENEID) %>% dplyr::mutate(WT_us.g=any(WT_us),Dex_us.g=any(Dex_us)) %>% 
  ungroup %>%
  dplyr::select(GENEID,WT_us.g,Dex_us.g) %>%
  unique %>%
  dplyr::select(-GENEID) %>%
  table 

data_master %>% 
  dplyr::select(PROMID,GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue) %>%
  unique %>%
  mutate(WT_us=(WT_us==0),Dex_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID,WT_us,Dex_us,Rhythmic) %>% 
  as_tibble %>% group_by(GENEID) %>% dplyr::mutate(WT_us.g=any(WT_us),Dex_us.g=any(Dex_us),Rhythmic.g=any(Rhythmic)) %>% 
  ungroup %>%
  dplyr::select(GENEID,WT_us.g,Dex_us.g,Rhythmic.g) %>%
  unique %>%
  group_by(WT_us.g,Dex_us.g) %>%
  dplyr::mutate(pct.rhy=sum(Rhythmic.g)*100/length(Rhythmic.g)) %>%
  dplyr::select(pct.rhy) %>%
  unique

#Transcripts
data_master %>% 
  dplyr::select(PROMID,GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue) %>%
  unique %>%
  mutate(WT_us=(WT_us==0),Dex_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(PROMID,GENEID,ENHANCERID,WT_us,Dex_us) %>% 
  as_tibble %>%
  #group_by(PROMID) %>% 
  group_by(PROMID,GENEID,ENHANCERID) %>% 
  dplyr::mutate(WT_us.g=any(WT_us),Dex_us.g=any(Dex_us)) %>% 
  ungroup %>%
  dplyr::select(PROMID,WT_us.g,Dex_us.g) %>%
  dplyr::select(-PROMID ) %>%
  table 

data_master %>% 
  dplyr::select(PROMID,GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue) %>%
  unique %>%
  mutate(WT_us=(WT_us==0),Dex_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID,WT_us,Dex_us,Rhythmic) %>% 
  as_tibble %>% group_by(GENEID) %>% dplyr::mutate(WT_us.g=any(WT_us),Dex_us.g=any(Dex_us),Rhythmic.g=any(Rhythmic)) %>% 
  ungroup %>%
  dplyr::select(GENEID,WT_us.g,Dex_us.g,Rhythmic.g) %>%
  group_by(WT_us.g,Dex_us.g) %>%
  dplyr::mutate(pct.rhy=sum(Rhythmic.g)*100/length(Rhythmic.g)) %>%
  dplyr::select(pct.rhy) %>%
  unique


#Transcripts (merged promoters)
data_master %>% 
  dplyr::select(GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue) %>%
  unique %>%
  mutate(WT_us=(WT_us==0),Dex_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID,WT_us,Dex_us) %>% 
  as_tibble %>% group_by(GENEID) %>% dplyr::mutate(WT_us.g=any(WT_us),Dex_us.g=any(Dex_us)) %>% 
  ungroup %>%
  dplyr::select(GENEID,WT_us.g,Dex_us.g) %>%
  dplyr::select(-GENEID) %>%
  table 

data_master %>% 
  dplyr::select(GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue) %>%
  unique %>%
  mutate(WT_us=(WT_us==0),Dex_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID,WT_us,Dex_us,Rhythmic) %>% 
  as_tibble %>% group_by(GENEID) %>% dplyr::mutate(WT_us.g=any(WT_us),Dex_us.g=any(Dex_us),Rhythmic.g=any(Rhythmic)) %>% 
  ungroup %>%
  dplyr::select(GENEID,WT_us.g,Dex_us.g,Rhythmic.g) %>%
  group_by(WT_us.g,Dex_us.g) %>%
  dplyr::mutate(pct.rhy=sum(Rhythmic.g)*100/length(Rhythmic.g)) %>%
  dplyr::select(pct.rhy) %>%
  unique


# phase plots -------------------------------------------------------------

(data_master %>% 
   dplyr::select(PROMID,GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
   unique %>%
   mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
   dplyr::select(GENEID,Rhythmic,Both_us,JTK_adjphase) %>% 
   as_tibble %>% group_by(GENEID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
   ungroup %>%
   dplyr::select(GENEID,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
   unique %>%
   filter(Rhythmic.g==TRUE) %>%
   filter(Both_us.g==TRUE) %>%
   ggplot(aes(x=JTK_adjphase))+
   geom_bar(color="black",fill="red",size=0.25,width=2) +
   theme_bw()+
   theme(panel.border=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
   coord_polar(start=-pi/12) +ggtitle("GR bound and Rhythmic")+
   scale_x_continuous(breaks=c(0,3,6,9,12,15,18,21))) +
  
  (data_master %>% 
     dplyr::select(PROMID,GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
     unique %>%
     mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
     dplyr::select(GENEID,Rhythmic,Both_us,JTK_adjphase) %>% 
     as_tibble %>% group_by(GENEID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
     ungroup %>%
     dplyr::select(GENEID,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
     unique %>%
     filter(Rhythmic.g==TRUE) %>%
     filter(Both_us.g==FALSE) %>%
     ggplot(aes(x=JTK_adjphase))+
     geom_bar(color="black",fill="red",size=0.25,width=2) +
     theme_bw()+
     theme(panel.border=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
     coord_polar(start=-pi/12)+ggtitle("Not GR bound and Rhythmic")) + 
    scale_x_continuous(breaks=c(0,3,6,9,12,15,18,21))-> pplots

ggsave(plot = pplots,"~/GR_Nov2019/RNAphase_promoterV2.png",width=5,height=3)


data_master %>% 
  dplyr::select(PROMID,GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID,Rhythmic,Both_us,JTK_adjphase) %>% 
  as_tibble %>% group_by(GENEID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(GENEID,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
  unique %>%
  filter(Rhythmic.g==TRUE) %>%
  filter(Both_us.g==FALSE) %>%
  dplyr::select(GENEID) %>% unlist -> rhythmic_noGR_promoter
  
data_master %>% filter(GENEID %in% rhythmic_noGR_promoter)  %>%
  as_tibble %>%
  dplyr::select(ENHANCERID) %>%
  unique %>%
  separate(ENHANCERID,into=c("seqnames","starts","ends")) %>%
  write.table(file="csv/rhythmic_noGR_promoter.bed", quote=F, sep="\t", row.names=F, col.names=F)

data_master %>% 
  dplyr::select(PROMID,GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID,Rhythmic,Both_us,JTK_adjphase) %>% 
  as_tibble %>% group_by(GENEID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(GENEID,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
  unique %>%
  filter(Rhythmic.g==TRUE) %>%
  filter(Both_us.g==TRUE) %>%
  dplyr::select(GENEID) %>% unlist -> rhythmic_GR_promoter

data_master %>% filter(GENEID %in% rhythmic_GR_promoter)  %>%
  as_tibble %>%
  dplyr::select(ENHANCERID) %>%
  unique %>%
  separate(ENHANCERID,into=c("seqnames","starts","ends")) %>%
  write.table(file="csv/rhythmic_GR_promoter.bed", quote=F, sep="\t", row.names=F, col.names=F)


data_master %>% 
  dplyr::select(PROMID,GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID,Rhythmic,Both_us,JTK_adjphase) %>% 
  as_tibble %>% group_by(GENEID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(GENEID,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
  unique %>%
  filter(Rhythmic.g==FALSE) %>%
  filter(Both_us.g==TRUE) %>%
  dplyr::select(GENEID) %>% unlist -> norhythmic_GR_promoter

data_master %>% filter(GENEID %in% norhythmic_GR_promoter)  %>%
  as_tibble %>%
  dplyr::select(ENHANCERID) %>%
  unique %>%
  separate(ENHANCERID,into=c("seqnames","starts","ends"),sep = "[.]") %>%
  write.table(file="csv/norhythmic_GR_promoter.bed", quote=F, sep="\t", row.names=F, col.names=F)
