

library(mclust)
library(tidyverse)
library(patchwork)
library(pracma)

load(file="Data/RNAChIPtogether_GRdist.RData")
load("Fig2/Data/paschodata.RData")

data_master_all_unorm %>% 
  dplyr::select(PROMID,To,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,RNA1:RNA24) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  as_tibble %>% group_by(To) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  pivot_longer(RNA1:RNA24,names_to = "Sample",values_to = "RNA") %>%
  mutate(Time=factor(Sample,levels=c("RNA1","RNA2","RNA3","RNA4","RNA5",
                                       "RNA6","RNA7","RNA8","RNA9","RNA10",
                                       "RNA11","RNA12","RNA13","RNA14","RNA15",
                                       "RNA16","RNA17","RNA18","RNA19","RNA20",
                                       "RNA21","RNA22","RNA23","RNA24"),labels=c(rep(seq(0, 20, by=4),each=4)))) -> ll
ll %>%
  dplyr::select(To,Rhythmic.g,Both_us.g,Time,RNA) %>%
  unique %>% 
  group_by(To,Time) %>%
  mutate(RNA.m=mean(RNA)) %>%
  dplyr::select(-RNA) %>% 
  unique %>%
  ungroup %>%
  spread(Time,RNA.m) -> rna_data

((rna_data %>% dplyr::select(-To,-Rhythmic.g,-Both_us.g)) !=0 ) %>%
 apply(1,sum) -> num_nonzeros

filtered_dat <- rna_data[num_nonzeros>2,]

filtered_dat %>%
  dplyr::select(-To,-Rhythmic.g,-Both_us.g) %>%
  t() %>%
  scale %>%
  t() %>%
  as_tibble -> m_dat



set.seed(102)
Mclust(m_dat[!is.na(m_dat[,1]),] ,modelNames = "EII",G = 10,verbose = FALSE) -> tt

tt$classification %>% table

cbind(filtered_dat[!is.na(m_dat[,1]),],Cluster=tt$classification) %>%  as_tibble %>%
  group_by(Cluster) %>%
  dplyr::mutate(num_clus=length(Cluster)) %>%
  ungroup %>%
  dplyr::mutate(Cluster=factor(Cluster,levels=unique(Cluster[order(num_clus)]))) %>%
  dplyr::select(-num_clus)-> clus_data

clus_data %>%
  gather(key,value,-Cluster,-To,-Rhythmic.g,-Both_us.g) %>%
  #mutate(Cluster=as.character(Cluster)) %>%
  mutate(key=as.numeric(as.character(key))) %>%
  group_by(key,Cluster,Both_us.g) %>%
  dplyr::mutate(m.val=median(value)) %>%
  ungroup %>%
  dplyr::select(key,Cluster,m.val,Both_us.g) %>% 
  unique -> clus_plot_data

clus_plot_data %>%
  filter(Both_us.g==F) %>%
  ggplot(aes(x=key,y=m.val,colour=Cluster)) +geom_line(size=1) +
  theme_bw()+ggtitle("No GR Binding")+ 
  ylim(c(0,1500))+ 
  ylab("Median RPKM") +xlab("Time")+
  theme(legend.position = "none")+
  scale_x_continuous(breaks=c(0,8,16,24,32,40))-> clust_noGR

clus_plot_data %>%
  filter(Both_us.g==T) %>%
  ggplot(aes(x=key,y=m.val,colour=Cluster)) +geom_line(size=1) +
  theme_bw()+ggtitle("GR Bound")+ 
  ylim(c(0,1500))+ 
  xlab("Time")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 5)))+
  scale_x_continuous(breaks=c(0,8,16,24,32,40))-> clust_GR

clus_data%>%
  dplyr::select(Both_us.g,Cluster) %>% table %>% as.data.frame() %>% 
  group_by(Both_us.g) %>% dplyr::mutate(pct.freq=100*Freq/sum(Freq)) %>%
  ungroup %>%
  dplyr::mutate(Both_us.g=factor(Both_us.g,labels=c("No GR Binding","GR bound"))) %>%
  ggplot(aes(y=pct.freq,x=Cluster,fill=Both_us.g))+geom_bar(position="dodge",stat="identity",alpha=0.5) + 
  theme_bw() + ylab("Percentage of genes")+
  theme(legend.title = element_blank(),
        legend.position = c(0.3,.9),
        legend.background = element_blank())+
  scale_fill_manual(values=c("black", "red", "#56B4E9"))-> bar1

clus_data%>%
  dplyr::select(Both_us.g,Cluster) %>% table %>% as.data.frame() %>% 
  group_by(Both_us.g) %>% dplyr::mutate(pct.freq=100*Freq/sum(Freq)) %>%
  ungroup %>%
  dplyr::mutate(Both_us.g=factor(Both_us.g,labels=c("No GR Binding","GR bound"))) %>%
  dplyr::select(Cluster,Freq) %>%
  mutate(Freq=1) %>% unique %>%
  ggplot(aes(x=Cluster,fill=Cluster,y=Freq)) + geom_bar(position="dodge",stat="identity") +
  theme_void() + 
  guides(fill="none") -> legend_bar

layout <- "
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
BBBCCDD
"

bar1 + legend_bar + clust_noGR + clust_GR+ plot_layout(design = layout) -> clust_plot
ggsave(plot = clust_plot,filename = "clust_plot_fitz.png",width=10,height=5)

################


filtered_dat %>%
  dplyr::select(`0`,`4`,`8`,`12`,`16`,`20`) %>%
  t %>%
  detrend(tt = "linear") %>%
  t %>%
  as_tibble -> filtered_dat_detrend

filtered_dat_detrend %>%
  t %>%
  scale %>%
  t %>% as_tibble -> m_dat_detrend

set.seed(106)
Mclust(m_dat_detrend[!is.na(m_dat_detrend[,1]),] ,modelNames = "EII",G = 10,verbose = FALSE) -> tt2

tt2$classification %>% table

cbind(filtered_dat[!is.na(m_dat_detrend[,1]),],Cluster=tt2$classification) %>%  as_tibble %>%
  group_by(Cluster) %>%
  dplyr::mutate(num_clus=length(Cluster)) %>%
  ungroup %>%
  dplyr::mutate(Cluster=factor(Cluster,levels=unique(Cluster[order(num_clus)]))) %>%
  dplyr::select(-num_clus)-> clus_data

clus_data %>%
  gather(key,value,-Cluster,-To,-Rhythmic.g,-Both_us.g) %>%
  #mutate(Cluster=as.character(Cluster)) %>%
  mutate(key=as.numeric(as.character(key))) %>%
  group_by(key,Cluster,Both_us.g) %>%
  dplyr::mutate(m.val=median(value)) %>%
  ungroup %>%
  dplyr::select(key,Cluster,m.val,Both_us.g) %>% 
  unique -> clus_plot_data

clus_plot_data %>%
  filter(Both_us.g==F) %>%
  ggplot(aes(x=key,y=m.val,colour=Cluster)) +geom_line(size=1) +
  theme_bw()+ggtitle("No GR Binding")+ 
  ylim(c(0,1500))+
  ylab("Median RPKM") +xlab("Time")+
  theme(legend.position = "none")+
  scale_x_continuous(breaks=c(0,8,16,24,32,40))-> clust_noGR

clus_plot_data %>%
  filter(Both_us.g==T) %>%
  ggplot(aes(x=key,y=m.val,colour=Cluster)) +geom_line(size=1) +
  theme_bw()+ggtitle("GR Bound")+ 
  ylim(c(0,1500))+
  xlab("Time")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 5)))+
  scale_x_continuous(breaks=c(0,8,16,24,32,40))-> clust_GR

clus_data%>%
  dplyr::select(Both_us.g,Cluster) %>% table %>% as.data.frame() %>% 
  group_by(Both_us.g) %>% dplyr::mutate(pct.freq=100*Freq/sum(Freq)) %>%
  ungroup %>%
  dplyr::mutate(Both_us.g=factor(Both_us.g,labels=c("No GR Binding","GR bound"))) %>%
  ggplot(aes(y=pct.freq,x=Cluster,fill=Both_us.g))+geom_bar(position="dodge",stat="identity",alpha=0.5) + 
  theme_bw() + ylab("Percentage of genes")+
  theme(legend.title = element_blank(),
        legend.position = c(0.3,.9),
        legend.background = element_blank())+
  scale_fill_manual(values=c("black", "red", "#56B4E9"))-> bar1

clus_data%>%
  dplyr::select(Both_us.g,Cluster) %>% table %>% as.data.frame() %>% 
  group_by(Both_us.g) %>% dplyr::mutate(pct.freq=100*Freq/sum(Freq)) %>%
  ungroup %>%
  dplyr::mutate(Both_us.g=factor(Both_us.g,labels=c("No GR Binding","GR bound"))) %>%
  dplyr::select(Cluster,Freq) %>%
  mutate(Freq=1) %>% unique %>%
  ggplot(aes(x=Cluster,fill=Cluster,y=Freq)) + geom_bar(position="dodge",stat="identity") +
  theme_void() + 
  guides(fill="none") -> legend_bar

layout <- "
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
AAACCDD
BBBCCDD
"

bar1 + legend_bar + clust_noGR + clust_GR+ plot_layout(design = layout) -> clust_detrend
ggsave(plot = clust_detrend,filename = "clust_detrend_fitz.png",width=10,height=5)
  