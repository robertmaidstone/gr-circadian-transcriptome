library(tidyverse)
library(patchwork)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

load("Data/enh_withnearestprom_data_reads.RData") #originally from NearestEnhancer folder

enh_withnearestprom_data %>%
  mutate(GR=(GR_lou_dex==TRUE)|(GR_lou_WT==TRUE)) -> data_m

data_m %>% dplyr::select(ENHANCERID)  %>%
  separate(ENHANCERID, c("seqnames", "starts", "ends"),sep="[.]") %>% unique -> enhancers

mm10 = TxDb.Mmusculus.UCSC.mm10.knownGene 

makeGRangesFromDataFrame(enhancers,seqnames.field = "seqnames",start.field = "starts",end.field = "ends") -> enh_granges

peakAnno <- annotatePeak(enh_granges, tssRegion=c(-2000, 400),
                         TxDb=mm10, annoDb="org.Mm.eg.db")

###############

peakAnno %>% as_tibble %>% mutate(ENHANCERID=paste(seqnames,start,end,sep=".")) %>%
  mutate(location=(str_split(annotation,pattern = " ") %>% purrr::map(.f = function(x){x[1]}) %>% unlist)) %>%
  dplyr::select(ENHANCERID,distanceToTSS,location) -> peakAnno_t

data_m %>% mutate(Prom=closest_prom_dist==0) %>%
  merge(peakAnno_t,by="ENHANCERID") -> data_master_m


data_master_m %>% 
  mutate(GR=c("No GR binding","GR binding")[GR+1]) %>%
  dplyr::select(ENHANCERID,Prom,distanceToTSS,GR,location) %>% 
  unique %>% 
  filter(Prom==F) %>% 
  filter(location!="Promoter")%>% 
  ggplot(aes(x=abs(distanceToTSS)/1000,fill=GR))+geom_density(alpha=0.5) + theme_bw() + scale_fill_manual(values = c("red","black")) +
  theme(legend.title = element_blank(),plot.title = element_text(size = 12)) +
  scale_x_log10(limits=c(0.1,2000),breaks=c(0.1,1,10,100,1000),labels=c("0.1","1","10","100","1000"))+ 
  xlab("Distance to TSS (kb)") + ylab("Density")+
  geom_text(label="p<0.01",x=2.9,y=.6)-> p_dist2TSS


data_master_m %>% 
  mutate(GR=c("No GR binding","GR binding")[GR+1]) %>%
  dplyr::select(ENHANCERID,Prom,location,distanceToTSS,GR) %>% 
  unique %>% 
  filter(Prom==F) %>% filter(location!="Promoter")%>% 
  dplyr::select(location,GR) %>%
  table %>% as_tibble %>%
  group_by(GR) %>%
  dplyr::mutate(freq=100*n/sum(n)) %>%
  ungroup %>%
  dplyr::mutate(location=factor(location,levels=unique(location[order(freq)]))) %>%
  ggplot(aes(x=location,fill=GR,y=freq))+
  geom_bar(alpha=0.5,color="black",position = "dodge",stat="identity") + 
  theme_bw() + scale_fill_manual(values = c("red","black")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.3,.85),
        legend.background = element_blank(),
        plot.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Annotated Genomic Location")+ylab("Percentage") -> p_barloc

p_barloc + p_dist2TSS + guides(fill=FALSE) -> loc_plots

ggsave(plot = loc_plots,"Fig2/Plots/locationofenhancers.png",width=5.5,height=3.5)

# t-test on log10 distance to TSS data
data_master_m %>% 
  mutate(GR=c("No GR binding","GR binding")[GR+1]) %>%
  dplyr::select(ENHANCERID,Prom,distanceToTSS,GR,location) %>% 
  unique %>% 
  filter(Prom==F) %>% 
  filter(location!="Promoter")%>% filter(GR=="No GR binding") %>% dplyr::select(distanceToTSS) %>% abs %>% log10 -> noGR_distancetoTSS

data_master_m %>% 
  mutate(GR=c("No GR binding","GR binding")[GR+1]) %>%
  dplyr::select(ENHANCERID,Prom,distanceToTSS,GR,location) %>% 
  unique %>% 
  filter(Prom==F) %>% 
  filter(location!="Promoter")%>% filter(GR=="GR binding") %>% dplyr::select(distanceToTSS) %>% abs %>% log10 -> GR_distancetoTSS

t.test(noGR_distancetoTSS,GR_distancetoTSS)
