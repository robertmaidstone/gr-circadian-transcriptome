load("Data/enh_prom_links_CHIC_3b.RData")
load("Data/enhancer_withRNA.RData")

enh_prom_links %>% mutate(PROMID=paste(prom_seqnames,prom_start,prom_end,sep=".")) %>% 
  merge(prom_regions%>% as.data.frame,by="PROMID") %>% 
  merge(.,RNA_n_jtk_2,by.y="To",by.x="GENEID") %>%
  mutate(ENHANCERID=paste(enh_seqnames,enh_start,enh_end,sep=".")) %>%
  merge(.,enh_withRNA,by.x="ENHANCERID",by.y="ENHANCERID") -> data_master_all

data_master_all %>% mutate(enhpromdist=Mod((prom_end-1100)-(enh_end-500))) %>% 
  dplyr::select(GENEID.y,enhpromdist) %>%
  group_by(GENEID.y) %>%
  mutate(min_dist=min(enhpromdist)) %>%
  dplyr::select(GENEID.y,min_dist) %>%
  unique -> ep_dists

load("Data/paschodata_enhancer.RData")

data_master_all_unorm %>% separate(ENHANCERID,into=c("enh_chr","enh_start","enh_end"),remove = F) %>%
  separate(PROMID,into=c("prom_chr","prom_start","prom_end"),remove = F) %>% 
  mutate(enhpromdist=Mod((as.numeric(prom_end)-1100)-(as.numeric(enh_end)-500))) %>% 
  dplyr::select(To,enhpromdist) %>% unique() %>%
  group_by(To) %>%
  mutate(min_dist=min(enhpromdist)) %>%
  dplyr::select(To,min_dist) %>% unique() -> dm_dists


par(mfrow=c(1,2))
dm_dists$min_dist %>% log10 %>% hist(ylim = c(0,2500),xlim=c(0,9))
ep_dists$min_dist%>% log10 %>% hist(ylim = c(0,2500),xlim=c(0,9))

rbind(cbind(dm_dists,link="Nearest"),
cbind(ep_dists,link="pC-HiC")) %>%
  filter(link=="Nearest") %>%
  ggplot(aes(x=min_dist)) + geom_histogram(fill="grey",colour="black") +
  scale_x_log10(limits=c(0.1,10000000),breaks=c(10,1000,100000),labels = label_number())+
  xlab("")+ ylab("Frequency") + 
  ggtitle("Nearest") + ylim(c(0,1500)) +theme_bw() -> p1

rbind(cbind(dm_dists,link="Nearest"),
      cbind(ep_dists,link="pC-HiC")) %>%
  filter(link=="pC-HiC") %>%
  ggplot(aes(x=min_dist)) + geom_histogram(fill="grey",colour="black") +
  scale_x_log10(limits=c(0.1,10000000),breaks=c(10,1000,100000),labels = label_number())+
  xlab("")+ ylab("Frequency") + 
  ggtitle("pC-HiC") + ylim(c(0,1500)) +theme_bw() -> p2

p1 +p2

png("Fig2/Plots/distance_compare.png", width = 1800, height = 900, res = 300)  # adjust size/res as needed
#grid.newpage()
grid.draw(
  patchworkGrob(p1 | p2)
)
grid.text("Distance to closest linked enhancer (bp)", y = unit(0.03, "npc"), gp = gpar(fontsize = 10))
dev.off()

###

enh_prom_links %>% mutate(PROMID=paste(prom_seqnames,prom_start,prom_end,sep=".")) %>% 
  merge(prom_regions%>% as.data.frame,by="PROMID") %>% 
  merge(.,RNA_n_jtk_2,by.y="To",by.x="GENEID") %>%
  mutate(ENHANCERID=paste(enh_seqnames,enh_start,enh_end,sep=".")) %>%
  merge(.,enh_withRNA,by.x="ENHANCERID",by.y="ENHANCERID") -> data_master_all

data_master_all %>% mutate(enhpromdist=Mod((prom_end-1100)-(enh_end-500))) %>% 
  dplyr::select(GENEID.y,enhpromdist) %>%
  group_by(GENEID.y) %>%
  mutate(min_dist=min(enhpromdist)) %>%
  dplyr::select(GENEID.y,min_dist) %>%
  unique -> ep_dists

data_master_all %>% 
  dplyr::select(GENEID.y) %>%
  mutate(ID=paste0(GENEID.y)) %>%
  unique-> chic_all

data_master_all_unorm %>% 
  dplyr::select(To) %>%
  mutate(ID=paste0(To)) %>%
  unique -> nearest_all


png("Fig2/Plots/geneven_compare.png", width = 900, height = 900, res = 300)  # adjust size/res as needed
venn.plot <- draw.pairwise.venn(
  length(c(nearest_all$To)), 
  length(c(chic_all$GENEID.y)), 
  length(intersect(c(nearest_all$ID), c(chic_all$ID))),
  category = c("", ""),
  fill = c(NA, NA),
  alpha = 0.5
)
dev.off()
