
load("Data/enh_prom_links_CHIC_3b.RData")
load("Data/enhancer_withRNA.RData")

enh_prom_links %>% mutate(PROMID=paste(prom_seqnames,prom_start,prom_end,sep=".")) %>% 
  merge(prom_regions%>% as.data.frame,by="PROMID") %>% 
  merge(.,RNA_n_jtk_2,by.y="To",by.x="GENEID") %>%
  mutate(ENHANCERID=paste(enh_seqnames,enh_start,enh_end,sep=".")) %>%
  merge(.,enh_withRNA,by.x="ENHANCERID",by.y="ENHANCERID") -> data_master_all

data_master_all %>% 
  mutate(JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
  dplyr::select(GENEID.y,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID.y,Rhythmic,Both_us,JTK_adjphase) %>% 
  as_tibble %>% group_by(GENEID.y) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(GENEID.y,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
  unique %>%
  filter(Rhythmic.g==TRUE) %>%
  filter(Both_us.g==TRUE) -> chic_TT

data_master_all %>% 
  mutate(JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
  dplyr::select(GENEID.y,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID.y,Rhythmic,Both_us,JTK_adjphase) %>% 
  as_tibble %>% group_by(GENEID.y) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(GENEID.y,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
  unique %>%
  filter(Rhythmic.g==TRUE) %>%
  filter(Both_us.g==FALSE) -> chic_TF

data_master_all %>% 
  mutate(JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
  dplyr::select(GENEID.y,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID.y,Rhythmic,Both_us,JTK_adjphase) %>% 
  as_tibble %>% group_by(GENEID.y) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(GENEID.y,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
  unique  -> chic_all

load("Data/paschodata_enhancer.RData")

data_master_all_unorm %>% 
  filter(prom_dist>1800) %>%
  dplyr::select(PROMID,To,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(To,Rhythmic,Both_us,JTK_adjphase) %>% 
  as_tibble %>% group_by(To) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(To,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
  unique %>%
  filter(Rhythmic.g==TRUE) %>%
  filter(Both_us.g==TRUE) -> nearest_TT

data_master_all_unorm %>% 
  filter(prom_dist>1800) %>%
  dplyr::select(PROMID,To,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(To,Rhythmic,Both_us,JTK_adjphase) %>% 
  as_tibble %>% group_by(To) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(To,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
  unique %>%
  filter(Rhythmic.g==TRUE) %>%
  filter(Both_us.g==FALSE) -> nearest_TF

data_master_all_unorm %>% 
  #filter(prom_dist>1800) %>%
  dplyr::select(PROMID,To,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(To,Rhythmic,Both_us,JTK_adjphase) %>% 
  as_tibble %>% group_by(To) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(To,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
  unique -> nearest_all

save(chic_TT,chic_TF,chic_all,nearest_TT,nearest_TF,nearest_all,file="Data/genelistsforcomparisonbetweenCHICandnearest.RData")

library(VennDiagram)
venn.plot <- draw.pairwise.venn(
  length(nearest_TT$To), 
  length(chic_TT$GENEID.y), 
  length(intersect(nearest_TT$To, chic_TT$GENEID.y)),
  category = c("nearest", "chic"),
  fill = c("blue", "red"),
  alpha = 0.5
)
venn.plot <- draw.pairwise.venn(
  length(nearest_TF$To), 
  length(chic_TF$GENEID.y), 
  length(intersect(nearest_TF$To, chic_TF$GENEID.y)),
  category = c("nearest", "chic"),
  fill = c("blue", "red"),
  alpha = 0.5
)

venn.plot <- draw.pairwise.venn(
  length(c(nearest_TF$To,nearest_TT$To)), 
  length(c(chic_TF$GENEID.y,chic_TT$GENEID.y)), 
  length(intersect(c(nearest_TF$To,nearest_TT$To), c(chic_TF$GENEID.y,chic_TT$GENEID.y))),
  category = c("nearest", "chic"),
  fill = c("blue", "red"),
  alpha = 0.5
)


venn.plot <- draw.pairwise.venn(
  length(c(nearest_all$To)), 
  length(c(chic_all$GENEID.y)), 
  length(intersect(c(nearest_all$To), c(chic_all$GENEID.y))),
  category = c("nearest", "chic"),
  fill = c("blue", "red"),
  alpha = 0.5
)

data_master_all_unorm %>% 
  filter(prom_dist>1800)  %>% 
  as_tibble %>% dplyr::select(To,ENHANCERID) %>% unique -> pairs_nn
data_master_all %>% as_tibble() %>% dplyr::select(GeneID,ENHANCERID) %>% unique() -> pairs_chic

intersect(pairs_nn$ENHANCERID,pairs_chic$ENHANCERID) -> ll

intersect((pairs_chic %>% filter(ENHANCERID %in% ll))$GeneID,(pairs_nn %>% filter(ENHANCERID %in% ll))$To) %>% unique %>% length
(pairs_nn %>% filter(ENHANCERID %in% ll))$To %>% unique %>% length
(pairs_chic %>% filter(ENHANCERID %in% ll))$GeneID %>% unique %>% length
venn.plot <- draw.pairwise.venn(
  1125, 
  1643, 
  408,
  category = c("nearest", "chic"),
  fill = c("blue", "red"),
  alpha = 0.5
)



(pairs_nn %>% mutate(tog=paste(To,ENHANCERID)) %>% filter(ENHANCERID %in% ll))$tog %>% unique %>% length
(pairs_chic %>% mutate(tog=paste(GeneID,ENHANCERID)) %>% filter(ENHANCERID %in% ll))$tog %>% unique %>% length
intersect((pairs_chic %>% mutate(tog=paste(GeneID,ENHANCERID)) %>% filter(ENHANCERID %in% ll))$tog,(pairs_nn %>% mutate(tog=paste(To,ENHANCERID)) %>% filter(ENHANCERID %in% ll))$tog) %>% unique %>% length
venn.plot <- draw.pairwise.venn(
  1327, 
  2670, 
  131,
  category = c("nearest", "chic"),
  fill = c("blue", "red"),
  alpha = 0.5
)



enh_prom_links %>% mutate(enhpromdist=Mod((prom_end-1100)-(enh_end-500))) -> ep_dists
data_master_all_unorm %>% separate(ENHANCERID,into=c("enh_chr","enh_start","enh_end"),remove = F) %>% separate(PROMID,into=c("prom_chr","prom_start","prom_end"),remove = F) %>% 
  mutate(enhpromdist=Mod((as.numeric(prom_end)-1100)-(as.numeric(enh_end)-500))) %>% dplyr::select(PROMID,ENHANCERID,enhpromdist) %>% unique() -> dm_dists


par(mfrow=c(1,2))
dm_dists$enhpromdist %>% log10 %>% hist(ylim = c(0,3500),xlim=c(0,9))
ep_dists$enhpromdist%>% log10 %>% hist(ylim = c(0,10000),xlim=c(0,9))


chic %>% mutate(enhpromdist=Mod((V3-V2)/2-(V6-V5)/2)) -> chic_dists

par(mfrow=c(3,1))
dm_dists$enhpromdist %>% log10 %>% hist(ylim = c(0,3500),xlim=c(0,9))
ep_dists$enhpromdist%>% log10 %>% hist(ylim = c(0,30000),xlim=c(0,9))
chic_dists$enhpromdist%>% log10 %>% hist(ylim = c(0,50000),xlim=c(0,9))


c(dm_dists$enhpromdist, ep_dists$enhpromdist)%>% log10 %>% hist(ylim = c(0,10000),xlim=c(0,9))



data_master_all %>% 
  dplyr::select(PROMID.y,GENEID.y,ENHANCERID) %>%
  mutate(ID=paste0(PROMID.y,GENEID.y,ENHANCERID)) %>%
  unique-> chic_all

data_master_all_unorm %>% 
  #filter(prom_dist>1800) %>%
  dplyr::select(PROMID,To,ENHANCERID) %>%
  mutate(ID=paste0(PROMID,To,ENHANCERID)) %>%
  unique -> nearest_all

venn.plot <- draw.pairwise.venn(
  length(c(nearest_all$To)), 
  length(c(chic_all$GENEID.y)), 
  length(intersect(c(nearest_all$ID), c(chic_all$ID))),
  category = c("nearest", "chic"),
  fill = c("blue", "red"),
  alpha = 0.5
)

data_master_all %>% 
  dplyr::select(GENEID.y,ENHANCERID) %>%
  mutate(ID=paste0(GENEID.y,ENHANCERID)) %>%
  unique-> chic_all

data_master_all_unorm %>% 
  dplyr::select(To,ENHANCERID) %>%
  mutate(ID=paste0(To,ENHANCERID)) %>%
  unique -> nearest_all

venn.plot <- draw.pairwise.venn(
  length(c(nearest_all$To)), 
  length(c(chic_all$GENEID.y)), 
  length(intersect(c(nearest_all$ID), c(chic_all$ID))),
  category = c("nearest", "chic"),
  fill = c("blue", "red"),
  alpha = 0.5
)


data_master_all %>% 
  dplyr::select(GENEID.y) %>%
  mutate(ID=paste0(GENEID.y)) %>%
  unique-> chic_all

data_master_all_unorm %>% 
  dplyr::select(To) %>%
  mutate(ID=paste0(To)) %>%
  unique -> nearest_all

venn.plot <- draw.pairwise.venn(
  length(c(nearest_all$To)), 
  length(c(chic_all$GENEID.y)), 
  length(intersect(c(nearest_all$ID), c(chic_all$ID))),
  category = c("nearest", "chic"),
  fill = c("blue", "red"),
  alpha = 0.5
)
