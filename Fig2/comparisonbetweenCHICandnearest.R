
load("Data/enh_prom_links_CHIC.RData")
load("Data/enhancer_withRNA.RData")

enh_prom_links %>% mutate(PROMID=paste(prom_seqnames,prom_start,prom_end,sep=".")) %>% 
  merge(TSS_ourproms%>% mutate(PROMID=paste(Chr,Start,End,sep=".")),by="PROMID") %>% 
  merge(.,RNA_n_jtk_2,by.y="To",by.x="GeneID") %>%
  mutate(ENHANCERID=paste(enh_seqnames,enh_start,enh_end,sep=".")) %>%
  merge(.,enh_withRNA,by.x="ENHANCERID",by.y="ENHANCERID") -> data_master_all

data_master_all %>% 
  mutate(JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
  dplyr::select(GeneID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GeneID,Rhythmic,Both_us,JTK_adjphase) %>% 
  as_tibble %>% group_by(GeneID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(GeneID,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
  unique %>%
  filter(Rhythmic.g==TRUE) %>%
  filter(Both_us.g==TRUE) -> chic_TT

data_master_all %>% 
  mutate(JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
  dplyr::select(GeneID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GeneID,Rhythmic,Both_us,JTK_adjphase) %>% 
  as_tibble %>% group_by(GeneID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(GeneID,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
  unique %>%
  filter(Rhythmic.g==TRUE) %>%
  filter(Both_us.g==FALSE) -> chic_TF

data_master_all %>% 
  mutate(JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
  dplyr::select(GeneID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GeneID,Rhythmic,Both_us,JTK_adjphase) %>% 
  as_tibble %>% group_by(GeneID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(GeneID,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
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
  filter(prom_dist>1800) %>%
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
  length(chic_TT$GeneID), 
  length(intersect(nearest_TT$To, chic_TT$GeneID)),
  category = c("nearest", "chic"),
  fill = c("blue", "red"),
  alpha = 0.5
)
venn.plot <- draw.pairwise.venn(
  length(nearest_TF$To), 
  length(chic_TF$GeneID), 
  length(intersect(nearest_TF$To, chic_TF$GeneID)),
  category = c("nearest", "chic"),
  fill = c("blue", "red"),
  alpha = 0.5
)

venn.plot <- draw.pairwise.venn(
  length(c(nearest_TF$To,nearest_TT$To)), 
  length(c(chic_TF$GeneID,chic_TT$GeneID)), 
  length(intersect(c(nearest_TF$To,nearest_TT$To), c(chic_TF$GeneID,chic_TT$GeneID))),
  category = c("nearest", "chic"),
  fill = c("blue", "red"),
  alpha = 0.5
)


venn.plot <- draw.pairwise.venn(
  length(c(nearest_all$To)), 
  length(c(chic_all$GeneID)), 
  length(intersect(c(nearest_all$To), c(chic_all$GeneID))),
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

