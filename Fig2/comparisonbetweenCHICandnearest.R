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

data_master_all %>% 
  mutate(JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
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

save(chic_TT,chic_TF,nearest_TT,nearest_TF,file="Data/genelistsforcomparisonbetweenCHICandnearest.RData")

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
