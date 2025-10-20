data_model %>%
  dplyr::mutate(Rhy=Rhythmic.g) %>%
  dplyr::mutate(GR=Both_us.g) %>%
  dplyr::select(GeneID,RNA=RNA_m,GR,Rhy) %>%
  group_by(GeneID) %>%
  dplyr::mutate(GR_any=any(GR)) %>%
  dplyr::mutate(GR_any=ifelse(is.na(GR_any),FALSE,GR_any)) %>%
  dplyr::mutate(m.RNA=mean(RNA)) %>% 
  dplyr::mutate(GRrhy=ifelse(GR_any&Rhy,"Rhythmic with \nGR binding at \nlinked enhancer",ifelse(GR_any,"Non rhythmic with \nGR binding at \nlinked enhancer",
                                                                                                 ifelse(Rhy,"Rhythmic with \nno GR binding at \nlinked enhancer","Non rhythmic with \nno GR binding at \nlinked enhancer")))) %>%
  dplyr::select(-RNA,-GR) %>%
  ungroup %>%
  unique %>%
  mutate(GR_any=c("No GR binding", "GR binding")[GR_any+1]) -> RNA_pCHiC


data_model %>%
  dplyr::mutate(Rhy=Rhythmic.g) %>%
  dplyr::mutate(GR=Both_us.g) %>%
  dplyr::select(GENEID,RNA=RNA_m,GR,Rhy) %>%
  group_by(GENEID) %>%
  dplyr::mutate(GR_any=any(GR)) %>%
  dplyr::mutate(GR_any=ifelse(is.na(GR_any),FALSE,GR_any)) %>%
  dplyr::mutate(m.RNA=mean(RNA)) %>% 
  dplyr::mutate(GRrhy=ifelse(GR_any&Rhy,"Rhythmic with \nGR binding at \nlinked enhancer",ifelse(GR_any,"Non rhythmic with \nGR binding at \nlinked enhancer",
                                                                                                 ifelse(Rhy,"Rhythmic with \nno GR binding at \nlinked enhancer","Non rhythmic with \nno GR binding at \nlinked enhancer")))) %>%
  dplyr::select(-RNA,-GR) %>%
  ungroup %>%
  unique %>%
  mutate(GR_any=c("No GR binding", "GR binding")[GR_any+1]) %>%
  dplyr::mutate(GeneID=as.character(GENEID)) %>% dplyr::select(-GENEID)-> RNA_nearest


rbind(cbind(RNA_nearest,link="Nearest"),
cbind(RNA_pCHiC,link="pC-HiC")) %>% 
  ggplot(aes(x=m.RNA,y=link)) +
  geom_jitter(colour="grey",alpha=.5,height=.2,shape=16) +
  geom_boxplot(fill=NA) +
  #scale_x_log10(name="RPKM") +
  # scale_y_continuous(name="Density")+
  #scale_colour_manual(values = c("red","black"))+
  ylab("") +
  theme_bw() +
  theme(legend.position = c(0.75, .92),legend.background = element_blank()) +
  guides(colour="none")  + scale_x_log10("Normalised Counts") -> p_compare

ggsave(p_compare,filename = "Fig2/Plots/rna_boxplots_compare.png",height = 3,width = 4)
