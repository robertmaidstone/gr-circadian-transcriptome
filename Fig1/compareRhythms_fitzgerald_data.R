

#compareRhythms phases

cbind(RNA_n[,],RNA_n[,-1]) -> RNA_CR

group_genotype <- factor(c(rep("WT",24),rep("HET",24)))
group_time=c(rep(seq(0, 20, by=4),each=4),rep(seq(0, 20, by=4),each=4))

group<- paste(group_genotype,group_time,sep="_")

exp_design <- data.frame(group=group_genotype,time=group_time)
levels(exp_design$group)# the first level is the reference group
exp_design$group <- relevel(exp_design$group, "WT") # choose NC as the correct reference group

######
setwd("~/compareRhythms_extra/")

source("compareRhythms_model_select_diffrhy_2.R")
source("utils.R")

(RNA_CR==0) %>% apply(1,sum) -> num_zeros 

(RNA_CR[,group_genotype=="WT"] == 0) %>% apply(1,sum) -> num_zeros_WT
(RNA_CR[,group_genotype=="HET"] == 0) %>% apply(1,sum) -> num_zeros_HET

(num_zeros==48) %>% sum()
((num_zeros_HET==24)| (num_zeros_WT==24)) %>% sum()

data_forCR <- RNA_CR[!((num_zeros_HET==24)| (num_zeros_WT==24)),]


data_forCR[,-1]<-mean(count_sums) *data_forCR[,-1]

#data_forCR <- FC_RNA_norm[!((num_zeros_HET==24)| (num_zeros_WT==24)),]

CR_out <- compareRhythms_model_select_diffrhy(data=data_forCR,
                                                 exp_design=exp_design,
                                                 schwarz_wt_cutoff = 0.6,
                                                 just_classify = FALSE,
                                                 diffrhy = TRUE)
table(CR_out$results$category) #with diffrhy category
table(CR_out$results$cat_star) #with diffrhy category split and added to loss, gain and change


# convert phases to 24hr clock
CR_out$results %>% as_tibble() %>%
  mutate(WT_phase=ifelse(WT_phase < 0, WT_phase /(2*pi)*24 + 24, WT_phase /(2*pi)*24)) %>%
  mutate(WT_phase=ifelse(WT_phase>23,WT_phase-24,WT_phase)) %>%
  dplyr::select(id,category,WT_phase) -> CR_phases

##plotting

CR_phases %>% filter(category=="same") %>%# filter(group=="WT_phase") %>%
  mutate(WT_phase=ifelse(WT_phase>23,WT_phase-24,WT_phase)) %>%
  ggplot(aes(x=WT_phase)) + geom_bar(color="black",fill="grey",size=0.25,width=1) +
  theme_bw()+
  theme(panel.border=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
  #ylim(0,505)  +
  coord_polar(start=-pi/12) +
  scale_x_binned(breaks= c(-1,1,3,5,7,9,11,13,15,17,19,21,23))

read.csv("Data/conv_pascho.txt",sep = "\t") ->conv_1
tibble(geneID=data_forCR[,1],CRID=CR_out$mppa$geneID) %>% merge(CR_phases,by.x="CRID",by.y="id") %>%
  merge(conv_1 %>% dplyr::select(From,To),by.x = "geneID",by.y="From",all.x = T)  -> ll


load("Data/RNAChIPtogether_GRdist.RData")

data_master %>% head
merge(ll,data_master,by.x="To",by.y="GENEID",all=TRUE) -> data_master_all



(data_master_all %>% 
    dplyr::select(PROMID,To,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,category,WT_phase) %>%
    unique %>%
    mutate(Both_us=(Dex_us==0),Rhythmic=(category=="same")) %>% 
    dplyr::select(To,Rhythmic,Both_us,WT_phase) %>% 
    as_tibble %>% group_by(To) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
    ungroup %>%
    dplyr::select(To,Rhythmic.g,Both_us.g,WT_phase) %>%
    unique %>%
    filter(Rhythmic.g==TRUE) %>%
    filter(Both_us.g==TRUE) %>%
    mutate(WT_phase=2*round(WT_phase/2)) %>%
    ggplot(aes(x=WT_phase))+
    geom_bar(color="black",fill="red",size=0.25,width=2) +
    theme_bw()+
    theme(panel.border=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
    coord_polar(start=-pi/12) +ggtitle("GR bound and Rhythmic")+
    scale_x_continuous(breaks=c(0,3,6,9,12,15,18,21))+ labs(subtitle = "p < 0.01")) +
  (data_master_all %>% 
     dplyr::select(PROMID,To,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,category,WT_phase) %>%
     unique %>%
     mutate(Both_us=(Dex_us==0),Rhythmic=(category=="same")) %>% 
     dplyr::select(To,Rhythmic,Both_us,WT_phase) %>% 
     as_tibble %>% group_by(To) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
     ungroup %>%
     dplyr::select(To,Rhythmic.g,Both_us.g,WT_phase) %>%
     unique %>%
     filter(Rhythmic.g==TRUE) %>%
     filter(Both_us.g==FALSE) %>%
     mutate(WT_phase=2*round(WT_phase/2)) %>%
     ggplot(aes(x=WT_phase))+
     geom_bar(color="black",fill="red",size=0.25,width=2) +
     theme_bw()+
     theme(panel.border=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
     coord_polar(start=-pi/12) +ggtitle("GR bound and Rhythmic")+
     scale_x_continuous(breaks=c(0,3,6,9,12,15,18,21))+ labs(subtitle = "p < 0.01")) -> pplots

#
data_master_all %>% 
  dplyr::select(PROMID,To,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,category,WT_phase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=(category=="same")) %>% 
  dplyr::select(To,Rhythmic,Both_us,WT_phase) %>% 
  as_tibble %>% group_by(To) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(To,Rhythmic.g,Both_us.g,WT_phase) %>%
  unique -> phase_plot_data

phase_plot_data %>%
  filter(Rhythmic.g==TRUE) %>%
  filter(Both_us.g==FALSE) -> rhynoGR_vec
rhynoGR_vec$WT_phase %>% circular(units="hours") %>% rayleigh.test()

phase_plot_data %>%
  filter(Rhythmic.g==TRUE) %>%
  filter(Both_us.g==TRUE) -> rhyyesGR_vec
rhyyesGR_vec$WT_phase %>% circular(units="hours") %>% rayleigh.test()

watson.two.test(rhyyesGR_vec$WT_phase %>% circular(units="hours"),
                rhynoGR_vec$WT_phase %>% circular(units="hours"))

ggplot()+
  geom_segment(aes(y = 1,yend= 2, x=2.5))+
  geom_segment(aes(y = 1,yend= 2, x=4.5))+
  geom_segment(aes(y = 1, x=2.5,xend=4.5)) + 
  geom_text(aes(label="p<0.01",x=3.5,y=0.7)) + #from watson test
  theme_void() +
  xlim(c(1.7,5.3)) + ylim(c(-1.5,2)) -> p_pval

pplots/p_pval + plot_layout(heights=c(2,.6)) -> p_phase_p

ggsave("Fig1/Plots/p_phase_p_CR.png",p_phase_p,height =4,width = 6.5)

# Venn diagram

data_master_all %>% 
  dplyr::select(PROMID,To,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,category,WT_phase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=(category=="same")) %>% 
  dplyr::select(To,Rhythmic,Both_us) %>% 
  as_tibble %>% group_by(To) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(To,Rhythmic.g,Both_us.g) %>%
  unique %>%
  dplyr::select(-To) -> venn_data

library(VennDiagram) 
venn.diagram(list(`GR Bound` = which(venn_data$Both_us.g),Rhythmic = which(venn_data$Rhythmic.g)), 
             fill = c("white", "white"), 
             category.names = c("" , "" ),
             alpha = c(0.5, 0.5),
             lwd =1,
             cat.cex=0.6,
             disable.logging = TRUE,
             height=1500,width=1500, "Fig1/Plots/venn_diagram.png")
