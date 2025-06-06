library(edgeR)
library(MetaCycle)
library(tidyverse)
#library(DODR)
library(pracma)
library(ggforce)
library(patchwork)

setwd("C://Users/mqbssrmj/Documents/GRKO_analysis_rpkm/")
#setwd("~/GRKO_analysis/")

load("Data/featurecounts_out.RData")
readxl::read_excel("RM_timecourse_final_sample_list.xlsx") -> meta_data

# create groups -----------------------------------------------------------

group_genotype <- factor(meta_data$group_n)
group_time <- factor(meta_data$timepoint)
group_sex <- factor(meta_data$sex)

group<- paste(group_genotype,group_time,group_sex,sep="_")

exp_design <- data.frame(group=group_genotype[-c(49,50,51,52)],time=as.numeric(as.character(group_time[-c(49,50,51,52)]),rep=1))
exp_design$group %>% as.character() %>% as.factor -> exp_design$group
levels(exp_design$group)# the first level is the reference group
exp_design$group <- relevel(exp_design$group, "WT") # choose NC as the correct reference group

# y <- DGEList(counts=(fc$counts[,-c(51,52)]),group=group[-c(51,52)]) #remove controls
# cpm_counts <- cpm(y)
# 
# countCheck <- cpm_counts > 1
# keep <- which(rowSums(countCheck) >= 2)
# y <- y[keep,]
# cpm_counts <- cpm(y)

y <- DGEList(counts=fc$counts[,-c(49,50,51,52)],genes=fc$annotation,group=group[-c(49,50,51,52)])
calcNormFactors(y)-> y
y = estimateCommonDisp(y, verbose=TRUE)
rpkm_y<- rpkm(y)
FC_RNA_norm <- rpkm_y

dim(rpkm_y)
dim(exp_design)

#save(rpkm_y,exp_design,file = "~/JueunInternship/GR_data.Rdata")

raw_counts <- fc$counts[,-c(49,50,51,52)]

#save(raw_counts,exp_design,file = "~/JueunInternship/GR_data_rawcounts.Rdata")


## JTK Naive (double plot test - 5 ploting!!!)

# jtk_data_wt <- as.data.frame(FC_RNA_norm[,group_genotype[-c(51,52)]=="WT"])
# names(jtk_data_wt) <- paste0((names(jtk_data_wt) %>% strsplit("_") %>% unlist %>% .[(1:24)*4 - 3]),"_WT")
# jtk_data_wt %>%
#   mutate(ID=rownames(jtk_data_wt)) %>% 
#   dplyr::select(ID,everything()) -> jtk_data_wt
# timepoints <- as.numeric(as.character(group_time[group_genotype=="WT"]))
#timepoints <- timepoints1 + c(0,24,48,72,96)

# jtk_data_het <- as.data.frame(FC_RNA_norm[,group_genotype[-c(51,52)]=="HET"])
# jtk_data_het %>% .[,-c(25,26)] -> jtk_data_het
# names(jtk_data_het) <- paste0((names(jtk_data_het) %>% strsplit("_") %>% unlist %>% .[(1:24)*4 - 3]),"_CIA")
# jtk_data_het %>%
#   mutate(ID=rownames(jtk_data_het)) %>% 
#   dplyr::select(ID,everything()) -> jtk_data_het
# timepoints1 <- as.numeric(as.character(group_time[group_genotype=="HET"][-c(25,26)]))
# 
# jtk_data_wt_detrend <-t(detrend(t(jtk_data_wt[,-1]),tt="linear"))
# jtk_data_wt_detrend %>% as_tibble %>%
#   mutate(ID=jtk_data_wt$ID) %>% dplyr::select(ID,everything()) -> jtk_data_wt_detrend
# 
# jtk_data_het_detrend <-t(detrend(t(jtk_data_het[,-1]),tt="linear"))
# jtk_data_het_detrend %>% as_tibble %>%
#   mutate(ID=jtk_data_het$ID) %>% dplyr::select(ID,everything()) -> jtk_data_het_detrend

######
setwd("~/compareRhythms_extra/")

source("compareRhythms_model_select_diffrhy_2.R")
source("utils.R")

(FC_RNA_norm==0) %>% apply(1,sum) -> num_zeros 

(FC_RNA_norm[,group_genotype[-c(49,50,51,52)]=="WT"] == 0) %>% apply(1,sum) -> num_zeros_WT
(FC_RNA_norm[,group_genotype[-c(49,50,51,52)]=="HET"] == 0) %>% apply(1,sum) -> num_zeros_HET

(num_zeros==50) %>% sum()
((num_zeros_HET==24)| (num_zeros_WT==24)) %>% sum()

data_forCR <- FC_RNA_norm[!((num_zeros_HET==24)| (num_zeros_WT==24)),]

liver_out <- compareRhythms_model_select_diffrhy(data=data_forCR,
                                                 exp_design=exp_design,
                                                 schwarz_wt_cutoff = 0.6,
                                                 just_classify = FALSE,
                                                 diffrhy = TRUE)
# convert phases to 24hr clock
liver_out$results %>% as_tibble() %>%
  mutate(WT_phase=ifelse(WT_phase < 0, WT_phase /(2*pi)*24 + 24, WT_phase /(2*pi)*24)) %>%
  mutate(HET_phase=ifelse(HET_phase < 0, HET_phase /(2*pi)*24 + 24, HET_phase /(2*pi)*24)) %>%
  dplyr::select(id,cat_star,WT_phase,HET_phase) %>%
  pivot_longer(c(WT_phase,HET_phase),names_to = "group",values_to="phase") %>% 
  mutate(phase=ifelse(phase>23,phase-24,phase))-> liver_phases

##################

write.table(data_forCR[,exp_design$group=="WT"], file="data_forJTK.txt",
            sep="\t", quote=FALSE, row.names=TRUE,col.names = FALSE)

Meta_Output_WT <- meta2d(infile="data_forJTK.txt", filestyle="txt",
                      outputFile = FALSE, timepoints=seq(0, 46, by=2),
                      cycMethod=c("JTK"), outIntegration="both",minper = 24,maxper = 24,
                      adjustPhase="predictedPer")

Meta_Output_WT$JTK %>% mutate(LAG=ifelse(LAG>23,LAG-24,LAG)) -> Meta_Output_WT$JTK

Meta_Output_WT$JTK %>% filter(ADJ.P<0.05) %>% dim
Meta_Output_WT$JTK %>% filter(ADJ.P<0.05) %>% ggplot(aes(x=LAG)) + geom_histogram(fill="red",colour="black",breaks=seq(-1,23,by=2)) -> p1

liver_phases %>% filter(cat_star%in%c("same","loss")) %>% filter(group=="WT_phase") %>% dim
liver_phases %>% filter(cat_star%in%c("same","loss")) %>% filter(group=="WT_phase") %>% ggplot(aes(x=phase)) + geom_histogram(fill="red",colour="black",breaks=seq(-1,23,by=2)) -> p2

p1 + p2
(p1 + p2) * coord_polar(start=-2*pi/24) * scale_x_continuous(limits = c(-1,23),breaks=c(0,2,4,6,8,10,12,14,16,18,20,22)) * theme_minimal() * theme(axis.title = element_blank(),axis.text.y = element_blank())


p1+coord_polar()


merge(liver_phases,Meta_Output_WT$JTK,by.x="id",by.y="CycID") -> comb_data

comb_data %>%
  mutate(phase=floor(LAG/2)*2) %>%
  mutate(phase=ifelse(phase>23,phase-24,phase)) %>%
  group_by(cat_star,group,phase) %>%
  mutate(freq=n()) %>%
  ungroup %>%
  group_by(cat_star,group)%>%
  mutate(prop=freq/n()) %>%
  ungroup %>%
  dplyr::select(-id,-LAG,-BH.Q,-ADJ.P,-PER,-AMP) %>%
  unique %>%
  filter(cat_star!="arrhy") %>%
  mutate(cat_group=paste0(cat_star,group)) %>%
  filter(cat_group%in%c("sameWT_phase",
                        "lossWT_phase",
                        "changeWT_phase",
                        "changeHET_phase",
                        "gainHET_phase"))  -> phase_plot_data


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

comb_data %>% filter(cat_star%in%c("same")) %>% filter(group=="WT_phase") %>% ggplot(aes(x=phase)) + geom_histogram(fill=cbPalette[2] ,colour="black",breaks=seq(-1,23,by=2)) -> p1
comb_data %>% filter(cat_star%in%c("same")) %>% filter(group=="WT_phase") %>% ggplot(aes(x=LAG)) + geom_histogram(fill=cbPalette[2],colour="black",breaks=seq(-1,23,by=2)) -> p2

(p1 + p2) * coord_polar(start=-2*pi/24) * scale_x_continuous(limits = c(-1,23),breaks=c(0,2,4,6,8,10,12,14,16,18,20,22)) * theme_minimal() * theme(axis.title = element_blank(),axis.text.y = element_blank())

comb_data %>% filter(cat_star%in%c("loss")) %>% filter(group=="WT_phase") %>% ggplot(aes(x=phase)) + geom_histogram(fill=cbPalette[3],colour="black",breaks=seq(-1,23,by=2)) -> p1_l
comb_data %>% filter(cat_star%in%c("loss")) %>% filter(group=="WT_phase") %>% ggplot(aes(x=LAG)) + geom_histogram(fill=cbPalette[3],colour="black",breaks=seq(-1,23,by=2)) -> p2_l

(p1_l + p2_l) * coord_polar(start=-2*pi/24) * scale_x_continuous(limits = c(-1,23),breaks=c(0,2,4,6,8,10,12,14,16,18,20,22)) * theme_minimal() * theme(axis.title = element_blank(),axis.text.y = element_blank())

comb_data %>% filter(cat_star%in%c("gain")) %>% filter(group=="HET_phase") %>% ggplot(aes(x=phase)) + geom_histogram(fill=cbPalette[4],colour="black",breaks=seq(-1,23,by=2)) -> p1_g
comb_data %>% filter(cat_star%in%c("gain")) %>% filter(group=="HET_phase") %>% ggplot(aes(x=LAG)) + geom_histogram(fill=cbPalette[4],colour="black",breaks=seq(-1,23,by=2)) -> p2_g

(p1_g + p2_g) * coord_polar(start=-2*pi/24) * scale_x_continuous(limits = c(-1,23),breaks=c(0,2,4,6,8,10,12,14,16,18,20,22)) * theme_minimal() * theme(axis.title = element_blank(),axis.text.y = element_blank())


(p1 + p1_l + p1_g)* coord_polar(start=-2*pi/24) * scale_x_continuous(limits = c(-1,23),breaks=c(0,2,4,6,8,10,12,14,16,18,20,22)) * theme_minimal() * theme(axis.title = element_blank(),axis.text.y = element_blank())


phase_plot_data %>%
  filter(cat_group%in%c("lossWT_phase")) -> loss_vec
rep(loss_vec$phase,loss_vec$freq) %>% circular(units="hours") %>% rayleigh.test()

phase_plot_data %>%
  filter(cat_group%in%c("sameWT_phase")) -> same_vec
rep(same_vec$phase,same_vec$freq) %>% circular(units="hours") %>% rayleigh.test()

phase_plot_data %>%
  filter(cat_group%in%c("gainHET_phase")) -> gain_vec
rep(gain_vec$phase,gain_vec$freq) %>% circular(units="hours") %>% rayleigh.test()

(p2 + ggtitle("Same rhythm in both") + theme(plot.title = element_text(hjust = 0.5),plot.caption = element_text(size=12))+ labs(subtitle = "p < 0.01") +
    p2_l + ggtitle("Loss of rhythm") + theme(plot.title = element_text(hjust = 0.5),plot.caption = element_text(size=12))+ labs(subtitle = "p = 0.19")+ 
    p2_g+ ggtitle("Gain of rhythm") + theme(plot.title = element_text(hjust = 0.5),plot.caption = element_text(size=12))+ labs(subtitle = "p = 0.01"))*
  coord_polar(start=-2*pi/24) * scale_x_continuous(limits = c(-1,23),breaks=c(0,2,4,6,8,10,12,14,16,18,20,22)) * theme_minimal() * theme(axis.title = element_blank(),axis.text.y = element_blank()) -> p_phase


watson.two.test(rep(loss_vec$phase,loss_vec$freq) %>% circular(units="hours") ,
                rep(same_vec$phase,same_vec$freq) %>% circular(units="hours"))

watson.two.test(rep(loss_vec$phase,loss_vec$freq) %>% circular(units="hours") ,
                rep(gain_vec$phase,gain_vec$freq) %>% circular(units="hours"))

watson.two.test(rep(gain_vec$phase,gain_vec$freq) %>% circular(units="hours") ,
                rep(same_vec$phase,same_vec$freq) %>% circular(units="hours"))

ggplot()+
  geom_segment(aes(y = 1,yend= 2, x=2.5))+
  geom_segment(aes(y = 1,yend= 2, x=4.8))+
  geom_segment(aes(y = 1, x=2.5,xend=4.8)) + 
  geom_text(aes(label="p<0.01",x=3.65,y=0.7)) +
  #
  geom_segment(aes(y = 1,yend= 2, x=5.2))+
  geom_segment(aes(y = 1,yend= 2, x=7.5))+
  geom_segment(aes(y = 1, x=5.2,xend=7.5)) + 
  geom_text(aes(label="p<0.01",x=6.35,y=0.7)) +
  #
  geom_segment(aes(y = -.5,yend= .5, x=2.5))+
  geom_segment(aes(y = -.5,yend= .5, x=7.5))+
  geom_segment(aes(y = -.5, x=2.5,xend=7.5)) + 
  geom_text(aes(label="p<0.01",x=5,y=-.8)) +
  #
  theme_void() +
  xlim(c(2,8)) + ylim(c(-1.5,2)) -> p_pval

p_phase/p_pval + plot_layout(heights=c(2,1)) -> p_phase_p

ggsave("plots/p_phase_p.png",p_phase_p,height =3.5,width = 6.5)
