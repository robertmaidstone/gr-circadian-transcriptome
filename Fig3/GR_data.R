library(edgeR)
library(MetaCycle)
library(tidyverse)
library(DODR)
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

source("compareRhythms_model_select_diffrhy.R")
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
table(liver_out$results$category) #with diffrhy category
table(liver_out$results$cat_star) #with diffrhy category split and added to loss, gain and change

# convert phases to 24hr clock
liver_out$results %>% as_tibble() %>%
  mutate(WT_phase=ifelse(WT_phase < 0, WT_phase /(2*pi)*24 + 24, WT_phase /(2*pi)*24)) %>%
  mutate(HET_phase=ifelse(HET_phase < 0, HET_phase /(2*pi)*24 + 24, HET_phase /(2*pi)*24)) %>%
  dplyr::select(id,cat_star,WT_phase,HET_phase) %>%
  pivot_longer(c(WT_phase,HET_phase),names_to = "group",values_to="phase") -> liver_phases

##plotting

liver_phases %>% filter(cat_star=="same") %>% filter(group=="WT_phase") %>%
  mutate(phase=ifelse(phase>23,phase-24,phase)) %>%
  ggplot(aes(x=phase)) + geom_bar(color="black",fill="grey",size=0.25,width=1) +
  theme_bw()+
  theme(panel.border=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
  #ylim(0,505)  + 
  coord_polar(start=-pi/12) + 
  scale_x_binned(breaks= c(-1,1,3,5,7,9,11,13,15,17,19,21,23)) 

liver_phases %>% filter(cat_star=="change") %>% 
  mutate(phase=ifelse(phase>23,phase-24,phase)) %>%
  ggplot(aes(x=phase)) + geom_bar(color="black",fill="grey",size=0.25,width=1) +
  theme_bw()+
  theme(panel.border=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
  #ylim(0,505)  + 
  coord_polar(start=-pi/12) + 
  scale_x_binned(breaks= c(-1,1,3,5,7,9,11,13,15,17,19,21,23)) + 
  facet_grid(~group)

liver_phases %>% filter(cat_star=="gain") %>% filter(group=="HET_phase") %>%
  mutate(phase=ifelse(phase>23,phase-24,phase)) %>%
  ggplot(aes(x=phase)) + geom_bar(color="black",fill="grey",size=0.25,width=1) +
  theme_bw()+
  theme(panel.border=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
  #ylim(0,505)  + 
  coord_polar(start=-pi/12) + 
  scale_x_binned(breaks= c(-1,1,3,5,7,9,11,13,15,17,19,21,23)) 

liver_phases %>% filter(cat_star=="loss") %>% filter(group=="WT_phase") %>%
  mutate(phase=ifelse(phase>23,phase-24,phase)) %>%
  ggplot(aes(x=phase)) + geom_bar(color="black",fill="grey",size=0.25,width=1) +
  theme_bw()+
  theme(panel.border=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
  #ylim(0,505)  + 
  coord_polar(start=-pi/12) + 
  scale_x_binned(breaks= c(-1,1,3,5,7,9,11,13,15,17,19,21,23))


liver_phases %>%  
  separate(group,into="group",extra = "drop") %>%
  mutate(category=paste(cat_star,group,sep="_")) %>%
  filter(category%in%c("same_WT","change_WT","loss_WT","change_HET","gain_HET")) %>%
  mutate(category=ifelse(category=="same_WT","same",category)) %>%
  mutate(phase=ifelse(phase>23,phase-24,phase)) %>%
  ggplot(aes(x=phase,colour=category,fill=category)) + geom_point(stat="count") + geom_polygon(stat="count",alpha=0.2) +
  theme_bw()+
  theme(panel.border=element_blank(),axis.title.x=element_blank(),
        #axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank()
  )+
  #ylim(0,505)  + 
  coord_polar(start=-pi/12) + 
  scale_x_binned(breaks= c(-1,1,3,5,7,9,11,13,15,17,19,21,23)) -> p1

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

liver_phases %>%  
  separate(group,into="group",extra = "drop") %>%
  mutate(category=paste(cat_star,group,sep="_")) %>%
  filter(category%in%c("same_WT","change_WT","loss_WT","change_HET","gain_HET")) %>%
  mutate(phase=ifelse(phase>23,phase-24,phase)) %>%
  filter(category!="same_WT") %>%
  ggplot(aes(x=phase,colour=category,fill=category)) + geom_point(stat="count") + geom_polygon(stat="count",alpha=0.2) +
  theme_bw()+
  theme(panel.border=element_blank(),axis.title.x=element_blank(),
        #axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank()
  )+
  #ylim(0,505)  + 
  coord_polar(start=-pi/12) + 
  scale_x_binned(breaks= c(-1,1,3,5,7,9,11,13,15,17,19,21,23)) + 
  scale_fill_manual(values=gg_color_hue(5)[1:4]) + 
  scale_colour_manual(values=gg_color_hue(5)[1:4])-> p2

library(patchwork)
p1 + p2




######

#FC_RNA_norm[apply(FC_RNA_norm,1,sum)!=0,] -> FC_RNA_norm_filt
FC_RNA_norm[!((num_zeros_HET==24)| (num_zeros_WT==24)),] -> FC_RNA_norm_filt

liver_out$mppa

merge(liver_out$mppa,liver_out$results,by.x="geneID",by.y="id",all.x=TRUE) -> res_merge

res_merge %>% 
  mutate(cat_star=as.character(cat_star)) %>%
  mutate(cat_star=ifelse(cat_star=="gain","gain*",cat_star)) %>% 
  mutate(cat_star=ifelse(cat_star=="change","change*",cat_star)) %>% 
  mutate(cat_star=ifelse(cat_star=="loss","loss*",cat_star)) %>%
  mutate(WT_phase=ifelse(WT_phase < 0, WT_phase/(2*pi)*24 + 24, WT_phase/(2*pi)*24)) %>%
  mutate(HET_phase=ifelse(HET_phase < 0, HET_phase/(2*pi)*24 + 24, HET_phase/(2*pi)*24)) -> res_merge

res_merge$category %>% table
res_merge$cat_star %>% table

plot_heat<-function(rhy_cat,lag_ord,lims=c(-2,2)){
  timepoints <- seq(0,48,by=2)
  FC_RNA_norm_filt %>% as_tibble(rownames = "geneID") %>% 
    merge(res_merge,by="geneID") %>% as_tibble %>% filter(cat_star==rhy_cat)%>% 
    gather(Sample,Value,colnames(FC_RNA_norm_filt)[group_genotype[-c(49,50,51,52)]=="WT"]) %>% dplyr::mutate(Time=factor(Sample,ordered = T,levels=unique(Sample))) %>%
    dplyr::mutate(Time=timepoints[as.numeric(Time)]) %>% group_by(geneID,Time) %>%
    dplyr::mutate(m.Value=mean(Value)) %>% as_tibble %>%
    dplyr::select(geneID,Time,m.Value,HET_phase,WT_phase) %>%
    unique %>%
    group_by(geneID) %>%
    dplyr::mutate(m_all=mean(m.Value)) %>%
    dplyr::mutate(var_all=var(m.Value)) %>%
    ungroup %>%
    dplyr::mutate(n.Value=(m.Value-m_all)/sqrt(var_all)) %>%
    dplyr::mutate(n.Value=ifelse(n.Value < lims[1],lims[1],ifelse(n.Value > lims[2],lims[2],n.Value)))-> mm
  
  if(lag_ord=="WT"){
    mm %>% dplyr::select(geneID,LAG=WT_phase) %>% unique -> mm_ordering
  }else{
    mm %>% dplyr::select(geneID,LAG=HET_phase) %>% unique -> mm_ordering
  }
  mm %>% dplyr::mutate(geneID=factor(geneID,ordered=T,levels=mm_ordering$geneID[order(mm_ordering$LAG,decreasing = T)])) %>%
    ggplot(aes(x=Time,y=geneID,fill=n.Value))+geom_tile() +scale_fill_gradient2( low = "blue", high = "red", mid="white",limits=lims) + theme_void() ->  rhy_cia_n
  
  ########## het 
  FC_RNA_norm_filt %>% as_tibble(rownames = "geneID") %>% 
    merge(res_merge,by="geneID") %>% as_tibble %>% filter(cat_star==rhy_cat)%>% 
    gather(Sample,Value,colnames(FC_RNA_norm_filt)[group_genotype[-c(49,50,51,52)]=="HET"]) %>% dplyr::mutate(Time=factor(Sample,ordered = T,levels=unique(Sample))) %>%
    dplyr::mutate(Time=timepoints[as.numeric(Time)]) %>% group_by(geneID,Time) %>%
    dplyr::mutate(m.Value=mean(Value)) %>% as_tibble %>%
    dplyr::select(geneID,Time,m.Value,HET_phase,WT_phase) %>%
    unique %>%
    group_by(geneID) %>%
    dplyr::mutate(m_all=mean(m.Value)) %>%
    dplyr::mutate(var_all=var(m.Value)) %>%
    ungroup %>%
    dplyr::mutate(n.Value=(m.Value-m_all)/sqrt(var_all))%>%
    dplyr::mutate(n.Value=ifelse(n.Value < lims[1],lims[1],ifelse(n.Value > lims[2],lims[2],n.Value)))-> mm
  if(lag_ord=="WT"){
    mm %>% dplyr::select(geneID,LAG=WT_phase) %>% unique -> mm_ordering
  }else{
    mm %>% dplyr::select(geneID,LAG=HET_phase) %>% unique -> mm_ordering
  }
  mm %>% dplyr::mutate(geneID=factor(geneID,ordered=T,levels=mm_ordering$geneID[order(mm_ordering$LAG,decreasing = T)])) %>%
    ggplot(aes(x=Time,y=geneID,fill=n.Value))+geom_tile() +scale_fill_gradient2( low = "blue", high = "red", mid="white",limits=lims) + theme_void() -> rhy_cia_cia
  
  return(rhy_cia_n + guides(fill=FALSE)  + rhy_cia_cia)
}

plot_heat("gain*","HET",c(-2,2)) 
plot_heat("loss*","WT",c(-2,2))
plot_heat("same","WT",c(-2,2))
plot_heat("change*","WT",c(-2,2))

width_per_gene <- 0.005
#width_per_gene <- 0.08

#width_per_gene <- 0.015
p <- plot_heat("gain*","HET",c(-2,2))+ guides(fill="none")
ggsave("GRplots/gainstar.png",p,height = width_per_gene*(res_merge %>% filter(cat_star=="gain*") %>% dim %>% .[1]),width = 4)
p <- plot_heat("loss*","WT",c(-2,2))+ guides(fill="none")
ggsave("GRplots/lossstar.png",p,height = width_per_gene*(res_merge %>% filter(cat_star=="loss*") %>% dim %>% .[1]),width = 4)
#width_per_gene <- 0.0035
p <- plot_heat("same","HET",c(-2,2)) + guides(fill="none")
ggsave("GRplots/same.png",p,height = width_per_gene*(res_merge %>% filter(cat_star=="same") %>% dim %>% .[1]),width = 4)

width_per_gene <- 0.02
p <- plot_heat("change*","WT",c(-2,2))+ guides(fill="none")
ggsave("GRplots/changestar_WT.png",p,height = width_per_gene*(res_merge %>% filter(cat_star=="change*") %>% dim %>% .[1]),width = 4)
p <- plot_heat("change*","HET",c(-2,2))+ guides(fill="none")
ggsave("GRplots/changestar_HET.png",p,height = width_per_gene*(res_merge %>% filter(cat_star=="change*") %>% dim %>% .[1]),width = 4)


res_merge %>% dplyr::select(geneID,cat_star) %>%
  mutate(cat_star=ifelse(is.na(cat_star),"unassigned",cat_star)) %>%
  merge(FC_RNA_norm %>% as_tibble(rownames = "geneID") %>% dplyr::select(geneID),by="geneID",all=T) %>%
  mutate(cat_star=ifelse(is.na(cat_star),"filtered",cat_star)) %>%
  mutate(cat_star=factor(cat_star,levels=c("arrhy","same","gain*","loss*","change*","unassigned","filtered"),
                                     labels=c("Arrhythmic\n in both","Same rhythm\n in both","Gain of rhythm\n in GR-LKO","Loss of rhythm\n in GR-LKO","Altered rhythm\n in GR-LKO","Unassigned category\n by compareRhythms","Filtered due\n to low counts"))) %>%
  ggplot(aes(x=cat_star)) + geom_bar(fill="grey",colour="black") + geom_text(stat='count', aes(label=..count..), vjust=-1)+ theme_bw() +
  theme(axis.title = element_blank())

res_merge %>% dplyr::select(geneID,cat_star) %>%
  mutate(cat_star=ifelse(is.na(cat_star),"unassigned",cat_star)) %>%
  merge(FC_RNA_norm %>% as_tibble(rownames = "geneID") %>% dplyr::select(geneID),by="geneID",all=T) %>%
  mutate(cat_star=ifelse(is.na(cat_star),"filtered",cat_star)) %>%
  mutate(cat_star=factor(cat_star,levels=c("arrhy","same","gain*","loss*","change*","filtered","unassigned"),
                         labels=c("Arrhythmic\n in both","Same rhythm\n in both","Gain of rhythm\n in GR-LKO","Loss of rhythm\n in GR-LKO","Altered rhythm\n in GR-LKO","Filtered due\n to low counts","Unassigned\n category"))) %>%
  mutate(barcol=factor(ifelse(as.numeric(cat_star) %in% c(6,7),1,0))) %>%
  ggplot(aes(x=cat_star,fill=barcol,group=barcol)) + 
  geom_bar(colour="black",alpha=0.5) + geom_text(stat='count', aes(label=..count..), vjust=-1) + geom_vline(xintercept = 5.5,linetype="dashed") +
  theme_bw() +  theme(axis.title = element_blank(),axis.text.x=element_text(angle=30, hjust=0.5, vjust=.7)) + scale_fill_manual(values=c("black","grey")) + guides(fill=FALSE)


res_merge %>% dplyr::select(geneID,cat_star) %>%
  mutate(cat_star=ifelse(is.na(cat_star),"unassigned",cat_star)) %>%
  merge(FC_RNA_norm %>% as_tibble(rownames = "geneID") %>% dplyr::select(geneID),by="geneID",all=T) %>%
  mutate(cat_star=ifelse(is.na(cat_star),"filtered",cat_star)) %>%
  mutate(cat_star=factor(cat_star,levels=c("arrhy","same","gain*","loss*","change*","filtered","unassigned"),
                         labels=c("Arrhythmic\n in both","Same rhythm\n in both","Gain of rhythm\n in GR-LKO","Loss of rhythm\n in GR-LKO","Altered rhythm\n in GR-LKO","Filtered due\n to low counts","Unassigned\n category"))) %>%
  mutate(barcol=factor(ifelse(as.numeric(cat_star) %in% c(6,7),1,0))) %>%
  mutate(barxmin=ifelse(as.numeric(cat_star) %in% c(6,7),as.numeric(cat_star)+.5-.35,as.numeric(cat_star)-.35)) %>%
  mutate(barxmax=ifelse(as.numeric(cat_star) %in% c(6,7),as.numeric(cat_star)+.5+.35,as.numeric(cat_star)+.35)) %>%
  group_by(cat_star) %>%
  dplyr::mutate(freq=length(cat_star)) %>%
  ungroup %>%
  dplyr::select(cat_star,barcol,barxmin,barxmax,freq) %>% 
  unique %>%
  ggplot(aes(y=freq,fill=barcol,group=barcol)) + 
  geom_rect(aes(xmin=barxmin,xmax=barxmax,ymax=freq),ymin=0,colour="black",alpha=0.5,stat="identity") + geom_text(aes(x=barxmin+.35,label=freq), vjust=-1) + geom_vline(xintercept = 5.75,linetype="dashed") +
  theme_bw() +  theme(axis.title = element_blank(),axis.text.x=element_text(angle=30, hjust=0.5, vjust=.7)) +
  scale_fill_manual(values=c("black","grey")) + 
  scale_y_continuous(breaks=c(0,5000,10000),limits = c(0,11500)) +
  scale_x_continuous(breaks=c(1:5,6.5,7.5),labels=c("Arrhythmic\n in both","Same rhythm\n in both","Gain of rhythm\n in GR-LKO","Loss of rhythm\n in GR-LKO","Altered rhythm\n in GR-LKO","Filtered due\n to low counts","Unassigned\n category")) + guides(fill=FALSE) -> cat_bar

ggsave(cat_bar,filename = "GRplots/category_barchart.png",height = 3.5,width = 5.25)



mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "useast.ensembl.org")
mousemart_out<-getBM(attributes = c("entrezgene_id","mgi_symbol"), filters = "entrezgene_id", values = rownames(FC_RNA_norm) ,
                     mart = mouse, uniqueRows=T)


merge(res_merge,mousemart_out,by.x="geneID",by.y="entrezgene_id",all = TRUE) -> ll


ll %>% filter(cat_star=="change*") %>%dplyr::select(mgi_symbol) %>% view
ll %>% filter(cat_star=="loss*") %>%dplyr::select(mgi_symbol) %>% view



###################

res_merge %>% write.csv("~/GRanalysis_master/Fig 4/comparerhythms_out.csv")

read.csv("~/GRanalysis_master/Fig 4/conv_david.txt",row.names = NULL,sep="\t") -> conv_david

conv_david %>% head

res_merge %>% merge(conv_david,by.x="geneID",by.y="From",all.x = T) %>%write.csv("~/GRanalysis_master/Fig 4/comparerhythms_mergeDavid.csv")


#######


liver_phases %>%
  mutate(phase=round(phase/2,0)*2) %>%
  mutate(phase=ifelse(phase>23,phase-24,phase)) %>%
  group_by(cat_star,group,phase) %>%
  mutate(freq=n()) %>%
  ungroup %>%
  group_by(cat_star,group)%>%
  mutate(prop=freq/n()) %>%
  ungroup %>%
  dplyr::select(-id) %>%
  unique %>%
  filter(cat_star!="arrhy") %>%
  mutate(cat_group=paste0(cat_star,group)) %>%
  filter(cat_group%in%c("sameWT_phase",
                        "lossWT_phase",
                        "changeWT_phase",
                        "changeHET_phase",
                        "gainHET_phase")) -> phase_plot_data
  
phase_plot_data %>%
  ggplot(aes(x=phase,y=prop)) + 
  geom_bar(stat="identity",color="black",fill="grey",size=0.25,width=2) +
  theme_bw()+
  theme(panel.border=element_blank(),axis.title.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
  coord_polar(start=-pi/12) +
  scale_x_continuous(breaks= seq(0,24,by=2)) + facet_wrap(vars(cat_group))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
(phase_plot_data %>%
    filter(cat_group=="sameWT_phase") %>%
    ggplot(aes(x=phase,y=prop)) + 
    geom_bar(stat="identity",color="black",fill=cbPalette[2],size=0.25,width=2) +
    theme_bw()+
    theme(panel.border=element_blank(),axis.title.x=element_blank(),
          axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
    coord_polar(start=-pi/12) +
    scale_x_continuous(breaks= seq(0,24,by=2)) + ggtitle("Same rhythm in both") + theme(plot.title = element_text(hjust = 0.5),plot.caption = element_text(size=12))+ labs(subtitle = "p < 0.01"))+
  (phase_plot_data %>%
     filter(cat_group=="lossWT_phase") %>%
     ggplot(aes(x=phase,y=prop)) + 
     geom_bar(stat="identity",color="black",fill=cbPalette[3],size=0.25,width=2) +
     theme_bw()+
     theme(panel.border=element_blank(),axis.title.x=element_blank(),
           axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
     coord_polar(start=-pi/12) +
     scale_x_continuous(breaks= seq(0,24,by=2))+ ggtitle("Loss of rhythm")+ theme(plot.title = element_text(hjust = 0.5),plot.caption = element_text(size=12)) + labs(subtitle = "p = 0.34"))+
  (phase_plot_data %>%
     filter(cat_group=="gainHET_phase") %>%
     ggplot(aes(x=phase,y=prop)) + 
     geom_bar(stat="identity",color="black",fill=cbPalette[4],size=0.25,width=2) +
     theme_bw()+
     theme(panel.border=element_blank(),axis.title.x=element_blank(),
           axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
     coord_polar(start=-pi/12) +
     scale_x_continuous(breaks= seq(0,24,by=2))+ ggtitle("Gain of rhythm")+ theme(plot.title = element_text(hjust = 0.5),plot.caption = element_text(size=12))+ labs(subtitle = "p < 0.01")) -> p_phase

ggsave("GRplots/p_phase.png",p_phase,height =5,width = 6)

###

phase_plot_data %>%
  filter(cat_group%in%c("lossWT_phase")) -> loss_vec
rep(loss_vec$phase,loss_vec$freq) %>% circular(units="hours") %>% rayleigh.test()

phase_plot_data %>%
  filter(cat_group%in%c("sameWT_phase")) -> same_vec
rep(same_vec$phase,same_vec$freq) %>% circular(units="hours") %>% rayleigh.test()

phase_plot_data %>%
  filter(cat_group%in%c("gainHET_phase")) -> gain_vec
rep(gain_vec$phase,gain_vec$freq) %>% circular(units="hours") %>% rayleigh.test()


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

ggsave("GRplots/p_phase_p.png",p_phase_p,height =3.5,width = 6.5)
