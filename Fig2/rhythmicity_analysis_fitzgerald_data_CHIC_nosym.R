library(tidyverse)
library(MetaCycle)
library(patchwork)
library(nnet)

read.csv("Data/pascho.csv",sep = ",",header = F) ->RNA_1

RNA_1 %>% head
tibble(Gene=RNA_1$V3,
       RNA1=as.numeric(RNA_1$V9),
       RNA2=as.numeric(RNA_1$V10),
       RNA3=as.numeric(RNA_1$V11),
       RNA4=as.numeric(RNA_1$V12),
       RNA5=as.numeric(RNA_1$V13),
       RNA6=as.numeric(RNA_1$V14),
       RNA7=as.numeric(RNA_1$V15),
       RNA8=as.numeric(RNA_1$V16),
       RNA9=as.numeric(RNA_1$V17),
       RNA10=as.numeric(RNA_1$V18),
       RNA11=as.numeric(RNA_1$V19),
       RNA12=as.numeric(RNA_1$V20),
       RNA13=as.numeric(RNA_1$V21),
       RNA14=as.numeric(RNA_1$V22),
       RNA15=as.numeric(RNA_1$V23),
       RNA16=as.numeric(RNA_1$V24),
       RNA17=as.numeric(RNA_1$V25),
       RNA18=as.numeric(RNA_1$V26),
       RNA19=as.numeric(RNA_1$V27),
       RNA20=as.numeric(RNA_1$V28),
       RNA21=as.numeric(RNA_1$V29),
       RNA22=as.numeric(RNA_1$V30),
       RNA23=as.numeric(RNA_1$V31),
       RNA24=as.numeric(RNA_1$V32)) -> RNA_rejig
##
RNA_rejig <- RNA_rejig[-1,]
RNA_rejig[,-1] %>% apply(2,sum) -> count_sums
RNA_n <- RNA_rejig
RNA_n[,-1] <- (RNA_rejig[,-1]/matrix(rep(count_sums,23733),nrow=23733,byrow=TRUE))  #normalised by dividing by total number of reads per sample
RNA_n[,-1] %>% apply(1,var) -> count_vars

# JTK cycle analysis ------------------------------------------------------


write.table(RNA_n, file="Data/RNA_n_forJTK.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

Meta_Output <- meta2d(infile="Data/RNA_n_forJTK.txt", filestyle="txt",
                      outputFile = FALSE, timepoints=rep(seq(0, 20, by=4),each=4),
                      cycMethod=c("JTK"), outIntegration="both",minper = 24,maxper = 24)

cbind(RNA_n,(Meta_Output$meta %>% dplyr::select(JTK_pvalue,JTK_adjphase))) -> RNA_n_jtk

RNA_n_jtk %>%
  filter(JTK_pvalue<0.05) %>%
  ggplot(aes(x=JTK_adjphase))+
  geom_bar(color="black",fill="red",size=0.25,width=2) +
  theme_bw()+
  theme(panel.border=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
  coord_polar(start=-pi/12) +ggtitle("GR bound and Rhythmic")+
  scale_x_continuous(breaks=c(0,3,6,9,12,15,18,21,24))

#
read.csv("Data/conv_pascho.txt",sep = "\t") ->conv_1

RNA_n_jtk %>% merge(conv_1 %>% dplyr::select(From,To),by.x = "Gene",by.y="From",all.x = T)  -> ll
RNA_rejig  %>% merge(conv_1 %>% dplyr::select(From,To),by.x = "Gene",by.y="From",all.x = T)  -> ll_unorm

ll  -> RNA_n_jtk_2

#####################

# library(biomaRt)
# 
# ensdb <- useEnsembl(biomart="genes", assembly=38, host="grch38.ensembl.org", path="/biomart")

RNA_n_jtk_2

library(Rsubread)
library(GenomicRanges)
# TSS <- promoterRegions("mm10", upstream=400, downstream=2000)
# 
# TSS[TSS$GeneID%in%RNA_n_jtk_2$To,] -> TSS_ourproms 
#   
# TSS_ourproms_Granges <- GRanges(seqnames= Rle(TSS_ourproms$Chr), ranges = IRanges(TSS_ourproms$Start, TSS_ourproms$End))

load("~/GRanalysis_master/GR_CircadianLiverTranscriptome/Data/promoter_data_final.RData")
data_master <- data_master_rev
data_master %>% dplyr::select(PROMID,GENEID) %>% separate(PROMID,into = c("chr","start","end"),remove = FALSE) %>% unique %>%
  mutate(start=as.numeric(start)) %>% mutate(end=as.numeric(end))-> prom_regions

prom_regions %>% filter(!is.na(start)) -> prom_regions

prom_regions <- makeGRangesFromDataFrame(prom_regions,
                                         keep.extra.columns = TRUE,
                                         seqnames.field = "chr",
                                         start.field = "start",
                                         end.field = "end")
TSS_ourproms_Granges <- prom_regions

chic<-read.csv("~/GRanalysis_master/GR_CircadianLiverTranscriptome/Data/All_step2_washU_text_3.txt", sep = "\t", header = F)
chic[,7]<- 1:length(chic[,1])
#Symmetric set is the complete one
chic_baits<-chic[,c(1:3,7)]
chic_otherends<-chic[,4:7]

#Bedfiles, GRanges Objects
chic_otherends_bed<-GRanges(seqnames= Rle(chic_otherends[,1]), ranges = IRanges(chic_otherends[,2], chic_otherends[,3]))
values(chic_otherends_bed) <- chic_otherends[,4]
chic_bait_bed<-GRanges(seqnames= Rle(chic_baits[,1]), ranges = IRanges(chic_baits[,2], chic_baits[,3]))
values(chic_bait_bed) <- chic_baits[,4]

chic_bait_bed<-easyLift::easyLiftOver(chic_bait_bed,map = "mm9_mm10")
chic_otherends_bed<-easyLift::easyLiftOver(chic_otherends_bed,map = "mm9_mm10")

prom_bait_overlaps<-findOverlaps(query = TSS_ourproms_Granges, subject = chic_bait_bed)

load("Data/enhancer_withRNA.RData")
enh_withRNA %>% dplyr::select(ENHANCERID) %>% unique %>% separate(ENHANCERID,into=c("seq","start","end")) %>%
  mutate(start=as.numeric(start)) %>% mutate(end=as.numeric(end))-> ourenhs
ourenhs_Granges <- GRanges(seqnames= Rle(ourenhs$seq), ranges = IRanges(ourenhs$start, ourenhs$end))

enh_otherend_overlaps<-findOverlaps(query = ourenhs_Granges, subject = chic_otherends_bed)

#to(prom...) returns baits from overlap. So chic_bait_bed[to(prom..)] returns bait regions in order. $X gives link to otherend

#chic_filt <- chic[(chic[,7] %in% chic_bait_bed[to(prom_bait_overlaps)]$X ) & (chic[,7] %in% chic_otherends_bed[to(enh_otherend_overlaps)]$X ),]

## rewrite chic_filt with direct indexing
# Get bait indices from overlaps with promoters
bait_indices <- subjectHits(prom_bait_overlaps)

# Get otherEnd indices from overlaps with enhancers
otherend_indices <- subjectHits(enh_otherend_overlaps)

# Get CHi-C row IDs from bait and otherEnd GRanges
bait_ids <- chic_bait_bed[bait_indices] %>% mcols() %>% unlist() #every interaction index that is linked to a bait that overlaps with a prom
otherend_ids <- chic_otherends_bed[otherend_indices] %>% mcols() %>% unlist() #every interaction index that is linked to a otherend that overlaps with a enh

# Filter CHi-C rows where the ID is found in both bait and otherEnd sets
chic_filt <- chic[chic[,7] %in% intersect(bait_ids, otherend_ids), ] #interactions that link to both prom and enh inour dataset
#########
enh_prom_links<-data.frame()
for(i in chic_filt[,7]){ #this is interaction index
  which(chic_bait_bed$X==i) -> cbbi
  which(chic_otherends_bed$X==i) -> coebi
  TSS_ourproms_Granges[from(prom_bait_overlaps[to(prom_bait_overlaps)==cbbi]),] %>% as.data.frame() -> temp_proms
  names(temp_proms) <- temp_proms %>% names %>% paste0("prom_",.)
  ourenhs_Granges[from(enh_otherend_overlaps[to(enh_otherend_overlaps)==coebi]),] %>% as.data.frame() -> temp_enh
  names(temp_enh) <- temp_enh %>% names %>% paste0("enh_",.)
  if((dim(temp_enh)[1]>=1)&(dim(temp_proms)[1]>=1)){
    if(dim(temp_proms)[1]>=2){
      for(j in 1:(dim(temp_proms)[1])){
        enh_prom_links<-rbind(enh_prom_links,cbind(temp_proms[j,],temp_enh))
      }
    }else{
  enh_prom_links<-rbind(enh_prom_links,cbind(temp_proms,temp_enh))
    }
  }
}

save(enh_prom_links,file = "Data/enh_prom_links_CHIC_3b.RData")
#####################

load("Data/enhancer_withRNA.RData")

# enh_withRNA %>% head
# merge(RNA_n_jtk_2,enh_withRNA,by.x="To",by.y="GENEID") -> data_master_all


# enh_prom_links %>% mutate(PROMID=paste(prom_seqnames,prom_start,prom_end,sep=".")) %>% 
# merge(TSS_ourproms%>% mutate(PROMID=paste(Chr,Start,End,sep=".")),by="PROMID") %>% 
#   merge(.,RNA_n_jtk_2,by.y="To",by.x="GeneID") %>%
#   mutate(ENHANCERID=paste(enh_seqnames,enh_start,enh_end,sep=".")) %>%
#   merge(.,enh_withRNA,by.x="ENHANCERID",by.y="ENHANCERID") -> data_master_all

enh_prom_links %>% mutate(PROMID=paste(prom_seqnames,prom_start,prom_end,sep=".")) %>% 
  merge(prom_regions%>% as.data.frame,by="PROMID") %>% 
  merge(.,RNA_n_jtk_2,by.y="To",by.x="GENEID") %>%
  mutate(ENHANCERID=paste(enh_seqnames,enh_start,enh_end,sep=".")) %>%
  merge(.,enh_withRNA,by.x="ENHANCERID",by.y="ENHANCERID") %>%
  mutate(GeneID=GENEID.y)-> data_master_all

# merge(ll_unorm,enh_withRNA,by.x="To",by.y="GENEID") -> data_master_all_unorm
# save(data_master_all_unorm,file = "Data/paschodata_enhancer.RData")

(data_master_all %>% 
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
    filter(Both_us.g==TRUE) %>%
    ggplot(aes(x=JTK_adjphase))+
    geom_bar(color="black",fill="red",size=0.25,width=2) +
    theme_bw()+
    theme(panel.border=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
    coord_polar(start=-pi/12) +ggtitle("GR bound and Rhythmic")+ labs(subtitle = "p < 0.01") + 
    scale_x_continuous(breaks=c(0,3,6,9,12,15,18,21))) +
  
  (data_master_all %>% 
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
     filter(Both_us.g==FALSE) %>%
     ggplot(aes(x=JTK_adjphase))+
     geom_bar(color="black",fill="red",size=0.25,width=2) +
     theme_bw()+
     theme(panel.border=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
     coord_polar(start=-pi/12)+ggtitle("Not GR bound and Rhythmic") + labs(subtitle = "p < 0.01")) +  
  scale_x_continuous(breaks=c(0,3,6,9,12,15,18,21))-> pplots

#
data_master_all %>% 
  mutate(JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
  dplyr::select(GeneID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GeneID,Rhythmic,Both_us,JTK_adjphase) %>% 
  as_tibble %>% group_by(GeneID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(GeneID,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
  unique -> phase_plot_data

phase_plot_data %>%
  filter(Rhythmic.g==TRUE) %>%
  filter(Both_us.g==FALSE) -> rhynoGR_vec
rhynoGR_vec$JTK_adjphase %>% circular(units="hours") %>% rayleigh.test()

phase_plot_data %>%
  filter(Rhythmic.g==TRUE) %>%
  filter(Both_us.g==TRUE) -> rhyyesGR_vec
rhyyesGR_vec$JTK_adjphase %>% circular(units="hours") %>% rayleigh.test()

watson.two.test(rhyyesGR_vec$JTK_adjphase %>% circular(units="hours"),
                rhynoGR_vec$JTK_adjphase %>% circular(units="hours"))

ggplot()+
  geom_segment(aes(y = 1,yend= 2, x=2.5))+
  geom_segment(aes(y = 1,yend= 2, x=4.5))+
  geom_segment(aes(y = 1, x=2.5,xend=4.5)) + 
  geom_text(aes(label="p<0.01",x=3.5,y=0.7)) +
  theme_void() +
  xlim(c(1.7,5.3)) + ylim(c(-1.5,2)) -> p_pval

pplots/p_pval + plot_layout(heights=c(2,.6)) -> p_phase_p

ggsave("Fig2/Plots/p_phase_p_chic.png",p_phase_p,height =4,width = 6.5)

#

data_master_all %>% 
  mutate(JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
  dplyr::select(GeneID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GeneID,Rhythmic,Both_us) %>% 
  as_tibble %>% group_by(GeneID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(GeneID,Rhythmic.g,Both_us.g) %>%
  unique %>%
  dplyr::select(-GeneID) %>%
  table 

# Venn Diagram ------------------------------------------------------------

data_master_all %>% 
  mutate(JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
  dplyr::select(PROMID.x,GeneID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GeneID,Rhythmic,Both_us) %>% 
  as_tibble %>% group_by(GeneID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(GeneID,Rhythmic.g,Both_us.g) %>%
  unique %>%
  dplyr::select(-GeneID) -> venn_data

library(VennDiagram) 
venn.diagram(list(`GR Bound` = which(venn_data$Both_us.g),Rhythmic = which(venn_data$Rhythmic.g)), 
             fill = c("white", "white"), 
             disable.logging = TRUE,
             category.names = c("" , "" ),
             alpha = c(0.5, 0.5),
             lwd =1,
             cat.cex=0.6,
             height=1500,width=1500, "Fig2/Plots/venn_diagram_chic.png")

### E ###
data_master_all %>% 
  mutate(JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
  mutate(RNA_m=data_master_all %>% dplyr::select(RNA1:RNA24) %>% apply(1,mean)) %>%
  dplyr::select(GeneID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,RNA_m) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GeneID,Rhythmic,Both_us,RNA_m) %>% 
  as_tibble %>% group_by(GeneID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup -> data_model

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
  mutate(GR_any=c("No GR binding", "GR binding")[GR_any+1]) %>%
  ggplot(aes(x=m.RNA,y=GRrhy)) +
  geom_jitter(colour="grey",alpha=.5,height=.2,shape=16) +
  geom_boxplot(fill=NA) +
  #scale_x_log10(name="RPKM") +
  # scale_y_continuous(name="Density")+
  #scale_colour_manual(values = c("red","black"))+
  ylab("") +
  theme_bw() +
  theme(legend.position = c(0.75, .92),legend.background = element_blank()) +
  guides(colour="none")  + scale_x_log10("Normalised Counts") -> RNA_plot

ggsave(RNA_plot,filename = "Fig2/Plots/rna_boxplots_fitz_chic.png",height = 3,width = 4)


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
  mutate(GR_any=c("No GR binding", "GR binding")[GR_any+1])-> ll


table(ll$Rhy,ll$GR_any)
res.aov2 <- aov(m.RNA ~ Rhy + GR_any, data = ll)
summary(res.aov2)


# homer bed files ---------------------------------------------------------

data_master_all %>% 
  mutate(GENEID=To,JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
  dplyr::select(PROMID,GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  as_tibble %>% group_by(GENEID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(ENHANCERID,Rhythmic.g,Both_us.g) %>%
  unique %>%
  filter(Rhythmic.g==TRUE)  %>% 
  filter(Both_us.g==TRUE) %>%
  dplyr::select(ENHANCERID) %>% separate(ENHANCERID,into = c("chrom","chromStart","chromEnd")) %>%
  unique -> temp_dat

temp_dat %>%
  mutate(id=1:(dim(temp_dat)[1])) %>%
  mutate(` `="") %>%
  mutate(strand=".") %>%
  write.table(file = "Data/GRboundRhythmicenhancers_forHomer.bed",sep = "\t",quote = FALSE,row.names = F)

# modelling ---------------------------------------------------------------

data_master_all %>% 
  mutate(GENEID=GeneID,JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
  mutate(RNA_m=data_master_all %>% dplyr::select(RNA1:RNA24) %>% apply(1,mean)) %>%
  dplyr::select(GENEID,ENHANCERID,RNA_m,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID,Rhythmic,Both_us,RNA_m,JTK_adjphase) %>% 
  as_tibble %>% group_by(GENEID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us))  %>%
  dplyr::select(GENEID,Rhythmic.g,Both_us.g,RNA_m,JTK_adjphase) %>%
  mutate(RNA_m_log10=log10(RNA_m)) %>% unique -> data_model

# logistic ----------------------------------------------------------------

glm(data = data_model %>% filter(RNA_m!=0),formula = "Rhythmic.g~RNA_m_log10 + Both_us.g",family = binomial(link="logit")) -> mod3
summary(mod3)

exp(cbind(coef(mod3)[2:3],confint.default(mod3,2:3))) %>%
  as_tibble(rownames="row")

# expression stratified not adjusted for expression  -------------------------------------------------------

data_model%>% filter(RNA_m!=0) %>% mutate(RNA_bin=cut(RNA_m_log10,quantile(data_model$RNA_m_log10,probs = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)),include.lowest=TRUE)) %>%
  mutate(RNA_bin=as.numeric(RNA_bin))-> data_model_2

quartOR <- c()
for(i in 1:10){
  glm(data = data_model_2 %>% filter(RNA_bin==i),formula = "Rhythmic.g~Both_us.g",family = binomial(link="logit")) -> mod3
  exp(cbind(coef(mod3)[2],confint.default(mod3,2))) %>%
    as_tibble(rownames="row") %>% .[1,] -> temp
  temp$row <- paste("quantile",i)
  quartOR<- rbind(quartOR,temp)
}
data_model_2 %>% ungroup %>% group_by(RNA_bin) %>% summarise(total=length(Rhythmic.g),numberRhy=sum(Rhythmic.g),numberGR=sum(Both_us.g),meanRNA=mean(RNA_m),medRNA=median((RNA_m))) %>%
  dplyr::select(-RNA_bin) %>% cbind(quartOR,.) %>%
  as_tibble  -> quantiletable

quantiletable %>%
  mutate(OR=ifelse(V1>100,NA,round(V1,digits=2))) %>%
  mutate(`95% CI`=paste(round(`2.5 %`,2),"-",round(`97.5 %`,2))) %>%
  dplyr::select(row,OR,`95% CI`,total,numberRhy,numberGR) %>% view

pd_width <- 0.6
quantiletable %>% 
  mutate(row=factor(row,ordered=T,levels=unique(quantiletable$row))) %>% 
  #filter(row!="quantile 1") %>%
  ggplot(aes(x=row,y=V1)) +
  geom_hline(aes(yintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbar(aes(ymax = `97.5 %`, ymin = `2.5 %`), size = .5, width = .4,
                position = position_dodge(width = pd_width)) +
  geom_point(position = position_dodge(width = pd_width)) + theme_bw() + ylab("OR") + xlab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))-> p_boxplot

ggsave(filename = "Fig2/Plots/supp_quants_chic.png",width=5,height=4)

### supp 2B 
prob_len <- 50
prob_seq <- seq(0,1,length=prob_len+1)

data_model%>% filter(RNA_m!=0, Both_us.g==1) -> data_model_GR

## quantiles separately for GR and no GR
data_model_GR %>% mutate(RNA_bin=cut(RNA_m_log10,quantile(data_model_GR$RNA_m_log10,probs = prob_seq),include.lowest=TRUE)) %>%
  mutate(RNA_bin=as.numeric(RNA_bin)) %>% mutate(GR="GR binding") -> data_model_GR
data_model%>% filter(RNA_m!=0, Both_us.g==0) -> data_model_noGR
data_model_noGR %>% mutate(RNA_bin=cut(RNA_m_log10,quantile(data_model_noGR$RNA_m_log10,probs = prob_seq),include.lowest=TRUE)) %>%
  mutate(RNA_bin=as.numeric(RNA_bin))%>% mutate(GR="No GR binding")-> data_model_noGR

rbind(data_model_GR,data_model_noGR) -> data_model_together

data_model_together %>% ungroup %>% group_by(RNA_bin,GR) %>% summarise(total=length(Rhythmic.g),numberRhy=sum(Rhythmic.g),numberGR=sum(Both_us.g),meanRNA=mean(RNA_m),medRNA=median((RNA_m))) %>%
  mutate(percentRhy=numberRhy/total) %>%
  dplyr::select(-RNA_bin) %>%
  as_tibble  -> quantiletable

##gam
quantiletable %>%
  ggplot(aes(x=medRNA,y=percentRhy,colour=GR,fill=GR)) + 
  geom_point(alpha=0.3) +
  scale_x_log10() +
  theme_bw() +
  geom_smooth(alpha=.3,method = "gam") +
  ylab("Proportion Rhythmic")+
  xlab("Median RNA for quantile") +
  theme(legend.title = element_blank())+
  scale_colour_manual(values = c("red","black"))+
  scale_fill_manual(values = c("red","black")) -> p

ggsave(filename = "Fig2/Plots/suppB_chic.png",p,width=6,height=4)
