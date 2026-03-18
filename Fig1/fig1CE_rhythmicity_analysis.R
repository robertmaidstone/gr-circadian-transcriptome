
library(MetaCycle)
library(patchwork)
library(nnet)
library(circular)

read.csv("Data/pascho.csv",sep = ",",header = F) ->RNA_1


library(tidyverse)

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

#
read.csv("Data/conv_pascho.txt",sep = "\t") ->conv_1

RNA_n_jtk %>% merge(conv_1 %>% dplyr::select(From,To),by.x = "Gene",by.y="From",all.x = T)  -> ll

RNA_rejig  %>% merge(conv_1 %>% dplyr::select(From,To),by.x = "Gene",by.y="From",all.x = T)  -> ll_unorm

ll  -> RNA_n_jtk_2

load("Data/RNAChIPtogether_GRdist.RData")

data_master %>% head
merge(RNA_n_jtk_2,data_master,by.x="To",by.y="GENEID") -> data_master_all
merge(ll_unorm,data_master,by.x="To",by.y="GENEID") -> data_master_all_unorm
save(data_master_all_unorm,file = "Data/paschodata.RData")

(data_master_all %>% 
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
    filter(Both_us.g==TRUE) %>%
    ggplot(aes(x=JTK_adjphase))+
    geom_bar(color="black",fill="red",size=0.25,width=2) +
    theme_bw()+
    theme(panel.border=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
    coord_polar(start=-pi/12) +ggtitle("GR bound and Rhythmic")+
    scale_x_continuous(breaks=c(0,3,6,9,12,15,18,21))+ labs(subtitle = "p < 0.01")) + #pvalue from rayleigh test below
  
  (data_master_all %>% 
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
     filter(Both_us.g==FALSE) %>%
     ggplot(aes(x=JTK_adjphase))+
     geom_bar(color="black",fill="red",size=0.25,width=2) +
     theme_bw()+
     theme(panel.border=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
     coord_polar(start=-pi/12)+ggtitle("Not GR bound and Rhythmic")+ labs(subtitle = "p < 0.01")) + #pvalue from rayleigh test below
  scale_x_continuous(breaks=c(0,3,6,9,12,15,18,21))-> pplots

#
data_master_all %>% 
  mutate(JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
  dplyr::select(PROMID,To,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(To,Rhythmic,Both_us,JTK_adjphase) %>% 
  as_tibble %>% group_by(To) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(To,Rhythmic.g,Both_us.g,JTK_adjphase) %>%
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
  geom_text(aes(label="p<0.01",x=3.5,y=0.7)) + #from watson test
  theme_void() +
  xlim(c(1.7,5.3)) + ylim(c(-1.5,2)) -> p_pval

pplots/p_pval + plot_layout(heights=c(2,.6)) -> p_phase_p

ggsave("Fig1/Plots/p_phase_p.png",p_phase_p,height =4,width = 6.5)

# Venn diagram

data_master_all %>% 
  mutate(GENEID=To,JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
  dplyr::select(PROMID,GENEID,ENHANCERID,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID,Rhythmic,Both_us) %>% 
  as_tibble %>% group_by(GENEID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us)) %>% 
  ungroup %>%
  dplyr::select(GENEID,Rhythmic.g,Both_us.g) %>%
  unique %>%
  dplyr::select(-GENEID) -> venn_data

library(VennDiagram) 
venn.diagram(list(`GR Bound` = which(venn_data$Both_us.g),Rhythmic = which(venn_data$Rhythmic.g)), 
             fill = c("white", "white"), 
             category.names = c("" , "" ),
             alpha = c(0.5, 0.5),
             lwd =1,
             cat.cex=0.6,
             disable.logging = TRUE,
             height=1500,width=1500, "Fig1/Plots/venn_diagram.png")

# modelling to control for expression

data_master_all %>% 
  mutate(GENEID=To,JTK_pvalue=JTK_pvalue.x,JTK_adjphase=JTK_adjphase.x) %>%
  mutate(RNA_m=data_master_all %>% dplyr::select(RNA1:RNA24) %>% apply(1,mean)) %>%
  dplyr::select(PROMID,GENEID,ENHANCERID,RNA_m,WT_us=(GRlou_WT_us),Dex_us=GRlou_dex_us,JTK_pvalue,JTK_adjphase) %>%
  unique %>%
  mutate(Both_us=(Dex_us==0),Rhythmic=JTK_pvalue<0.05) %>% 
  dplyr::select(GENEID,Rhythmic,Both_us,RNA_m,JTK_adjphase) %>% 
  as_tibble %>% group_by(GENEID) %>% dplyr::mutate(Rhythmic.g=any(Rhythmic),Both_us.g=any(Both_us))  %>%
  dplyr::select(GENEID,Rhythmic.g,Both_us.g,RNA_m,JTK_adjphase) %>%
  mutate(RNA_m_log10=log10(RNA_m)) %>% unique -> data_model


# logistic (figure 1E) ----------------------------------------------------------------

glm(data = data_model %>% filter(RNA_m!=0),formula = "Rhythmic.g~RNA_m_log10 + Both_us.g",family = binomial(link="logit")) -> mod3
summary(mod3)

exp(cbind(coef(mod3)[2:3],confint.default(mod3,2:3))) %>%
  as_tibble(rownames="row")

# expression stratified (Supplemental 1B) -------------------------------------------------------

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
  dplyr::select(row,OR,`95% CI`,total,numberRhy,numberGR) #%>% view

pd_width <- 0.6
quantiletable %>% 
  mutate(row=factor(row,ordered=T,levels=unique(quantiletable$row))) %>% 
  filter(row!="quantile 1") %>% ggplot(aes(x=row,y=V1)) +
  geom_hline(aes(yintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbar(aes(ymax = `97.5 %`, ymin = `2.5 %`), size = .5, width = .4,
                position = position_dodge(width = pd_width)) +
  geom_point(position = position_dodge(width = pd_width)) + theme_bw() + ylab("OR") + xlab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))-> p_boxplot

ggsave(filename = "Fig1/Plots/supp1B.png",width=5,height=4)

# quantiletable %>%
#   mutate(row=factor(row,ordered=T,levels=unique(quantiletable$row))) %>% 
#   ggplot(aes(x=medRNA,y=numberRhy)) + geom_point() + scale_x_log10()


# proportion rhythmic (supplemental 1A) -----------------------------------

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
  geom_smooth(alpha=.3,method = "gam",fullrange=TRUE) +
  ylab("Proportion Rhythmic")+
  xlab("Median RNA for quantile") +
  theme(legend.title = element_blank())+
  scale_colour_manual(values = c("red","black"))+
  scale_fill_manual(values = c("red","black")) -> p

ggsave(filename = "Fig1/plots/supp1A.png",p,width=6,height=4)
