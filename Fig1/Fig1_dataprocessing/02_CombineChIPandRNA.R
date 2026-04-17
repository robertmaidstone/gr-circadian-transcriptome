library(tidyverse)
library(edgeR)
library(pracma)
setwd("Data")


RSS_fit<-c()
RSS_test<-c()
RSS_fit_int<-c()
RSS_test_int<-c()
coeffs_sep<-c()

for(timepoint in c(0,4,8,12,16,20)){
  load(paste0("CT",timepoint,"data_nc.RData"))
  
  READS_filt<-eval(parse(text=paste0("READS_CT",timepoint,"_filt")) )
  
  ############## filter out   blacklisted #############
  
  #BFs<-read.csv("ALL_BF_no_rep.out",sep="")
  
  # intersect using bedtools
  # bedtools intersect -bed -wa -a Sites.bed -b mm10.blacklist.bed > Out_CT0.bed 
  BlackRegions<- read.csv("BL_proms.bed",sep="\t",header = F)
  
  !(paste(READS_filt$Chr,READS_filt$Start,READS_filt$End) %in% paste(BlackRegions$V1,BlackRegions$V2,BlackRegions$V3)) -> bl_logic
  READS_filt %>%
    filter(bl_logic) -> READS_filt

  ####################################################################
  ################# Convert histone marks to CPM #####################
  c(4744198,6970058,4492392,5671641,4828132,11796353) -> h3k27ac_filt_tags
  c(24299934,15842142,17336065,19634619,23683261,18270194) -> h3k4me3_filt_tags
  
  READS_filt$H3K27ac<-(READS_filt$H3K27ac*1000000)/h3k27ac_filt_tags[which(timepoint== c(0,4,8,12,16,20))]
  READS_filt$H3K4me3<-(READS_filt$H3K4me3*1000000)/h3k4me3_filt_tags[which(timepoint== c(0,4,8,12,16,20))]
  ####################################################################
  
 # cbind(READS_filt %>% dplyr::select(ID,Chr,Start,End,maxBF) , log(READS_filt %>% dplyr::select(-ID,-Chr,-Start,-End,-maxBF)) ) %>%
 #   as_tibble() ->READS_filt
  
  assign(paste0("READS_CT",timepoint,"_filt"),READS_filt)
  
}



READS_CT0_filt
READS_CT0_filt %>% 
  dplyr::rename(seqnames=Chr,starts=Start,ends=End) %>%
  dplyr::select(c(2,3,4)) %>%
  write.table(file="ChIP_regions.bed", quote=F, sep="\t", row.names=F, col.names=F)

##########################################################

load("~/Gene_Regulatory_Networks_Data/RNAseqFC.RData")
FC_RNA_temp<-FC_RNA
FC_RNA$counts %>% apply(2,sum) -> count_sums

y <- DGEList(counts=FC_RNA$counts,genes=FC_RNA$annotation)
calcNormFactors(y)-> y
y = estimateCommonDisp(y, verbose=TRUE)
rpkm_y<- rpkm(y)
FC_RNA_norm <- rpkm_y


FC_RNA_detrend <- t(detrend(t(FC_RNA_norm),tt="linear"))

##### JTK Cycle analysis of RNA counts (after normalisation)

library(MetaCycle)
write.table(cbind(rownames(FC_RNA_detrend),FC_RNA_detrend), file="FC_RNA_edgeRnorm_forJTK.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

Meta_Output <- meta2d(infile="FC_RNA_edgeRnorm_forJTK.txt", filestyle="txt",
       outputFile = FALSE, timepoints=seq(0, 44, by=4),minper=24,maxper=24,
       cycMethod=c("JTK"), outIntegration="both")

cbind(FC_RNA_norm,(Meta_Output$meta %>% dplyr::select(JTK_pvalue,JTK_BH.Q, JTK_period, JTK_adjphase, JTK_amplitude))) -> FC_RNA_norm

############################################################

library("TxDb.Mmusculus.UCSC.mm10.knownGene")
mm10 = TxDb.Mmusculus.UCSC.mm10.knownGene  ## a shorter name, just forconvenience
columns(mm10)
GeneIDs<-as.character(FC_RNA$annotation$GeneID)

## Now at this point you have not really told me what you wanted.
## So lets start with just how you can get information about yourgenesof interest.
## Lets start by just gettting out some of the transcript information for these genes:
res2 = select(mm10, keys=GeneIDs,
              columns=c("TXNAME","TXSTRAND","TXCHROM","TXSTART","TXEND"),
              keytype="GENEID")
res2

## From your description it sounds like you wanted a GRanges objectwiththe ranges for your promoters.
## In that case you could do it like this:
proms = promoters(mm10)
## Then you can extract the txnames from before...
txnms  = unique(res2[,"TXNAME"])
## And then subset to only the promoters that have the names that you  want here
myproms = proms[mcols(proms)[,'tx_name'] %in% txnms,]
myproms <- myproms[(myproms %>% start())>0,]
tibble(seqnames=as.character(seqnames(myproms)),starts=start(myproms),ends=end(myproms)) %>% 
  write.table(file="RNA_promoters_edgeRnorm.bed", quote=F, sep="\t", row.names=F, col.names=F)
###########################################################

############# intersect using bedtools
############# bedtools intersect -bed -wa -loj -a ChIP_regions.bed -b RNA_promoters_edgeRnorm.bed > ChIP_RNAproms_intersect.bed #######################
############# bedtools intersect -bed -wa -loj -a RNA_promoters_edgeRnorm.bed -b ChIP_regions.bed > RNAproms_ChIP_intersect.bed ##################
rhythmic_prom_int<-read.csv(file="ChIP_RNAproms_intersect.bed",sep="\t",header = F)
rhythmic_prom_int %>% head
sum(rhythmic_prom_int$V4!=".")

rhythmic_prom_int %>% filter(rhythmic_prom_int$V4!=".") %>% 
  dplyr::select(c(1,2,3)) %>% unique() %>% dim

prom_rhythmic_int<-read.csv(file="RNAproms_ChIP_intersect.bed",sep="\t",header = F)
prom_rhythmic_int %>% head
sum(prom_rhythmic_int$V4!=".")

prom_rhythmic_int %>% filter(prom_rhythmic_int$V4!=".") %>%
  dplyr::select(c(1,2,3)) %>% unique() %>% dim



myproms_db <- tibble(seqnames=as.character(seqnames(myproms)),starts=start(myproms),ends=end(myproms),txname=myproms$tx_name)
temp<- merge(res2,myproms_db,by.x="TXNAME",by.y="txname")
temp %>% mutate(PROMID=paste(seqnames,starts,ends,sep=".")) -> temp
prom_rhythmic_int %>% mutate(PROMID=paste(V1,V2,V3,sep=".")) -> prom_rhythmic_int
temp2<-merge(temp,prom_rhythmic_int,by="PROMID")

temp2 %>% filter(V4!=".") %>% 
  mutate(ENHANCERID=paste(V4,V5,V6,sep=".")) %>% 
  dplyr::select(GENEID,PROMID,ENHANCERID) -> prom_enh_pairs

as_tibble(FC_RNA_norm) %>% mutate(GENEID=rownames(FC_RNA_norm)) %>%
  dplyr::rename(RNA0=RNA1_sorted.bam,
                RNA4=RNA2_sorted.bam,
                RNA8=RNA3_sorted.bam,
                RNA12=RNA4_sorted.bam,
                RNA16=RNA5_sorted.bam,
                RNA20=RNA6_sorted.bam,
                RNA24=RNA7_sorted.bam,
                RNA28=RNA8_sorted.bam,
                RNA32=RNA9_sorted.bam,
                RNA36=RNA10_sorted.bam,
                RNA40=RNA11_sorted.bam,
                RNA44=RNA12_sorted.bam) -> RNA_data

merge(prom_enh_pairs,RNA_data,by="GENEID") -> paired_RNA_data

paired_RNA_data %>% gather(Time,RNA,-GENEID,-PROMID,-ENHANCERID,-JTK_pvalue,-JTK_BH.Q, -JTK_period, -JTK_adjphase, -JTK_amplitude) %>%
  mutate(Time= factor(Time,labels=c(0,12,16,20,24,28,32,36,4,40,44,8))) -> paired_RNA_data

rbind(cbind(READS_CT0_filt,Time=0),
      cbind(READS_CT4_filt,Time=4),
      cbind(READS_CT8_filt,Time=8),
      cbind(READS_CT12_filt,Time=12),
      cbind(READS_CT16_filt,Time=16),
      cbind(READS_CT20_filt,Time=20),
      cbind(READS_CT0_filt,Time=24),
      cbind(READS_CT4_filt,Time=28),
      cbind(READS_CT8_filt,Time=32),
      cbind(READS_CT12_filt,Time=36),
      cbind(READS_CT16_filt,Time=40),
      cbind(READS_CT20_filt,Time=44)) %>% mutate(ENHANCERID=paste(Chr,Start,End,sep=".")) -> READS_temp

merge(paired_RNA_data,READS_temp,by=c("ENHANCERID","Time")) -> tempM
tempM %>% dplyr::select(-ID,-Chr,-Start,-End) -> data_master


head(data_master)
dim(data_master)
dim(data_master)[1]/12

data_master <- data_master %>% unique
head(data_master)
dim(data_master)
dim(data_master)[1]/12


save(data_master, file="RNAChIPtogether.RData")


