#source("https://bioconductor.org/biocLite.R")
#biocLite("DiffBind")


library(DiffBind)
setwd("~/rds-loudon-ray-rattray/rds-shared-loudon-ray-rattray/Takahashi_data/")

bed_to_granges <- function(file){
  df <- read.table(file,
                   header=F,
                   stringsAsFactors=F)
  
  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }
  
  if(length(df)<3){
    stop("File has less than 3 columns")
  }
  
  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  library("GenomicRanges")
  
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}

bed_to_granges("Promoter_Counting/RNA_promoters.bed") -> promoter_regions


ddmm<-dba(sampleSheet="samplesheet_TFs_nocontrol.csv")
ddmm_count<-dba.count(ddmm,peaks = promoter_regions,minOverlap=1)
save(ddmm_count,file="Promoter_Counting/promoter_counts_nov19.RData")