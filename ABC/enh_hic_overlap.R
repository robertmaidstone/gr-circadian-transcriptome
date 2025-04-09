

# Assuming your GRanges objects are:
# enhancers: GRanges object of enhancers
# hic_regions: GRanges object of Hi-C regions

library(GenomicRanges)

load("~/GRanalysis_master/Fig2/data/enhancer_master.RData")

data_master %>% dplyr::select(ENHANCERID) %>% separate(ENHANCERID,into = c("chr","start","end")) %>% unique %>%
  mutate(start=as.numeric(start)) %>% mutate(end=as.numeric(end))-> enh_regions
enh_regions_gr <- GenomicRanges::makeGRangesFromDataFrame(enh_regions)

hic_regions <- rtracklayer::import("~/GRanalysis_master/ABC/GSE104129_mm9_5000_ord.bed/5AM_5000_ord-2.bed", format = "BED")
head(hic_regions)

# Load the chain file
chain <- rtracklayer::import.chain("~/GRanalysis_master/ABC/mm9ToMm10.over.chain/mm9ToMm10.over.chain")

# Perform the LiftOver
hic_regions_mm10 <- rtracklayer::liftOver(hic_regions, chain)

# Flatten the result into a single GRanges object
hic_regions_mm10 <- unlist(hic_regions_mm10)


# Find overlaps
overlaps <- findOverlaps(enh_regions_gr, hic_regions_mm10)

# Extract the indices of overlapping Hi-C regions
hic_indices <- subjectHits(overlaps)

# Add a metadata column to the enhancers GRanges object
# e.g., Adding the names of the overlapping Hi-C regions
mcols(enh_regions_gr)$hic_region <- NA  # Initialize column
mcols(enh_regions_gr)$hic_region[queryHits(overlaps)] <- (hic_regions$name)[hic_indices]


######
##Proms

proms<-read.table(file="~/GRanalysis_master/Fig2/data/RNA_promoters_V2.bed")
names(proms)<-c("chr","start","end")

proms_gr <- GenomicRanges::makeGRangesFromDataFrame(proms)

overlaps <- findOverlaps(proms_gr, hic_regions_mm10)

# Extract the indices of overlapping Hi-C regions
hic_indices <- subjectHits(overlaps)

# Add a metadata column to the enhancers GRanges object
# e.g., Adding the names of the overlapping Hi-C regions
mcols(proms_gr)$hic_region <- NA  # Initialize column
mcols(proms_gr)$hic_region[queryHits(overlaps)] <- (hic_regions$name)[hic_indices]


as.numeric(enh_regions_gr$hic_region %>% unique) -> enh_hics
as.numeric(proms_gr$hic_region %>% unique) -> prom_hics
enh_hics[!is.na(enh_hics)] %>% sort -> enh_hics
prom_hics[!is.na(prom_hics)] %>% sort -> prom_hics

enh_hics2 <- paste0("bin_",enh_hics)
prom_hics2 <- paste0("bin_",prom_hics)

hic_sparse[enh_hics2,prom_hics2] -> hic_sparse_sml

save(enh_regions_gr,proms_gr,hic_sparse_sml,file = "~/GRanalysis_master/ABC/hic_data.RData")

data %>% filter(V1 %in% enh_hics) %>% filter(V2 %in% prom_hics) -> data_f

save(enh_regions_gr,proms_gr,data_f,file = "~/GRanalysis_master/ABC/hic_data_nonsparse.RData")
