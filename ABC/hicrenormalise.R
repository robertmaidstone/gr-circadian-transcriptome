library(tidyverse)
library(Matrix)
# Load the data into R
data <- read.table("~/GRanalysis_master/ABC/GSE104129_ZT10_5000_iced.matrix/5PM_5000_iced.matrix", header = FALSE)

# Determine the size of the matrix
names(data) <- c("bini","binj","interaction")

##
max_bin <- max(data$bini, data$binj)

hic_sparse <- sparseMatrix(i = data$bini, j = data$binj, x = data$interaction, dims = c(max_bin, max_bin))
# # Initialize the matrix
# hic_matrix <- matrix(0, nrow = max_bin, ncol = max_bin)
# # 
# # # Fill the matrix with interaction values
# for (i in 1:nrow(data)) {
#   hic_matrix[data$bini[i], data$binj[i]] <- data$interaction[i]
#   hic_matrix[data$binj[i], data$bini[i]] <- data$interaction[i]  # Ensure symmetry
# }

# Install and load the HiTC package if needed
# BiocManager::install("HiTC")
library(HiTC)

# Convert the matrix into Hi-C format

# hic_sparse <- as(Matrix(hic_matrix, sparse = TRUE), "dgCMatrix")
colnames(hic_sparse) <- rownames(hic_sparse) <- paste0("bin_", 1:max_bin)

granges_object <- rtracklayer::import("~/GRanalysis_master/ABC/GSE104129_mm9_5000_ord.bed/5AM_5000_ord-2.bed", format = "BED")
head(granges_object)

HTCdata<-HTCexp(intdata =hic_sparse,xgi = granges_object, ygi = granges_object)

# KR normalization
kr_norm <- HiTC::normICE(hic_sparse, method = "KR")
normalized_matrix <- as.matrix(kr_norm)

enhancer_bin <- 100  # Replace with actual enhancer bin
promoter_bin <- 150  # Replace with actual promoter bin

kr_contact <- normalized_matrix[enhancer_bin, promoter_bin]
print(paste("KR-normalized interaction frequency:", kr_contact))
