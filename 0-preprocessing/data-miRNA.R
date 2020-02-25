############################### Pancancer study (Hoadley et al. 2014) ##############################
##################################### miRNA data preprocessing #####################################

rm(list = ls())
setwd("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-DATA/data-PANCAN12")

library(cluster) 
library(data.table)
library(mclust)
library(pheatmap)
library(pracma)
library(stringr)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("miRBaseConverter")
library(miRBaseConverter)

## From Hoadley et al. (2014)

# Using Cluster 3 [1] we log-transformed and median-centered the 51 miRNA abundance profiles, and
# then hierarchically clustered only the rows using an absolute centered correlation and average
# linkage. We visualized the resulting matrix with Java Treeview [2].

# [1] de Hoon, M.J., et al., Open source clustering software.Bioinformatics, 2004. 20(9): p. 1453-4
# [2] Saldanha, A.J., Java Treeview--extensible visualization of microarray data.Bioinformatics,
#     2004. 20(17): p. 3246-8.

############################################# Load data ############################################

# The data can be downloaded from: https://www.synapse.org/Portal.html#!Synapse:syn2491366
# The clusters can be downloaded from: https://www.synapse.org/Portal.html#!Synapse:syn2027079

miRNAseqSubtypes <-
  fread("data-TCGA/miRNA.k15.txt", header = TRUE, sep = "\t")

miRNAseq <- read.csv(
  "data-TCGA/PanCan.miRNAseq.RPM.215-MIMATs-most-variant-25pc.4229-samples.NMF-input.BCGSC.20140603.csv",
  head = TRUE,
  row.names = 1
)

dim(miRNAseq) # 215 variables X 4229 samples
miRNAseq <- t(miRNAseq)

rownames(miRNAseq) <-
  gsub(".", "-", rownames(miRNAseq), fixed = TRUE)
# Don't shorten rownames, as they become non-unique

################## Select 51 miRNA abundance profiles used by Hoadley et al (2014) #################

few_miRNAs <- read.csv("data-TCGA/51mirnas.csv", head = FALSE)
few_miRNAs <- unlist(few_miRNAs, use.names = FALSE)

accessions <- str_sub(colnames(miRNAseq), start = -12)
names <- miRNA_AccessionToName(accessions, targetVersion = "v22")
colnames(miRNAseq) <- gsub("hsa-", "", names$TargetName)

few_miRNAs %in% colnames(miRNAseq)
not_found_indices <- which(!few_miRNAs %in% colnames(miRNAseq))
not_found <- few_miRNAs[not_found_indices]
for (i in not_found) {
  corresponding_name <-
    which(grepl(i, colnames(miRNAseq), fixed = TRUE))
  colnames(miRNAseq)[corresponding_name] <-
    substr(colnames(miRNAseq)[corresponding_name],
           1, nchar(colnames(miRNAseq)[corresponding_name]) -
             3)
}

miRNAseq <- miRNAseq[, few_miRNAs]

########################################### Data scaling ###########################################

# Log-transform and median-center data
miRNAseqLog2 <- log2(miRNAseq + 1)
miRNAseqLog2MedianCentred <- miRNAseqLog2
for (i in 1:ncol(miRNAseqLog2MedianCentred)) {
  miRNAseqLog2MedianCentred[, i] <- miRNAseqLog2MedianCentred[, i] -
    median(miRNAseqLog2MedianCentred[, i])
}

myCorMat <- cor(t(miRNAseqLog2MedianCentred))

subtypes   <- miRNAseqSubtypes$Cluster
#shortnames <- sapply( miRNAseqSubtypes$Sample, substr, 1, 12)
cancerNames <- miRNAseqSubtypes$Disease_code
cancerNames[cancerNames == "OV"] <- "OVCA"
names(subtypes) <-
  paste(cancerNames, miRNAseqSubtypes$Sample, sep = "_")
names(subtypes)[which(names(subtypes) == "OVCA_TCGA-09-0366-01A-01R")] <-
  c(
    "OVCA_TCGA-09-0366-01A-01R",
    "OVCA_TCGA-09-0366-01A-01R-1",
    "OVCA_TCGA-09-0366-01A-01R-2",
    "OVCA_TCGA-09-0366-01A-01R-3"
  )
inds <- match(names(subtypes), rownames(myCorMat))

annotation_row2 <-
  data.frame(Cluster = as.factor(miRNAseqSubtypes$Cluster))
rownames(annotation_row2) <- names(subtypes)

final_miRNA_data <- miRNAseqLog2MedianCentred
# R is smart enough not to make a copy if the variable is the same

############################################## Cluster #############################################

# Compute all the pairwise dissimilarities (distances) between observations in the data set
miRNA_dist <- myDist <- daisy(myCorMat)
save(miRNA_dist, final_miRNA_data, file = "data/preprocessed-miRNA-data.RData")

# Compute agglomerative hierarchical clustering of the dataset
ar <- agnes(miRNA_dist)

average_silh <- seq(0, 0, length = 20)
for (kval in 2:20) {
  si <- silhouette(cutree(ar, k = kval), miRNA_dist)
  summ <- summary(si)
  average_silh[kval] <- summ$avg.width
}

png("figures/miRNAseq_silhouette.png",
    width = 700,
    height = 500)
plot(
  2:20,
  average_silh[2:20],
  type = "b",
  xlab = "Clusters",
  ylab = "Average silhouette",
  cex = 1.2,
  cex.lab = 1.5,
  cex.axis = 1.2
)
dev.off()

print(which.max(average_silh))

finalClustering <- cutree(ar, k = which.max(average_silh))

TCGAclusters <- miRNAseqSubtypes$Cluster
names(TCGAclusters) <- miRNAseqSubtypes$V4

rwCrMat <-  rownames(myCorMat)
rwCrMat <-
  names(finalClustering) <- gsub("[.]", "-", substr(rwCrMat, 6, 20))
dupl <- which(duplicated(substr(rwCrMat, 6, 20)))
rwCrMat <- rwCrMat[-dupl]

finalClustering <- finalClustering[-dupl]

########################################### Plot clusters ##########################################

annotations <- data.frame(H = as.factor(finalClustering),
                          C = as.factor(as.vector(TCGAclusters[names(finalClustering)])))
rownames(annotations) <- rwCrMat
annotations <- annotations[order(annotations$H),]

adjustedRandIndex(annotations$C, annotations$H)

rownames(final_miRNA_data) <- substr(rownames(final_miRNA_data), 6, 20)
remove <- which(duplicated(rownames(final_miRNA_data)))
final_miRNA_data <- final_miRNA_data[-remove, ]

# rownames(final_miRNA_data) <- gsub("[.]", "-", rownames(final_miRNA_data))

negative_numbers <- linspace(min(final_miRNA_data, na.rm = TRUE), 0, n = 16)
positive_numbers <-
  linspace(0, max(final_miRNA_data, na.rm = TRUE), n = ceil((
    max(final_miRNA_data, na.rm = TRUE) / (negative_numbers[2] - negative_numbers[1])
  )))
col_breaks <-
  c(negative_numbers, positive_numbers[2:length(positive_numbers)])
my_colours <-
  colorRampPalette(c("#FF9900", "white"))(length(negative_numbers))
my_colours <-
  c(my_colours, colorRampPalette(c("white", "#146EB4"))(length(positive_numbers)))

png("figures/miRNAseq-hc.png",
    height = 100 + dim(final_miRNA_data)[1] / 6,
    width = 150 + 51 * 10)
pheatmap(
  final_miRNA_data[rownames(annotations),],
  color = my_colours,
  breaks = col_breaks,
  cluster_rows = F,
  annotation_row = annotations,
  show_colnames = F,
  show_rownames = F,
  fontsize = 20,
  cellwidth = 8
)
dev.off()

save(annotations, file = "data/clusters-miRNA.RData")

# These clusters do not match very well the clusters of Hoadley et al. (2014)
# Let's see if they match the tissue types instead

load("data/samples.RData")

################ Bonus: check adusted Rand index given by other numbers of clusters ################

for(kval in 2:20) {
  clustering <- cutree(ar, k = kval) # Use 15 clusts, as in paper
  names(clustering) <- rownames(myCorMat)
  subtypes   <- miRNAseqSubtypes$Cluster
  #shortnames <- sapply( miRNAseqSubtypes$Sample, substr, 1, 12)
  names(subtypes) <-
    paste(cancerNames, miRNAseqSubtypes$Sample, sep = "_")
  inds <- match(names(subtypes), names(clustering))
  clustering <- clustering[inds]
  all(names(clustering) == names(subtypes))
  
  print(kval)
  print(mclust::adjustedRandIndex(clustering, subtypes))
}

clusters15 <- cutree(ar, k = 15)

annotation_row2 <- data.frame(clusters = as.factor(clusters15))
rownames(annotation_row2) <- rownames(myCorMat)

# png(file = "clusters15.png")
# pheatmap(
#   myCorMat[sort.int(clusters15, index.return = T)$ix, sort.int(clusters15, index.return = T)$ix],
#   show_rownames = F,
#   show_colnames = F,
#   cluster_rows = F,
#   cluster_cols = F,
#   annotation_row = annotation_row2
# )
# dev.off()

########################## Bonus: Check ARI between tissues and clusters ###########################

annotations_full <- annotations
duplicates <- which(duplicated(substr(rownames(annotations_full), 1, 12)))
annotations_full <- annotations_full[-duplicates,]
rownames(annotations_full) <- substr(rownames(annotations_full), 1, 12)
annotations_full$Tissue <- anno_col[rownames(annotations_full),]$Tissue
annotations_full$COCA <- anno_col[rownames(annotations_full),]$COCA
adjustedRandIndex(annotations_full$C, annotations_full$Tissue) # 0.2476742
adjustedRandIndex(annotations_full$H, annotations_full$Tissue) # 0.5737194
adjustedRandIndex(annotations_full$C, annotations_full$COCA) # 0.2807347
adjustedRandIndex(annotations_full$H, annotations_full$COCA) # 0.6254234