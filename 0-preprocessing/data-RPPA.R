############################### Pancancer study (Hoadley et al. 2014) ##############################
######################### Reverse phase protein array data preprocessing ###########################

rm(list=ls())

library(pheatmap)
library(mclust)

## From Hoadley et al. (2014)
# We performed unsupervised clustering on the protein expression data. Pearson correlation was used 
# as the distance metric and Ward was used as the linkage algorithm. We identified eight robust
# clusters. Theeight clusters and their protein expression patterns are shown in Figure S1E. As
# expected, most of the clusters are driven by tumor type. A few notable exceptions include basal
# and Her2 breast samples, which don’t cluster near the luminal breast samples; bladder (BLCA)
# samples, which cluster mainly with breast basal and Her2 samples; and head and neck (HNSC),
# lung squamous (LUSC) and lung adenocarcinoma (LUAD) samples that fall into a single cluster.
# Colon (COAD) and rectal (READ) samples cluster together,indicating that their proteomic profiles
# are very similar. The RPPA cluster memberships have been used for downstream analysis, such as the
# Cluster of Clusters Analysis in the main text and are discussed further there.

############################################# Load data ############################################

# The data can be downloaded from: https://www.synapse.org/Portal.html#!Synapse:syn1759392 
# The clusters can be downloeaded from: https://www.synapse.org/Portal.html#!Synapse:syn1756922.

RPPA <- read.csv("data-TCGA/PanCan11_RBN_RPPA_without_Duplicates_20130325.csv", 
                 header = TRUE, row.names = 1)
clusters <- read.csv("data-TCGA/PanCan11_RBN_SimpleCluster_20130411.csv",
                     header = TRUE, row.names = 1)
rownames(clusters) <- gsub(".", "-", rownames(clusters), fixed = TRUE)

# RPPA <- t(RPPA)
dim(RPPA) # 136 variables X 3467 samples 
# RPPA <- t(RPPA)
rownames_RPPA <- rownames(RPPA)
colnames_RPPA <- colnames(RPPA)

# Remove the first five columns 
RPPA_new <- matrix(unlist(RPPA[,-c(1:5)]), ncol = 131)
rownames(RPPA_new) <- rownames_RPPA
colnames(RPPA_new) <- colnames_RPPA[-c(1:5)]
pheatmap(RPPA_new[1:5,1:4], cluster_rows = F, cluster_cols = F) # Just to check if it worked
pheatmap(as.matrix(RPPA[1:5,6:9]), cluster_rows = F, cluster_cols = F) # These two should be the same

table(colSums(is.na(RPPA_new)))
# RPPA <- RPPA[,-which(colSums(is.na(RPPA))>0)]

# RPPA_scaled <- scale(RPPA_new, center = TRUE, scale = TRUE)

############################################## Cluster #############################################
RPPA_cor <- cor(t(RPPA_new), method = "pearson")
dim(RPPA_cor)
sum(colSums(is.na(RPPA_cor)))

RPPA_dist <- 1-RPPA_cor

save(RPPA, RPPA_dist, file = "data/preprocessed-RPPA-data.RData")

hc <- hclust(as.dist(RPPA_dist), method = "ward")
our_clusters <- cutree(hc, k = 8)

length(our_clusters)
length(clusters$K_8)

adjustedRandIndex(our_clusters, as.vector(clusters$K_8))

########################################### Plot clusters ##########################################

annotations <- data.frame(H = as.factor(as.vector(clusters$K_8)),
                          C = as.factor(our_clusters))

RPPA_sort <- RPPA_new[hc$order,]

negative_numbers <- linspace(min(RPPA_sort), 0, n = 16)
positive_numbers <-
    linspace(0, max(RPPA_sort), n = ceil((
        max(RPPA_sort) / (negative_numbers[2] - negative_numbers[1])
    )))
col_breaks <-
    c(negative_numbers, positive_numbers[2:length(positive_numbers)])
my_colours <-
    colorRampPalette(c("#FF9900", "white"))(length(negative_numbers))
my_colours <-
    c(my_colours, colorRampPalette(c("white", "#146EB4"))(length(positive_numbers)))

# Plot heatmap
# anno_col_RPPA <- anno_col[,c(1,2,5)]
# anno_col_RPPA$RPPA <- as.factor(anno_col_RPPA$RPPA)

png("figures/RPPA-hc.png", height = 650, width = 100 + 131 * 10)
pheatmap(
    RPPA_sort[,sample(1:(dim(RPPA_sort)[2]), 100, replace = FALSE)],
    color = my_colours,
    breaks = col_breaks,
    cluster_rows = F,
    annotation_row = annotations,
    show_rownames = F,
    show_colnames = F,
    fontsize = 20,
    cellwidth = 8
)
dev.off()

save(annotations, file = "clusters-RPPA.RData")
 