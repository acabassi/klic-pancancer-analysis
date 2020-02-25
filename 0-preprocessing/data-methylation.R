################################ Pancancer study (Hoadley et al. 2014) ################################
################################### Methylation data preprocessing ####################################

rm(list = ls())

library(pheatmap)
library(mclust)

## From Hoadley et al. (2014)
# For the clustering analysis, we focused on CpG loci that are unmethylated in normal tissues.
# Therefore, we removed probes that showed methylation (median beta value > 0.2) in any of the 12
# matched normal tissue types included in the current study. After the aforementioned filters,
# 11,696 probes remained. As these loci are mostly within CpG islands that remain constitutively
# unmethylated in normal tissues, we dichotomized the beta values in the tumors at 0.3. Tumors with
# a beta value of 0.3 or greater are designated methylated and tumors with a beta value of lower
# than 0.3 are designated unmethylated. The dichotomization greatly ameliorated the effect of tumor
# sample purity on the clustering, and further removed most residual batch/platform effects that are
# primarily reflected in small variations near the two ends. We selected the 2,203 probes that were
# methylated in more than 10% of any of the tumor types or 50% of any of the well-defined subtypes
# for clustering. We used hierarchical clustering with Ward’s method on the Jaccard Distance, a
# distance measure that best suits binary data. The dendrogram was cut at different levels with the
# ‘cutree’ function in R and evaluated for associations with clinical data and the k=19 result was
# used.See FigureS1D.

############################################# Load data ############################################
# The data can be downloaded from https://www.synapse.org/Portal.html#!Synapse:syn2486658
# The clusters can be downloaded from https://www.synapse.org/Portal.html#!Synapse:syn1875816.

methylation <-
    read.csv(
        "data-TCGA/DNAmethylationClusteringMatrix.csv",
        header = TRUE,
        row.names = 1
    )
dim(methylation) # 2043 variables X 4923 samples
methylation <- t(methylation)
sum(duplicated(colnames(methylation)))
rownames(methylation) <-
    gsub(".", "-", rownames(methylation), fixed = TRUE)

# I am not going to remove the extra samples because in the plot they use them all
# methylation <- methylation[rownames(methylation)%in%rownames(anno_col),]

# Hierarchical clustering
methylation_dist <-
    proxy::dist(methylation, by_rows = TRUE, method = "Jaccard")
dim(as.matrix(methylation_dist))
sum(colSums(is.na(as.matrix(methylation_dist))))

save(methylation, methylation_dist, file = "data/preprocessed-methylation-data.RData")

hc <- hclust(methylation_dist, method = "ward.D")
our_clusters <- cutree(hc, k = 19)
rowOrder <- hc$order
methylation <- methylation[rowOrder,]

# Load clusters
meth_clusters <-
    read.csv(
        "data-TCGA/DNA_Methylation_Cluster_130519.csv",
        header = TRUE,
        row.names = 1
    )
duplicates <- which(duplicated(meth_clusters[, 1]))
meth_clusters[duplicates, 1]
# lots of duplicates :(
for (dup in duplicates) {
    others <- which(meth_clusters[, 1] == meth_clusters[dup, 1])
    for (i in 1:length(others)) {
        print(i)
        print(meth_clusters[dup, 2])
        print(meth_clusters[others[i], 2])
        #     if(meth_clusters[dup,2] != meth_clusters[others[i],2])
        #         stop('Same name, different clusters')
    }
}
# but they are all NAs, so we don't care
meth_clusters <- meth_clusters[-duplicates,]
rownames(meth_clusters) <- meth_clusters[, 1]

# check that there is the same number of elements in our clusters and in the given clusters
dim(meth_clusters)[1] - sum(is.na(meth_clusters[, 2])) == length(our_clusters)

names_meth_clusters <- meth_clusters[, 1]
meth_clusters <- as.vector(meth_clusters[, 2])
names(meth_clusters) <- names_meth_clusters

names(our_clusters) %in% names(meth_clusters)
meth_clusters <- meth_clusters[names(our_clusters)]

# Adjusted Rand Index
adjustedRandIndex(our_clusters, meth_clusters)

annotations <- data.frame(H = as.factor(meth_clusters),
                          C = as.factor(our_clusters))

my_colours <- c("white",  "#146EB4")
col_breaks <- c(0, 0.5, 1)

png(
    "figures/methylation-hc.png",
    height = 100 + dim(methylation)[1] / 4,
    width = 100 + 100 * 10
)
pheatmap(
    methylation[, 1:100],
    color = my_colours,
    breaks = col_breaks,
    legend_breaks = col_breaks,
    cluster_cols = T,
    cluster_rows = F,
    annotation_row = annotations,
    show_colnames = F,
    show_rownames = F,
    fontsize = 20,
    cellwidth = 8,
    cellheight = .25
)
dev.off()

save(annotations, file = "data/clusters-methylation.RData")