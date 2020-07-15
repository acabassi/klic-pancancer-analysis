################################ Pancancer study (Hoadley et al. 2014) #############################
################################### Copy number data preprocessing #################################

rm(list = ls())
setwd("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-DATA/data-PANCAN12")

library(circlize)
library(mclust)
library(pheatmap)
library(pracma)

## From Hoadley et al. (2014)

# For copy number based clustering, tumors were clustered based on thresholdedcopy number at
# reoccurring alteration peaks from GISTIC analysis. Tumors were hierarchical clustered in R based
# on Euclidean distance using Ward’s method. The number of cluster groups was chosen based on
# cophenetic distancesgenerated from clustering. Forcomparison of broad and focal alteration between
# cluster of cluster groups, frequency of alterations in each cluster group was compared to the
# average frequency of all other groups by chi squared tests with an added Bonferroni correction to
# control for multiple testings.See Figures S1CandS4A-C

# The data can be downloaded from: https://www.synapse.org/Portal.html#!Synapse:syn1710678
# The clusters can be downloaded from: https://www.synapse.org/Portal.html#!Synapse:syn1712142.

###################################### Load copy number clusters ###################################

clusters <-
    read.table(
        "data-TCGA/SCNA_Cluster_table.txt",
        sep = "\t",
        header = TRUE,
        row.names = 1
    )
rownames(clusters)
dim(clusters)
clusters8 <- clusters$named_k8_clusters
names(clusters8) <- rownames(clusters)

######################################## Load copy number data #####################################

SCNA <-
    read.table("data-TCGA/all_lesions.conf_95.pancan12.txt",
               sep = "\t",
               header = TRUE)
rownames(SCNA)
dim(SCNA) # 168 variables X 4944 samples
SCNA <- t(SCNA)
rownames(SCNA)
rownames(SCNA) <- gsub(".", "-", rownames(SCNA), fixed = TRUE)
sum(1 - rownames(SCNA) %in% rownames(clusters)) # There are 10 extra rows

# rownames(SCNA) <-  substring(rownames(SCNA), 1, 12)
# Let's leave this for later

# First 9 ROWS are:
# [1] "Unique.Name" [2] "Descriptor" [3] "Wide.Peak.Limits" [4] "Peak.Limits"
# [5] "Region.Limits" [6] "q.values" [8] "Broad.or.Focal" [9] "Amplitude.Threshold"
# [7] "Residual.q.values.after.removing.segments.shared.with.higher.peaks"
# And the last row contains only NAs
SCNA[4944,]

colnames(SCNA) <- SCNA[2,]
# We only need the second half of the matrix
sum(1 - colnames(SCNA[, 1:84]) %in% colnames(SCNA[, 85:168]))

SCNA <- SCNA[10:4943, 85:168]

class(SCNA) <- "numeric"
dim(SCNA)

# Data scaling?

colMeans(SCNA)
apply(SCNA, 2, var)

min(SCNA)
max(SCNA)

# Scaling is not mentioned in the paper by Hoadley et al. but, since it improves
# the correspondence between our clusterings, we do it here.
SCNAscaled <- scale(SCNA, center = TRUE, scale = TRUE) 

sum(1 - rownames(SCNA) %in% rownames(clusters)) # Now the row names are the same

# TEMPORARILY Remove duplicate name (should add A and D instead)
colnames(SCNA) <- gsub(" ", "", colnames(SCNA), fixed = TRUE)
colnames(SCNA)[which(colnames(SCNA) == "4p16.3")[2]] <- "4p16.3 2"

save(SCNAscaled, file = "data/preprocessed-SCNA-data.RData")

####################################### Hierarchical clustering ####################################

hc <-
    hclust(dist(SCNAscaled, method = "euclidean"), method = "ward.D")
ourClusters <- cutree(hc, k = 8)
clusters8 <- clusters8[names(ourClusters)]
adjustedRandIndex(as.factor(ourClusters), clusters8)
sum(1 - (names(ourClusters) == names(clusters8[names(ourClusters)])))

############################################# Plot heatmap #########################################

# Choose row order.
# Option 1: use order given by hierarchical clustering
rowOrder <- hc$order

# Option 2: sort by TCGA cluster
rowOrder <- c(
    names(clusters8)[which(clusters8 == "Iq")],
    names(clusters8)[which(clusters8 == "BRCA-LUAD+")],
    names(clusters8)[which(clusters8 == "COAD-READ")],
    names(clusters8)[which(clusters8 == "Squamous")],
    names(clusters8)[which(clusters8 == "High")],
    names(clusters8)[which(clusters8 == "Quiet")],
    names(clusters8)[which(clusters8 == "GBM")],
    names(clusters8)[which(clusters8 == "Kirc+")]
)

SCNA <- SCNA[rowOrder,]

annotations <-
    data.frame(H = as.factor(clusters8),
               C = as.factor(ourClusters))

negative_numbers <- linspace(min(SCNA), 0, n = 8)
positive_numbers <-
    linspace(0, max(SCNA), n = ceil((
        max(SCNA) / (negative_numbers[2] - negative_numbers[1])
    )))
col_breaks <-
    c(negative_numbers, positive_numbers[2:length(positive_numbers)])
my_colours <-
    colorRampPalette(c("#FF9900", "white"))(length(negative_numbers))
my_colours <-
    c(my_colours, colorRampPalette(c("white", "#146EB4"))(length(positive_numbers)))

png("figures/SCNA-hc.png",
    height = 650,
    width = 100 + 84 * 10)
pheatmap(
    SCNA,
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

save(annotations, file = "data/clusters-CN.RData")

################################# Bonus: plot using ComplexHeatmap #################################

# library(ComplexHeatmap)

# col_fun = colorRamp2(c(min(SCNA), 0, max(SCNA)), c("#146EB4", "white",  "#FF9900"))
# col_fun(seq(-3, 3))
#
# hoadley_palette <-
#     c(
#         "Squamous" = "#999999",
#         "Kirc+" = "#E69F00",
#         "High" = "#56B4E9",
#         "COAD-READ" = "#009E73",
#         "BRCA-LUAD+" = "#F0E442",
#         "GBM" = "#0072B2",
#         "Quiet" = "#D55E00",
#         "Iq" = "#CC79A7"
#     )
#
# library(pals)
# pal.bands(alphabet, alphabet2, cols25, glasbey, kelly, polychrome,
#           stepped, tol, watlington,
#           show.names=FALSE)
#
# clusters <- HeatmapAnnotation("Hoadley et al." = as.factor(clusters8),
#                               "Cabassi & Kirk" = as.factor(ourClusters),
#                               which = "row",
#                               col = list("Hoadley et al." = hoadley_palette))
# Heatmap(SCNA,
#         cluster_rows = F,
#         color = col_fun,
#         show_row_names = F,
#         right_annotation = clusters,
#         heatmap_legend_param = list(title = "Values"))
