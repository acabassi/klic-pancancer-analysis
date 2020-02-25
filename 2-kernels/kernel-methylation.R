############################### Pancancer study (Hoadley et al. 2014) ##############################
###################################### Methylation data kernel #####################################

rm(list=ls())

library(coca) # for consensusCluster()
library(klic)
library(pheatmap)

################################ Load methylation preprocessed data ################################

load("data/preprocessed-methylation-data.RData")

######################################## Load sample names #########################################

load("data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

########################################## Select samples ##########################################

# To check which part of the sample names we need
rownames(methylation)[which((grepl("BT-A20W", rownames(methylation))))]

# Check whether there are any duplicated sample names
sum(duplicated(substr(rownames(methylation), 1, 12))) # None

methylation <- methylation[-which(duplicated(substr(rownames(methylation), 1, 12))),]
methylation_dist <- as.matrix(methylation_dist)
methylation_dist <- methylation_dist[-which(duplicated(substr(rownames(methylation_dist), 1, 12))),
                                     -which(duplicated(substr(colnames(methylation_dist), 1, 12)))]

# Take only part of sample names that we need 
rownames(methylation) <- substr(rownames(methylation), 1, 12)
rownames(methylation_dist) <- colnames(methylation_dist) <- substr(rownames(methylation_dist), 1, 12)

# Select only those that are present in Hoadley et al.'s COCA clustering
which_ones <- rownames(methylation)[rownames(methylation)%in%samples]

methylation <- methylation[which_ones, ]
methylation_dist <- methylation_dist[which_ones, which_ones]

save(which_ones, file = "data/names-methylation.RData")

######################################## Generate kernel ###########################################

# cc <- consensusCluster(K = 19, B = 100, clMethod = "hc", dist = methylation_dist, hcMethod = "ward.D")
load("data/kernel-methylation.RData")
coph_corr <- copheneticCorrelation(cc)

rownames(cc) <- colnames(cc) <- rownames(methylation_dist)

myBlues <-
    colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)

# png("figures/cc-methylation.png", height = 800, width = 800)
# pheatmap(cc, annotations = annotations_COCA, color = myBlues, fontsize = 20, show_rownames = F,
#          show_colnames = F)
# dev.off()

rownames(cc) <- colnames(cc) <- which_ones

save(cc, file = "data/kernel-methylation.RData")
