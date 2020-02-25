############################### Pancancer study (Hoadley et al. 2014) ##############################
######################################### RPPA data kernel #########################################

rm(list=ls())

library(coca) # for consensusCluster()
library(klic)
library(pheatmap)

################################ Load copy number preprocessed data ################################

load("data/preprocessed-RPPA-data.RData")

######################################## Load sample names #########################################

load("data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

########################################## Select samples ##########################################

# To check which part of the sample names we need
rownames(RPPA)[which((grepl("BT-A20W", rownames(RPPA))))]

# Select only those that are present in Hoadley et al.'s COCA clustering
which_ones <- rownames(RPPA)[rownames(RPPA)%in%samples]

RPPA <- RPPA[which_ones, ]
RPPA_dist <- RPPA_dist[which_ones, which_ones]

save(which_ones, file = "data/names-RPPA.RData")

######################################## Generate kernel ###########################################

# cc <- consensusCluster(K = 8, B = 100, clMethod = "hc", dist = RPPA_dist, hcMethod = "ward.D")
load("data/kernel-rppa.RData")
coph_corr <- copheneticCorrelation(cc)

myBlues <-
    colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)

# png("figures/cc-RPPA.png", height = 800, width = 800)
# pheatmap(cc, annotations = annotations_COCA, color = myBlues, fontsize = 20)
# dev.off()

save(cc, file = "data/kernel-rppa.RData")
