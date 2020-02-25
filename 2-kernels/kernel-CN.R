############################### Pancancer study (Hoadley et al. 2014) ##############################
###################################### Copy number data kernel #####################################

rm(list=ls())

library(coca) # for consensusCluster()
library(klic)
library(pheatmap)

################################ Load copy number preprocessed data ################################

load("data/preprocessed-SCNA-data.RData")

######################################## Load sample names #########################################

load("data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

########################################## Select samples ##########################################

# To check which part of the sample names we need
rownames(SCNAscaled)[which((grepl("BT-A20W", rownames(SCNAscaled))))]

# Check whether there are any duplicated sample names
which(duplicated(substr(rownames(SCNAscaled), 1, 12))) # None

# Take only part of sample names that we need 
rownames(SCNAscaled) <- substr(rownames(SCNAscaled), 1, 12)

# Select only those that are present in Hoadley et al.'s COCA clustering
which_ones <- rownames(SCNAscaled)[rownames(SCNAscaled)%in%samples]

SCNAscaled <- SCNAscaled[which_ones, ]

save(which_ones, file = "data/names-CN.RData")

######################################## Generate kernel ###########################################

# cc <- consensusCluster(SCNAscaled, K = 8, B = 1000, clMethod = "hc", hcMethod = "ward.D")
load("data/kernel-cn.RData")

coph_corr <- copheneticCorrelation(cc)

myBlues <-
    colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)

# png("figures/cc-CN.png", height = 800, width = 800)
# pheatmap(cc, annotations = annotations_COCA, color = myBlues, fontsize = 20)
# dev.off()

rownames(cc) <- colnames(cc) <- which_ones

save(cc, file = "data/kernel-cn.RData")
