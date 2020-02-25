############################### Pancancer study (Hoadley et al. 2014) ##############################
######################################## mRNAseq data kernel #######################################

rm(list=ls())

library(coca) # for consensusCluster()
library(klic)
library(pheatmap)

################################# Load mRNAseq preprocessed data ###################################

load("data/preprocessed-mRNA-data.RData")

######################################## Load sample names #########################################

load("data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

########################################## Select samples ##########################################

# To check which part of the sample names we need
rownames(mRNAseq)[which((grepl("BT-A20W", rownames(mRNAseq))))]

# Check whether there are any duplicated sample names
sum(duplicated(substr(rownames(mRNAseq), 1, 12))) # None

mRNAseq <- mRNAseq[-which(duplicated(substr(rownames(mRNAseq), 1, 12))),]

# Take only part of sample names that we need 
rownames(mRNAseq) <- substr(rownames(mRNAseq), 1, 12)

# Select only those that are present in Hoadley et al.'s COCA clustering
which_ones <- rownames(mRNAseq)[rownames(mRNAseq)%in%samples]

mRNAseq <- mRNAseq[which_ones, ]

save(which_ones, file = "data/names-mRNA.RData")

######################################## Generate kernel ###########################################

# cc <- consensusCluster(mRNAseq, K = 16, B = 100, clMethod = "hc", hcMethod = "average")
load("data/kernel-mRNA.RData")
coph_corr <- copheneticCorrelation(cc)

myBlues <-
    colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)

# png("figures/cc-mRNA.png", height = 800, width = 800)
# pheatmap(cc, annotations = annotations_COCA, color = myBlues, fontsize = 20)
# dev.off()

rownames(cc) <- colnames(cc) <- which_ones

save(cc, file = "data/kernel-mRNA.RData")
