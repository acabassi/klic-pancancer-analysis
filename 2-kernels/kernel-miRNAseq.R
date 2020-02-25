################################ Pancancer study (Hoadley et al. 2014) ################################
######################################### miRNAseq data kernel ########################################

rm(list=ls())

library(coca) # for consensusCluster()
library(pheatmap)

################################# Load miRNAseq preprocessed data ##################################

load("data/preprocessed-miRNA-data.RData")
miRNA <- final_miRNA_data; rm(final_miRNA_data)

######################################### Load sample names ###########################################

load("data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

########################################### Select samples ############################################

# To check which part of the sample names we need
rownames(miRNA)[which((grepl("BT-A20W", rownames(miRNA))))]

# Check whether there are any duplicated sample names
sum(duplicated(substr(rownames(miRNA), 6, 17))) # None

miRNA <- miRNA[-which(duplicated(substr(rownames(miRNA), 6, 17))),]
miRNA_dist <- as.matrix(miRNA_dist)
miRNA_dist <- miRNA_dist[-which(duplicated(substr(rownames(miRNA_dist), 6, 17))), 
                         -which(duplicated(substr(rownames(miRNA_dist), 6, 17)))]

# Take only part of sample names that we need 
rownames(miRNA) <- substr(rownames(miRNA), 6, 17)
rownames(miRNA_dist) <- colnames(miRNA_dist) <- substr(rownames(miRNA_dist), 6, 17)

# Select only those that are present in Hoadley et al.'s COCA clustering
which_ones <- rownames(miRNA)[rownames(miRNA)%in%samples]

miRNA <- miRNA[which_ones, ]
miRNA_dist <- miRNA_dist[which_ones, which_ones]

save(which_ones, file = "data/names-miRNA.RData")

######################################### Generate kernel #############################################

cc <- consensusCluster(K = 7, B = 100, clMethod = "hc", dist = miRNA_dist, hcMethod = "average")

png("figures/cc-miRNA.png")
pheatmap(cc, annotations = annotations_COCA)
dev.off()

rownames(cc) <- colnames(cc) <- which_ones

save(cc, file = "data/kernel-miRNA.RData")
