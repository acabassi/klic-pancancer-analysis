############################
### Pancancer12 clusters ###
############################

rm(list=ls())

library(pheatmap)

#mRNAseq
mRNAseq <- read.table("PanCan12.3602-corrected-v3.Subtypes.K16.txt", header = T, sep = "\t", row.names = 1)
head(mRNAseq)
rownames(mRNAseq) <- substr(rownames(mRNAseq), 1, 12)
length(unique(rownames(mRNAseq)))
length(unique(substr(rownames(mRNAseq),1,12)))

#miRNAseq
miRNAseq <- read.table("miRNAk15.txt", header = TRUE, sep = "\t")
sum(1 - miRNAseq$Cluster==miRNAseq$Cluster.1) # the two clusters are the same
miRNAseq$Cluster[c(423, 424, 425, 426)] # same
miRNAseq$X[c(423, 424, 425, 426)] # same
miRNAseq$Sample[c(423, 424, 425, 426)] # different

#SCNA
SCNA <- read.table("SCNA_cluster_table.txt", header = TRUE, row.names = 1)
head(SCNA)
rownames(SCNA) <- substr(rownames(SCNA),1,12)

#Methylation
methylation <- read.csv("DNA_Methylation_Cluster_130519.csv", header = TRUE, row.names = 1)
head(methylation)
rownames(methylation) <- substr(rownames(methylation),1,12)

# RPPA
RPPA <- read.csv("PanCan11_RBN_SimpleCluster_20130411.csv", header = TRUE, row.names = 1)
head(RPPA)
rownames(RPPA) <- gsub(".", "-", rownames(RPPA), fixed = TRUE)

# Mutation
mutation <- read.table("Pancan_12_mutation_subtypes.tsv", sep = '\t', header = TRUE, row.names = 1)
head(mutation)

# COCA clusters
coca <- read.table("data/CofC.noMut.K13.Hoadley.20130523.txt", header = TRUE, row.names = 1)
head(coca)
length(unique(substr(rownames(coca),1,12)))
unique(substr(rownames(coca),13,15))
length(rownames(coca))
