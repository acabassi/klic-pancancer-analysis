################################ Pancancer study (Hoadley et al. 2014) #############################
###################################### Reproducing COCA clusters ###################################

rm(list=ls())

library(coca)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(mclust)

######################################### Load sample names ########################################

load("data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

######################################### Create COCA matrices #####################################

OURclusters <- matrix(NA, length(samples), 5)
TCGAclusters <- matrix(NA, length(samples), 5)
rownames(OURclusters) <- rownames(TCGAclusters) <- samples

###################################### Load copy number clusters ###################################

load("data/clusters-CN.RData")

# To check which part of the sample names we need
rownames(annotations)[which((grepl("BT-A20W", rownames(annotations))))]

# Check whether there are any duplicated sample names
which(duplicated(substr(rownames(annotations), 1, 12))) # None

# Take only part of sample names that we need 
rownames(annotations) <- substr(rownames(annotations), 1, 12)

# Select only those that are present in Hoadley et al.'s COCA clustering
which_ones <- rownames(annotations)[rownames(annotations)%in%samples]

# Save them
OURclusters[which_ones,1] <- annotations[which_ones,]$C
TCGAclusters[which_ones,1] <- annotations[which_ones,]$H
datasetNames <- c("CN")

###################################### Load methylation clusters ###################################

load("data/clusters-methylation.RData")

# To check which part of the sample names we need
rownames(annotations)[which((grepl("BT-A20W", rownames(annotations))))]

# Shorten names
short_names <- substr(rownames(annotations), 1, 12)

# Check whether there are any duplicated sample names
which(duplicated(short_names)) # Yes

# Remove duplicate names 
remove <- unique(short_names[which(duplicated(short_names))])

# Remove duplicated names
annotations <- annotations[!short_names%in%remove,]

# Take only part of sample names that we need 
rownames(annotations) <- short_names[!short_names%in%remove]

# Select only those that are present in Hoadley et al.'s COCA clustering
which_ones <- rownames(annotations)[rownames(annotations)%in%samples]

# Save them
OURclusters[which_ones,2] <- annotations[which_ones,]$clusters
TCGAclusters[which_ones,2] <- annotations[which_ones,]$TCGAclusters
datasetNames <- c(datasetNames, "Meth")

########################################## Load RPPA clusters ######################################

load("data/clusters-RPPA.RData")

# To check which part of the sample names we need
rownames(annotations)[which((grepl("BT-A20W", rownames(annotations))))]

# Select only those that are present in Hoadley et al.'s COCA clustering
which_ones <- rownames(annotations)[rownames(annotations)%in%samples]

# Save them
OURclusters[which_ones, 3] <- annotations[which_ones,]$clusters
TCGAclusters[which_ones, 3] <- annotations[which_ones,]$TCGAclusters
datasetNames <- c(datasetNames, "RPPA")

########################################## Load mRNA clusters ######################################

load("data/clusters-mRNAseq.RData")

# To check which part of the sample names we need
rownames(annotations)[which((grepl("BT-A20W", rownames(annotations))))]

# Shorten names
short_names <- substr(rownames(annotations), 1, 12)

# Check whether there are any duplicated sample names
which(duplicated(short_names)) # Yes

# Remove duplicate names 
remove <- unique(short_names[which(duplicated(short_names))])

# Remove duplicated names
annotations <- annotations[!short_names%in%remove,]

# Take only part of sample names that we need 
rownames(annotations) <- short_names[!short_names%in%remove]

# Select only those that are present in Hoadley et al.'s COCA clustering
which_ones <- rownames(annotations)[rownames(annotations)%in%samples]

# Save them
OURclusters[which_ones, 4] <- annotations[which_ones,]$C
TCGAclusters[which_ones, 4] <- annotations[which_ones,]$H
datasetNames <- c(datasetNames, "mRNA")

########################################## Load miRNA clusters #####################################

load("data/clusters-miRNA.RData")

# Check which part of the sample names we need
rownames(annotations)[which((grepl("BT-A20W", rownames(annotations))))]

# Shorten names
short_names <- substr(rownames(annotations), 1, 12)

# Check whether there are any duplicated sample names
which(duplicated(short_names)) # Yes

# Remove duplicate names 
remove <- unique(short_names[which(duplicated(short_names))])

# Remove duplicated names
annotations <- annotations[!short_names%in%remove,]

# Take only part of sample names that we need 
rownames(annotations) <- short_names[!short_names%in%remove]

# Select only those that are present in Hoadley et al.'s COCA clustering
which_ones <- rownames(annotations)[rownames(annotations)%in%samples]

# Save them
OURclusters[which_ones, 5] <- annotations[which_ones,]$C
TCGAclusters[which_ones, 5] <- annotations[which_ones,]$H
datasetNames <- c(datasetNames, "miRNA")

######################################## Plot matrix of clusters ###################################

# OUR clusters
OURexpandedMOC <- expandMOC(OURclusters, datasetNames = datasetNames)
attach(OURexpandedMOC)
plotMOC(moc, datasetIndicator, datasetNames, annotations = annotations_COCA, save = TRUE, 
        fileName = "figures/moc-OURS-NAs.png")

# TCGA clusters
TCGAexpandedMOC <- expandMOC(TCGAclusters, datasetNames = datasetNames)
attach(TCGAexpandedMOC)
plotMOC(moc, datasetIndicator, datasetNames, annotations = annotations_COCA, save = TRUE, 
        fileName = "figures/moc-TCGA-NAs.png")

################################################## COCA #############################################

# Extract MOC matrices
OURmoc <- OURexpandedMOC$moc
TCGAmoc <- TCGAexpandedMOC$moc

# Replace all NAs with zeroes 
OURmoc[is.na(OURmoc)] <- 0
TCGAmoc[is.na(TCGAmoc)] <- 0

# Apply coca() function from R package ``coca''
# OURcoca <- coca(OURmoc, maxK = 20, savePNG = TRUE, fileName = "OUR-coca", verbose = TRUE)
# TCGAcoca <- coca(TCGAmoc, maxK = 20, savePNG = TRUE, fileName = "TCGA-coca", verbose = TRUE)
# save(OURcoca, TCGAcoca, file = "data/coca-output.RData")
load("data/coca-output.RData")
 
OURcoca$K
TCGAcoca$K

# Apply ConsensusClusterPlus() function from R package ``ConsensusClusterPlus''
# OURcc <- ConsensusClusterPlus(t(OURmoc), maxK = 20, reps = 1000, title = "OURcc")
# TCGAcc <- ConsensusClusterPlus(t(TCGAmoc), maxK = 20, reps = 1000, title = "TCGAcc")
# save(OURcc, TCGAcc, file = "data/ConsensusClusterPlus-output.RData")
load("data/ConsensusClusterPlus-output.RData")

OURcoca_clusters <- OURcoca$clusterLabels 
names(OURcoca_clusters) <- rownames(OURmoc)
annotations_COCA$OURcoca <- factor("1", levels = c("1","2","3","4","5","6","7","8","9","10"))
annotations_COCA[names(OURcoca_clusters),]$OURcoca <- as.factor(OURcoca_clusters)
# annotations_COCA$OURcc <- 
#     factor("1", levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13"))
# annotations_COCA[names(OURcc[[13]]$consensusClass),]$OURcc <- as.factor(OURcc[[13]]$consensusClass) 
names(annotations_COCA) <- c("COCA - Hoadley et al.", "Tissue", "COCA clusters")
annotations_COCA <- annotations_COCA[,c(3, 1, 2)]
plotMOC(
    OURexpandedMOC$moc,
    OURexpandedMOC$datasetIndicator,
    OURexpandedMOC$datasetNames,
    annotations = annotations_COCA,
    save = TRUE,
    fileName = "figures/moc-OURS-coca.png"
)

TCGAcoca_clusters <- TCGAcoca$clusterLabels 
names(TCGAcoca_clusters) <- rownames(TCGAmoc)
annotations_COCA$TCGAcoca <- 
    factor("1", levels = c("1","2","3","4","5","6","7","8","9","10"))
annotations_COCA[names(TCGAcoca_clusters),]$TCGAcoca <- as.factor(TCGAcoca_clusters)
annotations_COCA$TCGAcc <- 
    factor("1", levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13"))
annotations_COCA[names(TCGAcc[[13]]$consensusClass),]$TCGAcc <- 
    as.factor(TCGAcc[[13]]$consensusClass) 
annotations_COCA$OURcoca <- NULL
annotations_COCA$OURcc <- NULL
plotMOC(TCGAexpandedMOC$moc, datasetIndicator, datasetNames, annotations = annotations_COCA, 
        save = TRUE, fileName = "figures/moc-TCGA-coca.png")

####################################### Adjusted Rand Index ########################################

# Add these back
annotations_COCA$OURcoca <- factor("1", levels = c("1","2","3","4","5","6","7","8","9","10"))
annotations_COCA[names(OURcoca_clusters),]$OURcoca <- as.factor(OURcoca_clusters)
annotations_COCA$OURcc <- 
    factor("1", levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13"))
annotations_COCA[names(OURcc[[13]]$consensusClass),]$OURcc <- as.factor(OURcc[[13]]$consensusClass) 

save(annotations_COCA, file = "data/coca-clusters.RData")

adjustedRandIndex(annotations_COCA$Tissue, annotations_COCA$COCA) # 0.8175554

adjustedRandIndex(annotations_COCA$Tissue, annotations_COCA$TCGAcoca) # 0.8027291
adjustedRandIndex(annotations_COCA$Tissue, annotations_COCA$TCGAcc) # 0.8059866
adjustedRandIndex(annotations_COCA$Tissue, annotations_COCA$OURcoca) # 0.794509
adjustedRandIndex(annotations_COCA$Tissue, annotations_COCA$OURcc) # 0.7258063

adjustedRandIndex(annotations_COCA$COCA, annotations_COCA$TCGAcoca) # 0.9466913 
adjustedRandIndex(annotations_COCA$COCA, annotations_COCA$TCGAcc) # 0.9513015
adjustedRandIndex(annotations_COCA$COCA, annotations_COCA$OURcoca) # 0.9263707
adjustedRandIndex(annotations_COCA$COCA, annotations_COCA$OURcc) # 0.8500842

##################################### Coincidences - OURcoca ####################################### 
bestK <- 10 

coincidences <- matrix(NA, 12, bestK)
tissue_name <- c("BLCA", "BRCA", "COAD", "GBM", "HNSC", "KIRC", "LAML", "LUAD", "LUSC", 
                 "OV", "READ",  "UCEC")
rownames(coincidences) <- tissue_name
colnames(coincidences) <- paste("cluster", 1:bestK)

TissueLabels <- annotations_COCA$Tissue
names(TissueLabels) <- rownames(annotations_COCA)

ClLabels <- annotations_COCA$OURcoca
names(ClLabels) <- rownames(annotations_COCA)

for(i in 1:12){
    tissue_i <- tissue_name[i]
    whichTissue_i <- names(TissueLabels)[which(TissueLabels == tissue_i)]
    for(j in 1:bestK){
        whichClusterK <- names(ClLabels)[which(ClLabels == j)]
        coincidences[i,j] <- sum(whichTissue_i %in% whichClusterK)
    }
}

col_fun = colorRamp2(c(0, max(coincidences)), c("white", "deeppink2"))

png("figures/coca-coincidences-OURcoca.png")
Heatmap(coincidences, 
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%d", as.integer(coincidences[i, j])), x, y, 
                      gp = gpar(fontsize = 10))
        },
        col = col_fun)
dev.off()

#################################### Coincidences - TCGAcoca ####################################### 

bestK <- 10 

coincidences <- matrix(NA, 12, bestK)
tissue_name <- c("BLCA", "BRCA", "COAD", "GBM", "HNSC", "KIRC", "LAML", "LUAD", "LUSC", 
                 "OV", "READ",  "UCEC")
rownames(coincidences) <- tissue_name
colnames(coincidences) <- paste("cluster", 1:bestK)

TissueLabels <- annotations_COCA$Tissue
names(TissueLabels) <- rownames(annotations_COCA)

ClLabels <- annotations_COCA$TCGAcoca
names(ClLabels) <- rownames(annotations_COCA)

for(i in 1:12){
    tissue_i <- tissue_name[i]
    whichTissue_i <- names(TissueLabels)[which(TissueLabels == tissue_i)]
    for(j in 1:bestK){
        whichClusterK <- names(ClLabels)[which(ClLabels == j)]
        coincidences[i,j] <- sum(whichTissue_i %in% whichClusterK)
    }
}

col_fun = colorRamp2(c(0, max(coincidences)), c("white", "deeppink2"))

png("figures/coca-coincidences-TCGAcoca.png")
Heatmap(coincidences, 
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%d", as.integer(coincidences[i, j])), x, y, 
                      gp = gpar(fontsize = 10))
        },
        col = col_fun)
dev.off()

###################################### Coincidences - OURcc ######################################## 

bestK <- 13

coincidences <- matrix(NA, 12, bestK)
tissue_name <- c("BLCA", "BRCA", "COAD", "GBM", "HNSC", "KIRC", "LAML", "LUAD", "LUSC", 
                 "OV", "READ",  "UCEC")
rownames(coincidences) <- tissue_name
colnames(coincidences) <- paste("cluster", 1:bestK)

TissueLabels <- annotations_COCA$Tissue
names(TissueLabels) <- rownames(annotations_COCA)

ClLabels <- annotations_COCA$OURcc
names(ClLabels) <- rownames(annotations_COCA)

for(i in 1:12){
    tissue_i <- tissue_name[i]
    whichTissue_i <- names(TissueLabels)[which(TissueLabels == tissue_i)]
    for(j in 1:bestK){
        whichClusterK <- names(ClLabels)[which(ClLabels == j)]
        coincidences[i,j] <- sum(whichTissue_i %in% whichClusterK)
    }
}

col_fun = colorRamp2(c(0, max(coincidences)), c("white", "deeppink2"))

png("figures/coca-coincidences-OURcc.png")
Heatmap(coincidences, 
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%d", as.integer(coincidences[i, j])), x, y, 
                      gp = gpar(fontsize = 10))
        },
        col = col_fun)
dev.off()

##################################### Coincidences - TCGAcc ######################################## 

bestK <- 13

coincidences <- matrix(NA, 12, bestK)
tissue_name <- c("BLCA", "BRCA", "COAD", "GBM", "HNSC", "KIRC", "LAML", "LUAD", "LUSC", 
                 "OV", "READ",  "UCEC")
rownames(coincidences) <- tissue_name
colnames(coincidences) <- paste("cluster", 1:bestK)

TissueLabels <- annotations_COCA$Tissue
names(TissueLabels) <- rownames(annotations_COCA)

ClLabels <- annotations_COCA$TCGAcc
names(ClLabels) <- rownames(annotations_COCA)

for(i in 1:12){
    tissue_i <- tissue_name[i]
    whichTissue_i <- names(TissueLabels)[which(TissueLabels == tissue_i)]
    for(j in 1:bestK){
        whichClusterK <- names(ClLabels)[which(ClLabels == j)]
        coincidences[i,j] <- sum(whichTissue_i %in% whichClusterK)
    }
}

col_fun = colorRamp2(c(0, max(coincidences)), c("white", "deeppink2"))

png("figures/coca-coincidences-TCGAcc.png")
Heatmap(coincidences, 
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%d", as.integer(coincidences[i, j])), x, y, 
                      gp = gpar(fontsize = 10))
        },
        col = col_fun)
dev.off()
