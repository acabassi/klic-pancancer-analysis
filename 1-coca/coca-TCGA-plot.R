############################
### Pancancer12 clusters ###
############################

rm(list = ls())
setwd("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-DATA/data-PANCAN12/")

library(pheatmap)
library(RColorBrewer)

COCAmatrix <- read.table("clusters/Table S1.tsv", header = TRUE, row.names = 1)
COCAmatrix <- COCAmatrix[1:3527,]

samples <- rownames(COCAmatrix)

nSCNA <- length(table(COCAmatrix$Copy_Number))
nMeth <- length(table(COCAmatrix$DNA_Methylation))
nRPPA <- length(table(COCAmatrix$RPPA))
nmRNA <- length(table(COCAmatrix$mRNA))
nmiRNA <- length(table(COCAmatrix$miRNA))
#nMutations <- length(table(COCAmatrix$Mutations))

n_clusters <- c(nSCNA, nMeth, nRPPA, nmRNA, nmiRNA) #, nMutations
max_n_clusters <- max(nSCNA, nMeth, nRPPA, nmRNA, nmiRNA) # , nMutations

N <- dim(COCAmatrix)[1]
n_datasets <- 5
dataset_name <- c("CN", "DNAmeth", "RPPA", "mRNA",  "miRNA")

labelMatrix <- array(NA, c(max_n_clusters, N, n_datasets))
for(i in 1:n_datasets){
  for(k in 1:n_clusters[i]){
    clusterLabel <- names(table(COCAmatrix[,i+2]))[k]
    labelMatrix[k,,i] <- (COCAmatrix[,i+2]==clusterLabel)*i
  }
}
#Convert label matrix from logic to numeric matrix
# labelMatrix <- labelMatrix*1

# Build MOC matrix
MOC <- rbind(labelMatrix[1:n_clusters[1],,1],
             labelMatrix[1:n_clusters[2],,2],
             labelMatrix[1:n_clusters[3],,3],
             labelMatrix[1:n_clusters[4],,4],
             labelMatrix[1:n_clusters[5],,5])

rownamesMOC <- rep(NA, sum(n_clusters[1:5]))
cont = 0
for(i in 1:n_datasets){
  for(k in 1:n_clusters[i]){
    cont = cont + 1
    rownamesMOC[cont] <- paste(dataset_name[i], as.character(k), sep="")
  }
}

dim(MOC)[1]
length(rownamesMOC)
rownames(MOC) <- rownamesMOC

colnames(MOC) <- rownames(COCAmatrix)

anno_col <- as.data.frame(COCAmatrix[,c(2,1)])
table(anno_col$COCA)
cluster_names <- c("1 - LUAD-enriched",
                   "2 - Squamous-like",
                   "3 - BRCA/Luminal",
                   "4 - BRCA/Basal",
                   "5 - KIRC",
                   "6 - UCEC",
                   "7 - COAD/READ",
                   "8 - BLCA",
                   "9 - OV",
                   "10 - GBM",
                   "11 - small-various",
                   "12 - small-various",
                   "13 - AML")
for(i in 1:13){
  anno_col$COCA[which(anno_col$COCA == i)] <- cluster_names[i]
}
anno_col$COCA <- as.factor(anno_col$COCA)
anno_col$COCA <- factor(anno_col$COCA,cluster_names)
typeof(anno_col$COCA)

# MOC <-factor(MOC, c(0,1,2,3,4,5))
pdf("COCA_Hoadley2014.pdf", width = 9, height = 10)
pheatmap(MOC, kmeans_k = NA, cluster_rows = TRUE, clustering_distance_rows = "binary",
         cluster_cols = FALSE, 
         show_colnames = FALSE,  annotation_col = anno_col,
         legend_breaks = c(0,1,2,3,4,5), legend_labels = c("0", dataset_name),
         color =  c("white",(brewer.pal(n = 5, name = "Set2"))))
dev.off()
#, 
#         legend = FALSE,

colSums(is.na(MOC))
plot(colSums(is.na(MOC)))

save(samples, anno_col, file = "samples.RData")
      