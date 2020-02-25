################################ Pancancer study (Hoadley et al. 2014) #############################
################################################# KLIC #############################################

rm(list = ls())

# library(devtools)
# install_github("acabassi/klic")
library(klic)
# install_github("acabassi/coca")
library(coca)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

library(circlize)
library(colorspace)

############################################# Load kernels #########################################

n_datasets <- 5

datasetNames <- c("CN", "Meth", "miRNA", "mRNA", "RPPA")

# load("data/kernel-cn.RData")
# load("data/names-CN.RData")
# rownames(cc) <- colnames(cc) <- which_ones
# ccCN <- cc
# 
# load("data/kernel-methylation.RData")
# load("data/names-methylation.RData")
# rownames(cc) <- colnames(cc) <- which_ones
# ccMeth <- cc
# 
# load("data/kernel-miRNA.RData")
# load("data/names-miRNA.RData")
# rownames(cc) <- colnames(cc) <- which_ones
# ccmiRNA <- cc
# 
# load("data/kernel-mRNA.RData")
# load("data/names-mRNA.RData")
# rownames(cc) <- colnames(cc) <- which_ones
# ccmRNA <- cc
# 
# load("data/kernel-RPPA.RData")
# load("data/names-RPPA.RData")
# rownames(cc) <- colnames(cc) <- which_ones
# ccRPPA <- cc
# rm(cc)

############################################# Spectrum shift #######################################

# ccCN    <- spectrumShift(ccCN,    verbose = TRUE)
# ccMeth  <- spectrumShift(ccMeth,  verbose = TRUE)
# ccmiRNA <- spectrumShift(ccmiRNA, verbose = TRUE)
# ccmRNA  <- spectrumShift(ccmRNA,  verbose = TRUE)
# ccRPPA  <- spectrumShift(ccRPPA,  verbose = TRUE)
# save(ccCN, ccMeth, ccmiRNA, ccmRNA, ccRPPA, file = "data/shiftedMatrices.RData")

load("data/shiftedMatrices.RData")

###################################### Create array of all matrices ################################

load("data/samples.RData")

# allSamples <- unique(c(rownames(ccCN), rownames(ccMeth), rownames(ccmiRNA), rownames(ccmRNA),
#                        rownames(ccRPPA)))
N <- length(samples)
CM <- array(NA, c(N, N, n_datasets))
dimnames(CM) <- list(samples, samples, datasetNames)

CM[rownames(ccCN),    colnames(ccCN),    1] <- ccCN
CM[rownames(ccMeth),  colnames(ccMeth),  2] <- ccMeth
CM[rownames(ccmiRNA), colnames(ccmiRNA), 3] <- ccmiRNA
CM[rownames(ccmRNA),  colnames(ccmRNA),  4] <- ccmRNA
CM[rownames(ccRPPA),  colnames(ccRPPA),  5] <- ccRPPA

###################################### Indicator of missing values #################################

missing <- matrix(TRUE, N, n_datasets)
dimnames(missing) <- list(samples, datasetNames)
missing[rownames(ccCN),    "CN"] <- FALSE
missing[rownames(ccMeth),  "Meth"] <- FALSE
missing[rownames(ccmiRNA), "miRNA"] <- FALSE
missing[rownames(ccmRNA),  "mRNA"] <- FALSE
missing[rownames(ccRPPA),  "RPPA"] <- FALSE

################################################# KLIC #############################################

# Use localised multiple kernel k-means to integrate the datasets

parameters <- list()
# Set the maximum number of iterations for localised multiple kernel k-means
parameters$iteration_count <- 250

# Maximum number of clusters considered
maxK <- 20

# Each iteration was done separately on the high-performance computing (HPC) cluster ---------------
# i <- 13 # for(i in 2:maxK){
# 
# # Use kernel k-means with K=i to find weights and cluster labels
# parameters$cluster_count <- i # set the number of clusters K
# lmk <- lmkkmeans_missingData(CM, parameters, missing, verbose = TRUE)
# 
# # Save cluster labels
# clLabels <- lmk$clustering
# 
# # Save weights
# clWeights <- lmk$Theta
# 
# # } # end of for loop
# --------------------------------------------------------------------------------------------------

# save(
#     KM,
#     clLabels,
#     clWeights,
#     file = paste(
#         "data/klic-output",
#         i,
#         "_",
#         parameters$iteration_count,
#         "iterations.RData",
#         sep = ""
#     )
# )

########################## Put together kernel matrices and cluster labels #########################

KM_ALL <- array(NA, c(N, N, maxK - 1))
clLabels_ALL <- array(NA, c(maxK - 1, N))
clWeights_ALL <- array(NA, c(N, n_datasets, maxK - 1))

for (i in 2:maxK) {
   cat("Number of clusters: ", i, "\n")

   # Load output of KLIC, computed using the HPC cluster
   load(
      paste(
         "data/klic-output",
         i,
         "_",
         parameters$iteration_count,
         "iterations.RData",
         sep = ""
      )
   )
   rm(KM)
   
   # Compute weighted matrix
   KM  <- matrix(0, N, N)
   for (j in 1:dim(CM)[3]) {
      avail_j <- 1 - missing[, j]
      index_avail_j <- which(avail_j==1)
      KM[index_avail_j, index_avail_j] <-
         KM[index_avail_j, index_avail_j] + (clWeights[index_avail_j, j] %*% t(clWeights[index_avail_j, j])) * CM[index_avail_j, index_avail_j, j]
   }

   KM_ALL[, , i - 1] <- KM
   clLabels_ALL[i - 1,] <- clLabels
   clWeights_ALL[, , i - 1] <- clWeights
}

dimnames(KM_ALL) <-
    list(rownames(CM), colnames(CM), as.character(2:maxK))
dimnames(clLabels_ALL) <- list(as.character(2:maxK), rownames(CM))

save(KM_ALL, clLabels_ALL, file = paste("data/klic-maxK", maxK, ".RData", sep = ""))
save(clLabels_ALL, file = paste("data/klic-clLabels-maxK", maxK, ".RData", sep = ""))
save(clWeights_ALL, file = paste("data/klic-theta-", maxK, ".RData", sep = ""))

# Check output
randomK <- sample(1:(maxK-1), 1)
observations <- sample(1:(dim(KM_ALL)[1]), 1000, replace=FALSE)
HKM <-
   Heatmap(KM_ALL[observations, observations, randomK],
           show_row_names = FALSE,
           show_column_names = FALSE)
Hclusters <-
   Heatmap(clLabels_ALL[randomK, observations],
           show_row_names = FALSE,
           show_column_names = FALSE)

HKM + Hclusters

######################################### Maximise silhouette ######################################

maxSil <- maximiseSilhouette(KM_ALL, clLabels_ALL, maxK = maxK)
bestK <- maxSil$K


png(
  "figures/klic_silhouette.png",
  width = 500,
  height = 500
)
plot(
  2:20,
  maxSil$silhouette,
  type = "b",
  xlab = "Clusters",
  ylab = "",
  cex = 1.2,
  cex.lab = 1.5,
  cex.axis = 1.2
)
dev.off()

#################################### Plot weighted kernel matrix ###################################

# Let's consider 10 clusters instead of 7, in order to compare to COCA

n_clusters <- 10

my_blues <-
  colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)

HKM <- Heatmap(
  KM_ALL[, , n_clusters - 1]/max(KM_ALL[,,n_clusters-1]),
  col = my_blues,
  show_row_names = FALSE,
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title = "Weighted kernel",
    title_gp = gpar(fontsize = 28),
    labels_gp = gpar(fontsize = 28),
    direction = "vertical",
    nrow = 1,
    title_position = "lefttop-rot",
    legend_height = unit(10, "cm")
  )
)

clusters_palette <- qualitative_hcl(n_clusters, palette = "Dynamic")
names(clusters_palette) <- as.character(1:10)
annotations <- data.frame(Clusters = as.factor(clLabels_ALL[n_clusters - 1, ]))
rownames(annotations) <- colnames(clLabels_ALL[])

Hclusters <-
  HeatmapAnnotation(
    Clusters = annotations$Clusters,
    col = list(Clusters = clusters_palette),
    which = "row",
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      title = "Clusters",
      labels_gp = gpar(fontsize = 28),
      title_gp = gpar(fontsize = 28),
      # nrow = 1,
      direction = "vertical",
      title_position = "lefttop-rot",
      grid_height = unit(1, "cm")
    )
    
  )

# png("figures/weightedKernel_bestK.png", height = 800, width = 850)
# draw(
#   HKM + Hclusters,
#   merge_legend = TRUE,
#   heatmap_legend_side = "right",
#   annotation_legend_side = "right"
# )
# dev.off()

############################################ Plot weighs ####################################¢######

weights_palette <- rev(sequential_hcl(10, palette = "BluGrn"))
dimnames(clWeights_ALL) <- list(rownames(missing), 
                                c("DNA copy number",
                                  "DNA methylation",
                                  "mRNA expression",
                                  "mRNA expression",
                                  "Protein expression"), as.character(2:20))

annotations_with_names <- annotations$Clusters
names(annotations_with_names) <- rownames(annotations)

Hclusters <-
  HeatmapAnnotation(
    Clusters = annotations_with_names,
    col = list(Clusters = clusters_palette),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      title = "Clusters",
      labels_gp = gpar(fontsize = 22),
      title_gp = gpar(fontsize = 28),
      nrow = 1,
      direction = "horizontal"
      # title_position = "lefttop-rot",
      # grid_width = unit(, "cm")
    )
    
  )

multiplier <- matrix(NA, dim(missing)[1], dim(missing)[2])
multiplier[which(missing==FALSE)] <- 1

Hweights <-
  Heatmap(
    t(clWeights_ALL[, ,n_clusters-1])*t(multiplier),
    col = weights_palette,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    name = "Weights",
    bottom_annotation = Hclusters,
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 22),
      title_gp = gpar(fontsize = 28),
      direction = "horizontal",
      nrow = 1,
      legend_height = unit(4.5, "cm")
      # title_position = "topcenter"
    ),
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 22),
    # column_names_centered = TRUE,
    width = unit(20, "cm"),
    height = unit(5, "cm")
  )

png("figures/klic-weights.png", height = 250, width = 800)
draw(Hweights, merge_legends = TRUE, heatmap_legend_side = "bottom")
dev.off()


################################### Coincidence matrix / Best K ####################################

coincidences <- matrix(NA, 12, bestK)

tissue_name <-
    c(
        "BLCA",
        "BRCA",
        "COAD",
        "GBM",
        "HNSC",
        "KIRC",
        "LAML",
        "LUAD",
        "LUSC",
        "OV",
        "READ",
        "UCEC"
    )
rownames(coincidences) <- tissue_name
colnames(coincidences) <- paste("KLIC cluster", 1:bestK)

TissueLabels <- anno_col$Tissue
names(TissueLabels) <- rownames(anno_col)

ClLabels <- clLabels_ALL[bestK - 1, ]
names(ClLabels) <- rownames(CM)

for (i in 1:12) {
    tissue_i <- tissue_name[i]
    whichTissue_i <-
        names(TissueLabels)[which(TissueLabels == tissue_i)]
    for (j in 1:bestK) {
        whichClusterK <- names(ClLabels)[which(ClLabels == j)]
        coincidences[i, j] <- sum(whichTissue_i %in% whichClusterK)
    }
}

col_fun = colorRamp2(c(0, max(coincidences)), c("white", "deeppink2"))

png("figures/klic-coincidences-bestK.png",
    height = 800, 
    width = 900)
Heatmap(
    coincidences,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%d", as.integer(coincidences[i, j])), x, y,
                  gp = gpar(fontsize = 28))
    },
    col = col_fun,
    row_names_gp = gpar(fontsize = 22),
    column_names_gp = gpar(fontsize = 22),
    heatmap_legend_param = list(
      title = "Weighted kernel",
      title_gp = gpar(fontsize = 28),
      labels_gp = gpar(fontsize = 28),
      direction = "vertical",
      nrow = 1,
      title_position = "lefttop-rot",
      legend_height = unit(10, "cm")
    )
)
dev.off()

# plotSimilarityMatrix2(KM_ALL[,,bestK-1], clusLabels = clLabels_ALL[bestK-1,])

################################### Coincidence matrix / K = 10 ####################################

bestK <- 10

coincidences <- matrix(NA, 12, bestK)
tissue_name <-
    c(
        "BLCA",
        "BRCA",
        "COAD",
        "GBM",
        "HNSC",
        "KIRC",
        "LAML",
        "LUAD",
        "LUSC",
        "OV",
        "READ",
        "UCEC"
    )
rownames(coincidences) <- tissue_name
colnames(coincidences) <- paste("KLIC cluster", 1:bestK)

TissueLabels <- anno_col$Tissue
names(TissueLabels) <- rownames(anno_col)

ClLabels <- clLabels_ALL[bestK - 1, ]
names(ClLabels) <- rownames(CM)

for (i in 1:12) {
    tissue_i <- tissue_name[i]
    whichTissue_i <-
        names(TissueLabels)[which(TissueLabels == tissue_i)]
    for (j in 1:bestK) {
        whichClusterK <- names(ClLabels)[which(ClLabels == j)]
        coincidences[i, j] <- sum(whichTissue_i %in% whichClusterK)
    }
}

png("figures/klic-coincidences-10Cl.png",
    height = 800,
    width = 900)
Heatmap(
  coincidences,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%d", as.integer(coincidences[i, j])), x, y,
              gp = gpar(fontsize = 28))
  },
  col = col_fun,
  row_names_gp = gpar(fontsize = 22),
  column_names_gp = gpar(fontsize = 22),
  heatmap_legend_param = list(
    title = "Weighted kernel",
    title_gp = gpar(fontsize = 28),
    labels_gp = gpar(fontsize = 28),
    direction = "vertical",
    nrow = 1,
    title_position = "lefttop-rot",
    legend_height = unit(10, "cm")
  )
)
dev.off()

################################### Coincidence matrix / K = 13 ####################################

bestK <- 13

coincidences <- matrix(NA, 12, bestK)
tissue_name <-
    c(
        "BLCA",
        "BRCA",
        "COAD",
        "GBM",
        "HNSC",
        "KIRC",
        "LAML",
        "LUAD",
        "LUSC",
        "OV",
        "READ",
        "UCEC"
    )
rownames(coincidences) <- tissue_name
colnames(coincidences) <- paste("cluster", 1:bestK)

TissueLabels <- anno_col$Tissue
names(TissueLabels) <- rownames(anno_col)

ClLabels <- clLabels_ALL[bestK - 1, ]
names(ClLabels) <- rownames(CM)

for (i in 1:12) {
    tissue_i <- tissue_name[i]
    whichTissue_i <-
        names(TissueLabels)[which(TissueLabels == tissue_i)]
    for (j in 1:bestK) {
        whichClusterK <- names(ClLabels)[which(ClLabels == j)]
        coincidences[i, j] <- sum(whichTissue_i %in% whichClusterK)
    }
}

png("figures/klic-coincidences-13Cl.png",
    height = 800, 
    width = 900)
Heatmap(
  coincidences,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%d", as.integer(coincidences[i, j])), x, y,
              gp = gpar(fontsize = 28))
  },
  col = col_fun,
  row_names_gp = gpar(fontsize = 22),
  column_names_gp = gpar(fontsize = 22),
  heatmap_legend_param = list(
    title = "Weighted kernel",
    title_gp = gpar(fontsize = 28),
    labels_gp = gpar(fontsize = 28),
    direction = "vertical",
    nrow = 1,
    title_position = "lefttop-rot",
    legend_height = unit(10, "cm")
  )
)
dev.off()
