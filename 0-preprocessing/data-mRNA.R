############################### Pancancer study (Hoadley et al. 2014) ##############################
###################################### mRNA data preprocessing #####################################

rm(list = ls())

library(data.table)
library(mclust)
library(pheatmap)
library(pracma)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)

## From Hoadley et al. (2014)

# Using the platform corrected mRNAseq data, genes were filtered for those present in 70% of samples
# and then the top 6,000 most variable genes were selected. ConsensusClusterPlus R-package was used
# to identify clusters in the data using 1000 iterations, 80% sample resampling from 2 to 20
# clusters (k2 to k20) using hierarchical clustering with average innerLinkage and finalLinkage
# and Pearson correlation as the similarity metric. Eleven main groups were identified when 16
# clusters were used (Figure S1A). These 11 groups were observed to be stable through the use of 20
# clusters (K20) and significant in pairwise comparisons of the 11 main clusters with SigClust [11].

############################################# Load data ############################################

# The data can be downloaded from: https://www.synapse.org/Portal.html#!Synapse:syn1715755
# The clusters can be downloaded from: https://www.synapse.org/Portal.html#!Synapse:syn1715788

mRNAseq <-
    fread(
        "data-TCGA/PanCan12.3602-corrected-v3.txt",
        header = TRUE,
        sep = "\t",
        verbose = TRUE
    )
mRNAseqSaved <- mRNAseq
dim(mRNAseq) # 16117 variables X 3602 samples (1st row is the sample name, so 16116)
colnamesMRNA <- mRNAseq[1, 2:ncol(mRNAseq)]
rownamesMRNA <- as.vector(t(mRNAseqSaved[2:nrow(mRNAseqSaved), 1]))
mRNAseq <- mRNAseq[, 2:ncol(mRNAseq)]

mRNAseq <- t(mRNAseq) # This should be a 3602 x 16117 matrix

# The first column contains the sample labels
rownames(mRNAseq) <- substring(mRNAseq[, 1], 1, 12)
mRNAseq <- mRNAseq[, -1]
colnames(mRNAseq)  <- rownamesMRNA
class(mRNAseq) <- "numeric"
#save(mRNAseq,file = "data/mRNAseq.RDa")

########################################### Load clusters ##########################################

mRNAseqSubtypes <-
    fread(
        "data-TCGA/PanCan12.3602-corrected-v3.Subtypes.K16.txt",
        header = TRUE,
        sep = "\t"
    )
mRNAseqSubtypesLongNames <- mRNAseqSubtypes$Sample
mRNAseqSubtypesShortNames <-
    sapply(mRNAseqSubtypesLongNames, substr, 1, 12)

#Check that we have the same sample names:
all(mRNAseqSubtypesShortNames == rownames(mRNAseq))
# [1] TRUE  # Yes, these are the right samples!  And they are in the same order as in mRNAseqSubtypes

########################################### Select genes ###########################################

# Take only genes that are present in at least 70% of the samples
hist(colSums(is.na(mRNAseq)))
abline(v = dim(mRNAseq)[1] * 0.7, col = "red")

mRNAseqNumericSaved <- mRNAseq
mRNAseq <- mRNAseqNumericSaved
mRNAseq <-
    mRNAseqNumericSaved[,
                        -which(colSums(is.na(mRNAseqNumericSaved)) > 
                                   (dim(mRNAseqNumericSaved)[1] * 0.3))]

# Calculate variances and take 6k genes with highest variance
var_nonzero <- function(x)
    var(x[!is.na(x)])
variances <- apply(mRNAseq, 2, var_nonzero)
hist(variances)
threshold <- tail(sort(variances), 6000)[1]
mRNAseq <- mRNAseq[, which(variances >= threshold)]

# save(mRNAseq, file = "data/preprocessed-mRNA-data.RData")

############################################## Cluster #############################################

# The below should reproduce clusters obtained in paper (up to stochastic variability in sampling)

# CAUTION: This takes a long time to run!
# ccOut <-
#     ConsensusClusterPlus(
#         mRNAseq,
#         maxK = 20,
#         innerLinkage = "average",
#         finalLinkage = "average",
#         distance = "pearson",
#         corUse = "pairwise.complete.obs",
#         reps = 10
#     )

# Check if the clusters obtained above when K = 16 match (approximately)
# the clustering defined by mRNAseqSubtypes$K16.  What is the ARI between the clusterings?
# How do the results of the above change as reps is increased?  In the paper they said they used 1000.
# It is the consensus matrix obtained as a result of using the above procedure with K = 16
# (i.e. ccOut[[16]]$consensusMatrix) that we should use as the input to KLIC.

load('data/ccOut.RData')

obs <- names(ccOut[[16]]$consensusClass)

# Remove duplicated observations
obs <- obs[-which(duplicated(obs))]

### Load COCA clusters ----------------------------------------------------
load("data/samples.RData")

COCAmatrix <- read.table("data-TCGA/Table S1.tsv", header = TRUE, row.names = 1)
COCAmatrix <- COCAmatrix[1:3527,] # Remove rows containing only NAs

anno_col <- cbind(anno_col, COCAmatrix[rownames(anno_col),-c(1,2)])
sample_labels <- rownames(anno_col)
### -----------------------------------------------------------------------

########################################### Plot clusters ##########################################

annotations <-
    data.frame(H = as.factor(COCAmatrix[obs, ]$mRNA),
               C = as.factor(ccOut[[16]]$consensusClass[obs]))
rownames(annotations) <- obs

annotations <- annotations[order(annotations$C),]



annotations <- annotations[-which(is.na(annotations$H)),]

adjustedRandIndex(annotations$C, annotations$H)
# 0.9167216

obs <- rownames(annotations)

save(annotations, file = "data/clusters-mRNAseq.RData")

negative_numbers <- linspace(min(mRNAseq, na.rm = TRUE), 0, n = 16)
positive_numbers <-
    linspace(0, max(mRNAseq, na.rm = TRUE), n = ceil((
        max(mRNAseq, na.rm = TRUE) / (negative_numbers[2] - negative_numbers[1])
    )))
col_breaks <-
    c(negative_numbers, positive_numbers[2:length(positive_numbers)])
my_colours <-
    colorRampPalette(c("#FF9900", "white"))(length(negative_numbers))
my_colours <-
    c(my_colours, colorRampPalette(c("white", "#146EB4"))(length(positive_numbers)))

some_columns <- sample(1:(dim(mRNAseq)[2]), 100, replace = FALSE)

png("figures/mRNAseq-hc.png",
    height = 650,
    width = 100 + 100 * 10
)
pheatmap(
    mRNAseq[rownames(annotations), some_columns],
    color = my_colours,
    breaks = col_breaks,
    cluster_rows = F,
    annotation_row = annotations,
    show_rownames = F,
    show_colnames = F,
    fontsize = 20,
    cellwidth = 8
)
dev.off()
