################################ Pancancer study (Hoadley et al. 2014) #############################
######################################## Compare KLIC and COCA #####################################

rm(list=ls())

load("data/klic-clLabels-ALL.RData")
load("data/coca-clusters.RData")

klicK <- 10
cocaK <- 10 

coincidences <- matrix(NA, cocaK, klicK)

rownames(coincidences) <- paste("COCA cluster", 1:cocaK)
colnames(coincidences) <- paste("KLIC cluster", 1:klicK)

COCALabels <- annotations_COCA$OURcoca
names(COCALabels) <- rownames(annotations_COCA)

KLICLabels <- clLabels_ALL[klicK-1,]
names(KLICLabels) <- colnames(clLabels_ALL) 

for(i in 1:cocaK){
    whichCOCA_i <- names(COCALabels)[which(COCALabels == i)]
    cat(length(whichCOCA_i), '\n')
    for(j in 1:klicK){
        whichKLIC_j <- names(KLICLabels)[which(KLICLabels == j)]
        cat(length(whichKLIC_j), "\n")
        coincidences[i,j] <- sum(whichCOCA_i %in% whichKLIC_j)
    }
}

col_fun = colorRamp2(c(0, max(coincidences)), c("white", "deeppink2"))

png("figures/klic-coca-coincidences-10-10.png",
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
        title = "Coincidences",
        title_gp = gpar(fontsize = 28),
        labels_gp = gpar(fontsize = 28),
        direction = "vertical",
        nrow = 1,
        title_position = "lefttop-rot",
        legend_height = unit(10, "cm")
    )
)
dev.off()
