plotlargerMOC = function(moc, datasetIndicator, datasetNames = NULL, annotations = NULL,
                   clr = FALSE, clc = FALSE, save = FALSE, fileName = "moc.png",
                   showObsNames = FALSE, showClusterNames = FALSE,
                   annotation_colors = NA){
    
    moc = t(moc)
    
    # Get number of datasets
    M = length(table(datasetNames))
    
    # Get sum of the number of clusters in each dataset
    sumK = dim(moc)[1]
    
    # If dataset names are not provided
    if(is.null(datasetNames)){
        
        # Dataset names = dataset indicators
        datasetNames = as.character(datasetIndicator)
    }
    
    M = length(table(datasetNames))
    
    # Make dataset names unique by adding a different number at the end of
    # the dataset name for each cluster
    datasetNamesLong <- c()
    for(i in 1:sumK){
        
        nClustersLeft <- sum(datasetNames[i:sumK] == datasetNames[i])
        datasetNamesLong <- c(datasetNamesLong, paste(datasetNames[i],
                                                      nClustersLeft, sep = " "))
    }
    
    # We want every dataset to have a different colour :)
    for(i in 1:sumK){
        moc[i,] <- moc[i,]*datasetIndicator[i]
    }
    
    # Plot!
    rownames(moc) <- datasetNamesLong
    # if(!is.null(annotations) && is.null(rownames(annotations))){
    #     rownames(annotations) <- colnames(moc)
    # }
    
    
    if(M==2){
        mycols <- c("white", (RColorBrewer::brewer.pal(n = 2, name = "RdBu")))
        mycols <- mycols[c(1,2,4)]
    }else{
        mycols <- c("white", (RColorBrewer::brewer.pal(n = max(3,M), name = "Set3")))
    }
    
    if(save) grDevices::png(fileName, width = 1000, height = 850)
    
    pheatmap::pheatmap(moc,  legend = TRUE,
                       legend_breaks = 0:M,
                       legend_labels = c("", unique(datasetNames)),
                       color =  mycols,
                       cluster_rows = clr, clustering_distance_rows = "binary",
                       cluster_cols = clc, clustering_distance_cols = "binary",
                       annotation_col = annotations, show_colnames = showObsNames,
                       annotation_colors = annotation_colors,
                       show_rownames = showClusterNames, drop_levels = FALSE, na_col = "seashell2")
    
    if(save){
        grDevices::dev.off()
        warning('After saving a pheatmap plot to file, you sometimes have to repeat the
                `dev.off()` command in order to shut down the plotting device completely.')
    }
}