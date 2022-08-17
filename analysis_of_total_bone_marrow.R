# analysis of total bone marrow ####
if (load_intermediate_results == TRUE){
  TBM_combined <- readRDS(file = paths_of_objects$TBM_combined_path)
} else {
  set.seed(1)
  # load and process total bone marrow of WT
  TBM_WT <- load_process_data(path=path_TBM_WT, nFeature_lower_threshold=700, nFeature_upper_threshold=10000, percent.mt_threshold=7, including_normalization=TRUE)
  TBM_WT <- cluster_cells(object=TBM_WT, reduction = "pca", ndims=30, knn=30, create_neighbor_object=FALSE, resolution=0.8)
  TBM_WT <- remove_doublets(object=TBM_WT, last_PC=30, ndims=30)
  
  # load and process total bone marrow of KO
  TBM_KO <- load_process_data(path=path_TBM_KO, nFeature_lower_threshold=700, nFeature_upper_threshold=10000, percent.mt_threshold=7, including_normalization=TRUE)
  TBM_KO <- cluster_cells(object=TBM_KO, reduction = "pca", ndims=30, knn=30, create_neighbor_object=FALSE, resolution=0.8)
  TBM_KO <- remove_doublets(object=TBM_KO, last_PC=30, ndims=30)
  
  # combine WT and KO cells into one object
  TBM_combined <- merge_conditions(object_WT=TBM_WT, object_KO=TBM_KO)
  
  # cluster cells
  TBM_combined <- cluster_cells(object=TBM_combined, reduction = "pca", ndims=30, knn=30, create_neighbor_object=FALSE, resolution=0.8)
  
  # compute umap
  TBM_combined <- Seurat::RunUMAP(TBM_combined, dims = 1:30, return.model = TRUE, n.neighbors = 28, min.dist = 0.45, metric="cosine")#30, 0.3 # 28, 0.45
  DimPlot(TBM_combined, reduction = "umap", label = TRUE)+ NoLegend()
  
  # annotate cells by SingleR
  TBM_combined <- annotate_by_SingleR(object_to_annotate=TBM_combined, path_helper_functions="annotation_helper_functions.R")
  
  # adjust annotations to cluster boundaries and pool some populations
  TBM_combined <- get_coarse_annotation(object_to_annotate=TBM_combined)
  
  if (save_intermediate_results==TRUE){
    saveRDS(object = TBM_combined, file = paths_of_objects$TBM_combined_path)
  }
}

# create Fig1C
plot_Fig1C(TBM_combined, path_for_plot=folder_with_plots)

# compute log2(fold-changes) of cell type frequencies between conditions with 95% confidence intervalls; write pdf
quantify_compositional_changes_of_cell_types(TBM_combined, path_for_plot=paste(folder_with_plots, "changes_of_cell_type_frequencies.pdf", sep="/"), pseudocount=0.00)













# clean work space
elements_to_remove <- setdiff(ls(), lsf.str())
elements_to_remove <- c(elements_to_remove, "elements_to_remove")
elements_to_remove <- elements_to_remove[which(!(elements_to_remove %in% c("path_HSC_KO", "path_HSC_WT", "path_Linminus_KO", 
                                                                           "path_Linminus_WT", "path_LSK_KO", "path_LSK_WT",
                                                                           "path_MPP_KO", "path_MPP_WT", "path_ST_KO", 
                                                                           "path_ST_WT", "path_TBM_KO", "path_TBM_WT",
                                                                           "folder_with_plots", "replicate", "save_intermediate_results",
                                                                           "load_intermediate_results", "paths_of_objects", "replicate_parameters",
                                                                           "s.genes", "g2m.genes")))]
rm(list = elements_to_remove)
gc()
