# load required functions
source("script2.R", local=FALSE)

path_to_working_directory <- getwd()
folder_with_plots <- "plots_Camargo_data_20220811"
dir.create(paste(path_to_working_directory, folder_with_plots, sep="/"), showWarnings = FALSE)


# load the count matrices of MPP subsets published by Camargo in 2018 (GSE90742) ####

# load MPP2
MPP2 <- read.csv(file="MPP_subsets_Camargo_GSE90742/MPP2.raw_umifm_counts.csv", sep=",")
unique_cell_names <- paste(MPP2[,2], MPP2[,3], MPP2[,1], sep = "_")
MPP2_expression <- t(MPP2[,6:ncol(MPP2)])
colnames(MPP2_expression) <- unique_cell_names
MPP2_meta_data <- MPP2[,1:5]
rownames(MPP2_meta_data) <- unique_cell_names
MPP2_seurat <- CreateSeuratObject(MPP2_expression, project = "SeuratProject", assay = "RNA",
                                  min.cells = 0, min.features = 0, names.field = 2,
                                  names.delim = "_", meta.data = MPP2_meta_data)
rm(MPP2_expression)
rm(MPP2_meta_data)
rm(MPP2)

# load MPP3
MPP3 <- read.csv(file="MPP_subsets_Camargo_GSE90742/MPP3.raw_umifm_counts.csv", sep=",")
unique_cell_names <- paste(MPP3[,2], MPP3[,3], MPP3[,1], sep = "_")
MPP3_expression <- t(MPP3[,6:ncol(MPP3)])
colnames(MPP3_expression) <- unique_cell_names
MPP3_meta_data <- MPP3[,1:5]
rownames(MPP3_meta_data) <- unique_cell_names
MPP3_seurat <- CreateSeuratObject(MPP3_expression, project = "SeuratProject", assay = "RNA",
                                  min.cells = 0, min.features = 0, names.field = 2,
                                  names.delim = "_", meta.data = MPP3_meta_data)
rm(MPP3_expression)
rm(MPP3_meta_data)
rm(MPP3)

# load MPP4
MPP4 <- read.csv(file="MPP_subsets_Camargo_GSE90742/MPP4.raw_umifm_counts.csv", sep=",")
unique_cell_names <- paste(MPP4[,2], MPP4[,3], MPP4[,1], sep = "_")
MPP4_expression <- t(MPP4[,6:ncol(MPP4)])
colnames(MPP4_expression) <- unique_cell_names
MPP4_meta_data <- MPP4[,1:5]
rownames(MPP4_meta_data) <- unique_cell_names
MPP4_seurat <- CreateSeuratObject(MPP4_expression, project = "SeuratProject", assay = "RNA",
                                  min.cells = 0, min.features = 0, names.field = 2,
                                  names.delim = "_", meta.data = MPP4_meta_data)
rm(MPP4_expression)
rm(MPP4_meta_data)
rm(MPP4)

# merge MPP subsets ####
MPP_subsets <- merge(MPP2_seurat, y = c(MPP3_seurat, MPP4_seurat), add.cell.ids = NULL, project = "U22")
rm(MPP2_seurat)
rm(MPP3_seurat)
rm(MPP4_seurat)

# remove all cells that did not pass quality filters in Camargo's publication ####
MPP_subsets <- MPP_subsets[,which(MPP_subsets$pass_filter==1)]

# annotate cells using SingleR and remove those that are annotated as mature cell types ####
MPP_subsets <- annotate_by_SingleR(object_to_annotate=MPP_subsets, path_helper_functions="annotation_helper_functions.R")
#saveRDS(object=MPP_subsets, file = "/Volumes/addition_storage/R_objects_Usp22/annotated_MPP_3_4_Camargo_20220609.rds")
#saveRDS(object=MPP_subsets, file = "/Volumes/addition_storage/R_objects_Usp22/annotated_MPP_2_3_4_Camargo_20220609.rds")
# also saved at /Users/metz/Nextcloud/backup_Usp22_data_and_Robjects/R_objects_Usp22/annotated_MPP_2_3_4_Camargo_20220609.rds
MPP_subsets <- MPP_subsets[,which(!(MPP_subsets$SingleR_labels %in% c("NK", "T", "T.DN", "preT", "proB", "B", "Baso", "Eo", "DC", "Mono", "EryBl", "InfMono", "Neut", "CFUE")))]

# normalize data ####
DefaultAssay(MPP_subsets) <- "RNA"
MPP_subsets <- Seurat::SCTransform(MPP_subsets, verbose = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)

# create diffusion map using the 45 marker genes used in our manuscript ####
common_features <- readRDS(file = "features_from_FateID_20220311")
common_features <- common_features[which(common_features %in% rownames(MPP_subsets))]

MPP_combined <- MPP_subsets
filtered_gene_cell_matrix <- t(as.matrix(MPP_combined@assays[["SCT"]]@data[common_features,]))
colnames(filtered_gene_cell_matrix) <- paste(rep("GENE", ncol(filtered_gene_cell_matrix)), as.character(c(1:ncol(filtered_gene_cell_matrix))), sep="_")
filtered_gene_cell_matrix_object <- CreateDimReducObject(
  embeddings = filtered_gene_cell_matrix,
  loadings = new(Class = "matrix"),
  projected = new(Class = "matrix"),
  assay = "SCT",
  stdev = numeric(),
  key = "GENE_",
  global = FALSE,
  jackstraw = NULL,
  misc = list()
)
MPP_combined@reductions[["filtered"]] <- filtered_gene_cell_matrix_object
MPP_combined.sce <- as.SingleCellExperiment(MPP_combined)

# remove cells that have the same values within the used marker genes
barcodes_to_remove <- MPP_combined.sce$barcode[which(duplicated(reducedDim(MPP_combined.sce,'FILTERED')))]
MPP_combined.sce <- MPP_combined.sce[, which(!(MPP_combined.sce$barcode %in% barcodes_to_remove))]
MPP_combined <- MPP_combined[, which(!(MPP_combined$barcode %in% barcodes_to_remove))]

# compute diffusion map based on selected features
dmap <- destiny::DiffusionMap(reducedDim(MPP_combined.sce,'FILTERED'))
reducedDim(MPP_combined.sce, 'DiffusionMap') <- dmap@eigenvectors

# visualize initial diffusion map
dmap_plot_unfiltered <- scatter_3d(MPP_combined.sce, dim.red.name = 'DiffusionMap', color.column="library_id", marker_size = 20, scene = '')
scatter_3d(MPP_combined.sce, dim.red.name = 'DiffusionMap', color.column="SingleR_labels", marker_size = 20, scene = '')

# assign each cell to a branch as done with the data of our manuscript ####
# create clusters
MPP_combined.seurat <- as.Seurat(MPP_combined.sce)
MPP_combined.seurat <- cluster_cells(MPP_combined.seurat, reduction = "DiffusionMap", ndims=3, knn=40, create_neighbor_object=FALSE, resolution=0.1)
MPP_combined.sce@colData@listData[["Cluster1"]] <- as.character(MPP_combined.seurat$seurat_clusters)
clusters_plot <- scatter_3d(MPP_combined.sce, dim.red.name = 'DiffusionMap', color.column="Cluster1", marker_size = 20, scene = '')

# select start and end clusters by cell type
SingleR_labels_marking_tip_clusters <- c("PreCFUE", "CLP", "GMP")
clusters_corresponding_to_SingleR_labels <- c()
for (index in (1:length(SingleR_labels_marking_tip_clusters))){
  current_label <- SingleR_labels_marking_tip_clusters[index]
  cluster_table <- table(MPP_combined.seurat$seurat_clusters[which(MPP_combined.seurat$SingleR_labels==current_label)])
  cluster_vector <- as.vector(cluster_table)
  names(cluster_vector) <- names(cluster_table)
  clusters_corresponding_to_SingleR_labels <- append(clusters_corresponding_to_SingleR_labels, names(cluster_vector[order(cluster_vector, decreasing = TRUE)])[1])
}
names(clusters_corresponding_to_SingleR_labels) <- SingleR_labels_marking_tip_clusters

# run slingshot
MPP_combined.sce@int_colData@listData[["reducedDims"]]@listData[["DiffusionMap_subset"]] <- MPP_combined.sce@int_colData@listData[["reducedDims"]]@listData[["DiffusionMap"]][,1:3]
MPP_combined.sce <- slingshot(MPP_combined.sce, clusterLabels = 'Cluster1', reducedDim = 'DiffusionMap_subset', start.clus=clusters_corresponding_to_SingleR_labels["PreCFUE"], end.clus=c(clusters_corresponding_to_SingleR_labels["CLP"], clusters_corresponding_to_SingleR_labels["GMP"]), allow.breaks=TRUE)

# assign curves to lineages and add pseudo-time to meta data
curve1_final_cluster <- MPP_combined.sce@int_metadata[["slingshot"]]@lineages[["Lineage1"]][length(MPP_combined.sce@int_metadata[["slingshot"]]@lineages[["Lineage1"]])]
curve1_final_label <- names(clusters_corresponding_to_SingleR_labels)[which(clusters_corresponding_to_SingleR_labels==curve1_final_cluster)]
pseudotime_values_with_NAs <- slingPseudotime(MPP_combined.sce, na = TRUE)
if (curve1_final_label=="CLP"){
  MPP_combined.sce@colData@listData[["GMP_pseudotime_na"]] <- pseudotime_values_with_NAs[,"curve2"]
  MPP_combined.sce@colData@listData[["CLP_pseudotime_na"]] <- pseudotime_values_with_NAs[,"curve1"]
} else {
  MPP_combined.sce@colData@listData[["GMP_pseudotime_na"]] <- pseudotime_values_with_NAs[,"curve1"]
  MPP_combined.sce@colData@listData[["CLP_pseudotime_na"]] <- pseudotime_values_with_NAs[,"curve2"]
}

# get pseudo-time of most central HSC
HSCs <- MPP_combined.sce[, which(MPP_combined.sce$SingleR_labels=="LTHSC")]
HSC_3d_coordinated <- HSCs@int_colData@listData[["reducedDims"]]@listData[["DiffusionMap_subset"]]
median_coordinate <- as.double(apply(HSC_3d_coordinated, 2, median))
distance_to_median_coordinate <- c()
for (index in (1:nrow(HSC_3d_coordinated))){
  x1 <- HSC_3d_coordinated[index, 1]
  x2 <- HSC_3d_coordinated[index, 2]
  x3 <- HSC_3d_coordinated[index, 3]
  distance_to_median_coordinate <- append(distance_to_median_coordinate, sqrt(((x1-median_coordinate[1])^2)+((x2-median_coordinate[2])^2)+((x3-median_coordinate[3])^2)))
}
central_HSC <- rownames(HSC_3d_coordinated)[which(distance_to_median_coordinate==min(distance_to_median_coordinate))]

# all cells with smaller pseudo-time than HSC get assigned to MEP trajectory (their CLP and GMP pseudo-times become NA)
indices_of_cells_in_MEP_trajectory <- which((MPP_combined.sce$CLP_pseudotime_na < MPP_combined.sce[,central_HSC]$CLP_pseudotime_na) | 
                                              (MPP_combined.sce$GMP_pseudotime_na < MPP_combined.sce[,central_HSC]$GMP_pseudotime_na))
MPP_combined.sce@colData@listData[["MEP_pseudotime_na"]] <- rep(NA, ncol(MPP_combined.sce))
MEP_pseudo_times <- c()
for (index in (1:length(indices_of_cells_in_MEP_trajectory))){
  CLP_time <- MPP_combined.sce$CLP_pseudotime_na[indices_of_cells_in_MEP_trajectory[index]]
  GMP_time <- MPP_combined.sce$GMP_pseudotime_na[indices_of_cells_in_MEP_trajectory[index]]
  MEP_pseudo_times <- append(MEP_pseudo_times, mean(c(CLP_time, GMP_time), na.rm = TRUE))
}
# pseudo-time values of MEP-trajectory are inverted to account for the directionality
MEP_pseudo_times_scaled <- (MEP_pseudo_times*(-1))+max(MEP_pseudo_times, na.rm=TRUE)
MPP_combined.sce$MEP_pseudotime_na[indices_of_cells_in_MEP_trajectory] <- MEP_pseudo_times_scaled
MPP_combined.sce$CLP_pseudotime_na[indices_of_cells_in_MEP_trajectory] <- NA
MPP_combined.sce$GMP_pseudotime_na[indices_of_cells_in_MEP_trajectory] <- NA

# cells are assigned to a branch if they are exclusively assigned to one slingshot trajectory
MPP_combined.sce@colData@listData[["branches"]] <- rep("undefined", ncol(MPP_combined.sce))
MPP_combined.sce$branches[which((!(is.na(MPP_combined.sce$CLP_pseudotime_na))) & 
                                  (is.na(MPP_combined.sce$GMP_pseudotime_na)) & 
                                  (is.na(MPP_combined.sce$MEP_pseudotime_na)))] <- "CLP"
MPP_combined.sce$branches[which((is.na(MPP_combined.sce$CLP_pseudotime_na)) & 
                                  (!(is.na(MPP_combined.sce$GMP_pseudotime_na))) & 
                                  (is.na(MPP_combined.sce$MEP_pseudotime_na)))] <- "GMP"
MPP_combined.sce$branches[which((is.na(MPP_combined.sce$CLP_pseudotime_na)) & 
                                  (is.na(MPP_combined.sce$GMP_pseudotime_na)) & 
                                  (!(is.na(MPP_combined.sce$MEP_pseudotime_na))))] <- "MEP"

MPP_combined@meta.data[["branches"]] <- MPP_combined.sce$branches


# plot diffusion map ####
if (file.exists("3d_scatter_parameters_Camargo.rds")==FALSE){
  plotcol <- MPP_combined.sce$branches
  plotcol[which(plotcol=="CLP")] <- "orange"
  plotcol[which(plotcol=="GMP")] <- "blue"
  plotcol[which(plotcol=="MEP")] <- "red"
  plotcol[which(plotcol=="undefined")] <- "grey"
  plot3d(reducedDims(MPP_combined.sce)$DiffusionMap_subset, col = plotcol, size = 4, alpha=0.5)
  answer <- readline("Please adjust the orientation of the plot with your mouse and confirm by writing done: ")
  if (answer=="done"){
    setting_3d_plot <- par3d(no.readonly=TRUE)
    saveRDS(object = setting_3d_plot, file = "3d_scatter_parameters_Camargo.rds")
  }
}

plotcol <- MPP_combined.sce$library_id
plotcol[which(plotcol=="MPP4")] <- "orange"
plotcol[which(plotcol=="MPP3")] <- "blue"
plotcol[which(plotcol=="MPP2")] <- "red"
plot3d(reducedDims(MPP_combined.sce)$DiffusionMap_subset, col = plotcol, size = 4, alpha=0.5)
setting_camargo <- readRDS(file="3d_scatter_parameters_Camargo.rds")
par3d(setting_camargo)
snapshot3d(filename = paste(folder_with_plots, "Camargo_subsets.png", sep="/"), 
           fmt = "png", width = 5000, height = 5000)
img <- image_read(paste(folder_with_plots, "Camargo_subsets.png", sep="/"))
image_write(img, path = paste(folder_with_plots, "Camargo_subsets.pdf", sep="/"), format = "pdf")
rgl.close()


plotcol <- MPP_combined.sce$branches
plotcol[which(plotcol=="CLP")] <- "orange"
plotcol[which(plotcol=="GMP")] <- "blue"
plotcol[which(plotcol=="MEP")] <- "red"
plotcol[which(plotcol=="undefined")] <- "grey"
plot3d(reducedDims(MPP_combined.sce)$DiffusionMap_subset, col = plotcol, size = 4, alpha=0.5)
setting_camargo <- readRDS(file="3d_scatter_parameters_Camargo.rds")
par3d(setting_camargo)
snapshot3d(filename = paste(folder_with_plots, "Camargo_subsets_branches.png", sep="/"), 
           fmt = "png", width = 5000, height = 5000)
img <- image_read(paste(folder_with_plots, "Camargo_subsets_branches.png", sep="/"))
image_write(img, path = paste(folder_with_plots, "Camargo_subsets_branches.pdf", sep="/"), format = "pdf")
rgl.close()
