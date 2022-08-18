# function for loading data and doing first processing steps
load_process_data <- function(path, nFeature_lower_threshold, nFeature_upper_threshold, percent.mt_threshold, including_normalization){
  object.data <- Seurat::Read10X(data.dir = path)
  object <- Seurat::CreateSeuratObject(counts = object.data, project = "U22", min.cells = 3, min.features = 200)
  object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, pattern = "^mt-")
  object[["ratio_nFeature_counts"]] <- object[["nFeature_RNA"]]/object[["nCount_RNA"]]
  object <- subset(object, subset = nFeature_RNA > nFeature_lower_threshold & nFeature_RNA < nFeature_upper_threshold & percent.mt < percent.mt_threshold)
  DefaultAssay(object = object) <- "RNA"
  if (including_normalization==TRUE){
    object <- Seurat::SCTransform(object, verbose = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)
    DefaultAssay(object = object) <- "SCT"
    object <- Seurat::RunPCA(object, features = Seurat::VariableFeatures(object = object))
  }
  rm(object.data)
  return(object)
}
# function for clustering
cluster_cells <- function(object, reduction = "pca", ndims, knn, create_neighbor_object, resolution){
  object <- Seurat::FindNeighbors(object, dims = 1:ndims, k.param = knn, return.neighbor=create_neighbor_object, reduction = reduction)
  object <- Seurat::FindClusters(object, resolution = resolution)
  return(object)
}
# function for doublet detection and removal
remove_doublets <- function(object, last_PC, ndims){
  sce_object <- as.SingleCellExperiment(object)
  sce_object@assays@data@listData[["logcounts"]] <- NULL
  dbl.out <- scDblFinder::scDblFinder(sce_object, clusters=sce_object$seurat_clusters, dims = ndims, includePCs=1:last_PC)
  object@meta.data[["doublets"]] <- dbl.out$scDblFinder.class
  object <- subset(object, subset = doublets == "singlet")
  rm(dbl.out)
  return(object)
}
# function to merge corresponding samples of two conditions
merge_conditions <- function(object_WT, object_KO){
  DefaultAssay(object_WT) <- "RNA"
  DefaultAssay(object_KO) <- "RNA"
  object_combined <- merge(object_WT, y = c(object_KO), add.cell.ids = c("WT", "KO"), project = "Usp22", merge.data=FALSE)
  identity_vector <- c(rep("WT", ncol(object_WT)), rep("KO", ncol(object_KO)))
  object_combined@meta.data[["identity"]] <- identity_vector
  object_combined <- Seurat::SCTransform(object_combined, verbose = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)
  DefaultAssay(object = object_combined) <- "SCT"
  object_combined <- Seurat::RunPCA(object_combined, features = Seurat::VariableFeatures(object = object_combined))
  return(object_combined)
}
# function for mapping query sample onto reference in reduced dimension
map_sample <- function(ref_object, query_object, ndims, reference_reduction="pca", normalization_method="SCT", label_to_transfer="seurat_clusters", reduction_model="umap"){
  anchors <- FindTransferAnchors(
    reference = ref_object,
    query = query_object,
    normalization.method = normalization_method,
    reference.reduction = reference_reduction,
    dims = 1:ndims
  )
  
  if (is.null(label_to_transfer)){
    data_to_transfer <- NULL
  } else{
    data_to_transfer <- list(celltype = label_to_transfer)
  }
  
  query_object <- MapQuery(
    anchorset = anchors,
    query = query_object,
    reference = ref_object,
    refdata = data_to_transfer,
    reference.reduction = reference_reduction,
    reduction.model = reduction_model
  )
  return(query_object)
}
# function for annotating cells using SingleR
annotate_by_SingleR <- function(object_to_annotate, path_helper_functions="annotation_helper_functions.R"){
  assay_to_use <- "RNA"
  source(path_helper_functions, local = TRUE)
  # 1.1. read reference data as raw counts
  ref_list <- .read_ref_raw_counts()
  ref_data_immgen <- ref_list[["ref_data_immgen"]]
  ref_data_haemoshere <- ref_list[["ref_data_haemoshere"]]
  haemoshere_sample_table <- ref_list[["haemoshere_sample_table"]]
  ref_data_immgen2 <- ref_list[["ref_data_immgen2"]]
  
  # translate ensemble gene names of haemosphere data
  rownames(ref_data_haemoshere) <- .translate_gene_ID(gene_ID_vector=rownames(ref_data_haemoshere))
  
  # 1.2. extract only samples of interest (cell types that were purified from total bone marrow are used if available; otherwise samples from other tissues are used)
  ref_data_haemoshere <- .clean_haemosphere_data(haemoshere_sample_table=haemoshere_sample_table,
                                                 ref_data_haemoshere=ref_data_haemoshere)
  ref_data_immgen_multiple_tissues <- .extract_T_NKT_cells_from_immgen_which_are_not_from_bone_marrow(ref_data_immgen=ref_data_immgen)
  ref_data_immgen <- ref_data_immgen[,colnames(ref_data_immgen)[which(grepl("BM", colnames(ref_data_immgen), fixed = TRUE))]]
  ref_data_immgen2 <- ref_data_immgen2[,colnames(ref_data_immgen2)[which(grepl("BM", colnames(ref_data_immgen2), fixed = TRUE))]]
  
  # 2. select genes common in all reference and test sets and use those to subset the reference data and test data
  object_with_same_genes <- .create_same_feature_space_for_all_samples(object_to_annotate, ref_data_immgen, ref_data_haemoshere, ref_data_immgen2, ref_data_immgen_multiple_tissues, assay_to_use)
  sc_data <- object_with_same_genes[["sc_data"]]
  ref_data_immgen <- object_with_same_genes[["ref_data_immgen"]]
  ref_data_immgen2 <- object_with_same_genes[["ref_data_immgen2"]]
  ref_data_haemoshere <- object_with_same_genes[["ref_data_haemoshere"]]
  ref_data_immgen_multiple_tissues <- object_with_same_genes[["ref_data_immgen_multiple_tissues"]]
  
  # 3. do TPM normalization: divide by gene length and scale afterwards; this is necessary since 10X scRNAseq data are "TPM-normalized by nature" since molecules are counted by UMIs
  list_with_normalized_objects <- .TPM_normalization(genes_for_normalization=rownames(sc_data),
                                                     ref_data_immgen=ref_data_immgen,
                                                     ref_data_haemoshere=ref_data_haemoshere,
                                                     ref_data_immgen2=ref_data_immgen2,
                                                     ref_data_immgen_multiple_tissues=ref_data_immgen_multiple_tissues)
  ref_data_immgen <- list_with_normalized_objects[["ref_data_immgen"]]
  ref_data_immgen2 <- list_with_normalized_objects[["ref_data_immgen2"]]
  ref_data_haemoshere <- list_with_normalized_objects[["ref_data_haemoshere"]]
  ref_data_immgen_multiple_tissues <- list_with_normalized_objects[["ref_data_immgen_multiple_tissues"]]
  
  # 4. add pseudocount(1) and do log2 transformation
  # (necessary to allow gene selection within SingleR)
  ref_data_immgen_log2 <- log2(ref_data_immgen+1)
  ref_data_immgen2_log2 <- log2(ref_data_immgen2+1)
  ref_data_haemoshere_log2 <- log2(ref_data_haemoshere+1)
  ref_data_immgen_multiple_tissues_log2 <- log2(ref_data_immgen_multiple_tissues+1)
  
  # 5. merge immgen2 data into immgen dataset to allow for better marker selection within SingleR
  ref_data_immgen_log2 <- cbind(ref_data_immgen_log2, ref_data_immgen2_log2)
  
  # 6. rename samples
  renamed_tables <- .rename_samples(ref_data_immgen_log2=ref_data_immgen_log2,
                                    ref_data_haemoshere_log2=ref_data_haemoshere_log2,
                                    ref_data_immgen_multiple_tissues_log2=ref_data_immgen_multiple_tissues_log2)
  ref_data_immgen_log2=renamed_tables[["ref_data_immgen_log2"]]
  ref_data_haemoshere_log2=renamed_tables[["ref_data_haemoshere_log2"]]
  ref_data_immgen_multiple_tissues_log2=renamed_tables[["ref_data_immgen_multiple_tissues_log2"]]
  
  # 7. some samples are removed to avoid to much overlap between populations (e.g. HSC, LSK, MMP2/3/4)
  col_indices_to_keep_immgen <- which(! (colnames(ref_data_immgen_log2) %in% c("GN", "MMP2","MMP3", "MMP4")))
  col_indices_to_keep_haemosphere <- which(! (colnames(ref_data_haemoshere_log2) %in% c("LSK")))
  ref_data_immgen_log2_clean <- ref_data_immgen_log2[,col_indices_to_keep_immgen]
  ref_data_haemoshere_log2_clean <- ref_data_haemoshere_log2[,col_indices_to_keep_haemosphere]
  
  # 8. running singleR:
  annotation_results <- SingleR(test = sc_data,
                                ref = list(Immgen=ref_data_immgen_log2_clean, Haemoshere=ref_data_haemoshere_log2_clean, AllTissues=ref_data_immgen_multiple_tissues_log2),
                                labels = list(colnames(ref_data_immgen_log2_clean), colnames(ref_data_haemoshere_log2_clean), colnames(ref_data_immgen_multiple_tissues_log2)))
  # write cell labels into Seurat object and plot on umap
  annotation_results_sorted <- annotation_results[colnames(object_to_annotate@assays[[assay_to_use]]@counts),]
  object_to_annotate@meta.data[["SingleR_labels"]] <- annotation_results_sorted$pruned.labels
  # remove numbers in labels
  labels_with_numbers <- which(substr(object_to_annotate@meta.data[["SingleR_labels"]], start = nchar(object_to_annotate@meta.data[["SingleR_labels"]])-1, stop = nchar(object_to_annotate@meta.data[["SingleR_labels"]])-1) == ".")
  for (index in (1:length(labels_with_numbers))){
    position <- labels_with_numbers[index]
    object_to_annotate@meta.data[["SingleR_labels"]][position] <- substr(object_to_annotate@meta.data[["SingleR_labels"]][position], start=1, stop=nchar(object_to_annotate@meta.data[["SingleR_labels"]][position])-2)
  }
  # pool all T-cells
  object_to_annotate$SingleR_labels[which(object_to_annotate$SingleR_labels %in% c("T4", "Tgd", "T8", "NKT"))] <- "T"
  return(object_to_annotate)
}


# function to annotate cells more coarse grained
get_coarse_annotation <- function(object_to_annotate){
  object_to_annotate$SingleR_labels[which(object_to_annotate$SingleR_labels %in% c("LTHSC", "STHSC", "MPP"))] <- "stem cells"
  object_to_annotate$SingleR_labels[which(object_to_annotate$SingleR_labels %in% c("PreCFUE"))] <- "CFUE"
  object_to_annotate$SingleR_labels[which(object_to_annotate$SingleR_labels %in% c("Eo", "Neut"))] <- "Granulocytes"
  object_to_annotate$SingleR_labels[which(object_to_annotate$SingleR_labels %in% c("B", "proB", "CLP"))] <- "B-cells"
  object_to_annotate$SingleR_labels[which(object_to_annotate$SingleR_labels %in% c("Mono", "InfMono"))] <- "Monocytes"
  #DimPlot(object = object_to_annotate, reduction = "umap", label = TRUE, group.by = "SingleR_labels", repel = TRUE)+NoLegend()
  #DimPlot(object = object_to_annotate, reduction = "umap", label = TRUE, group.by = "seurat_clusters", repel = TRUE)+NoLegend()
  most_freq_cell_types <- c()
  for (index in (1:length(unique(object_to_annotate$seurat_clusters)))){
    current_cluster <- unique(object_to_annotate$seurat_clusters)[index]
    current_labels <- object_to_annotate$SingleR_labels[which(object_to_annotate$seurat_clusters==current_cluster)]
    frequencies <- table(current_labels)
    freq_vector <- as.vector(frequencies)
    names(freq_vector) <- names(frequencies)
    most_freq_cell_types <- append(most_freq_cell_types, names(freq_vector[order(freq_vector, decreasing = TRUE)])[1])
  }
  names(most_freq_cell_types) <- unique(object_to_annotate$seurat_clusters)
  B_cell_clusters <- names(most_freq_cell_types)[which(most_freq_cell_types=="B-cells")]
  object_to_annotate$SingleR_labels[which((!(object_to_annotate$seurat_clusters %in% B_cell_clusters)) & (object_to_annotate$SingleR_labels %in% c("B-cells")))] <- "CLP"
  T_cell_clusters <- names(most_freq_cell_types)[which(most_freq_cell_types=="T")]
  object_to_annotate$SingleR_labels[which(object_to_annotate$seurat_clusters %in% T_cell_clusters)] <- "T"
  NK_cell_clusters <- names(most_freq_cell_types)[which(most_freq_cell_types=="NK")]
  object_to_annotate$SingleR_labels[which(object_to_annotate$seurat_clusters %in% NK_cell_clusters)] <- "NK"
  return(object_to_annotate)
}
# function to visualize markers of clusters
visualize_top_marker_genes_of_clusters <- function(object_with_clusters, path_for_plot, number_of_top_markers_per_cluster=10, p_val_ad_thres=0.05, assay="SCT", slot="data", test.use = "wilcox"){
  markers <- FindAllMarkers(object_with_clusters, assay = assay, slot = slot, test.use = test.use)
  clusters_HSC <- as.character(unique(object_with_clusters$seurat_clusters))
  union_of_top_markers <- c()
  for (index in (1:length(clusters_HSC))){
    current_cluster <- clusters_HSC[index]
    markers_current_cluster <- markers[which((markers$cluster==current_cluster) &
                                               (markers$p_val_adj<p_val_ad_thres)),]
    sorted_genes <- markers_current_cluster[order(markers_current_cluster$avg_log2FC, decreasing = TRUE), "gene"]
    top_10 <- sorted_genes[1:number_of_top_markers_per_cluster]
    union_of_top_markers <- append(union_of_top_markers, top_10)
  }
  union_of_top_markers <- unique(union_of_top_markers)
  dotplot <- DotPlot(object = object_with_clusters, features = union_of_top_markers, assay = "SCT")+ coord_flip() + theme(axis.text.x = element_text(angle = 90))
  ggsave(
    path_for_plot,
    plot = dotplot,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 10,
    height = 15,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
}
# function creates a dotplot showing the expression of G2M and S-phase genes in cells classified by cell cycle phase
visualize_cell_cycle_gene_expression_in_LTHSCs <- function(LTHSC_combined=LTHSC_combined, folder_with_plots=folder_with_plots){
  LTHSC_combined$new_Phase <- factor(LTHSC_combined$new_Phase, levels=c("G1", "G1S", "S", "G2M"))
  dotplot <- DotPlot(object = LTHSC_combined, group.by = "new_Phase", features = c(s.genes, g2m.genes), assay = "SCT", dot.scale = 3)+ coord_flip() + theme(axis.text.x = element_text(angle = 90, size = 10), axis.text.y = element_text(size = 6))
  ggsave(
    paste(folder_with_plots, "LTs_cycle_genes_dotplot.pdf", sep="/"),
    plot = dotplot,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 10,
    height = 18,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
}
# function to identify cluster with cycling cells
cell_cycle <- function(object_WT, object_KO_with_ref_labels, path_for_plot){
  cyclin.genes <- rownames(object_WT)[grep("^Ccn[abde][0-9]$", rownames(object_WT))]
  
  # do heatmap with Pearson residuals
  object_WT@meta.data[["Cluster"]] <- object_WT$cell_cycle
  DefaultAssay(object_WT) <- "RNA"
  object_WT <- Seurat::SCTransform(object_WT, verbose = FALSE, residual.features=c(cyclin.genes, c("Mki67", "Top2a")))
  heat_WT <- DoHeatmap(object_WT, features = c(cyclin.genes, c("Mki67", "Top2a")), slot = "scale.data", assay = "SCT", group.by = "Cluster", combine = FALSE, raster = FALSE)
  
  object_KO_with_ref_labels@meta.data[["Cluster"]] <- object_KO_with_ref_labels$predicted.celltype
  DefaultAssay(object_KO_with_ref_labels) <- "RNA"
  object_KO_with_ref_labels <- Seurat::SCTransform(object_KO_with_ref_labels, verbose = FALSE, residual.features=c(cyclin.genes, c("Mki67", "Top2a")))
  heat_KO <- DoHeatmap(object_KO_with_ref_labels, features = c(cyclin.genes, c("Mki67", "Top2a")), slot = "scale.data", assay = "SCT", group.by = "Cluster", combine = FALSE, raster = FALSE)
  
  # do heatmaps with raw counts
  DefaultAssay(object_WT) <- "RNA"
  heat_WT_raw <- DoHeatmap(object_WT, features = c(cyclin.genes, c("Mki67", "Top2a")), slot = "counts", assay = "RNA", group.by = "Cluster", combine = FALSE, raster = FALSE)
  
  DefaultAssay(object_KO_with_ref_labels) <- "RNA"
  heat_KO_raw <- DoHeatmap(object_KO_with_ref_labels, features = c(cyclin.genes, c("Mki67", "Top2a")), slot = "counts", assay = "RNA", group.by = "Cluster", combine = FALSE, raster = FALSE)
  
  ggsave(
    paste(path_for_plot, "heat_map_WT.pdf", sep="/"),
    plot = heat_WT[[1]],
    device = NULL,
    path = NULL,
    scale = 1,
    width = 10,
    height = 10,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
  
  
  ggsave(
    paste(path_for_plot, "heat_map_KO.pdf", sep="/"),
    plot = heat_KO[[1]],
    device = NULL,
    path = NULL,
    scale = 1,
    width = 10,
    height = 10,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
  
  ggsave(
    paste(path_for_plot, "heat_map_WT_raw.pdf", sep="/"),
    plot = heat_WT_raw[[1]],
    device = NULL,
    path = NULL,
    scale = 1,
    width = 10,
    height = 10,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
  
  ggsave(
    paste(path_for_plot, "heat_map_KO_raw.pdf", sep="/"),
    plot = heat_KO_raw[[1]],
    device = NULL,
    path = NULL,
    scale = 1,
    width = 10,
    height = 10,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
  
  return(list(heat_WT=heat_WT, heat_KO=heat_KO, heat_WT_raw=heat_WT_raw, heat_KO_raw=heat_KO_raw))
}
# function for plotting umaps of LT-HSCs for cell cycle analysis
plot_umaps_cell_cycle <- function(WT, KO, path_for_plot){
  WT_plot <- DimPlot(object = WT, reduction = "umap", label = FALSE, group.by = "seurat_clusters", combine = FALSE, pt.size = 0.1)
  WT_plot <- WT_plot[[1]] + labs(color = "clusters")
  WT_plot[["labels"]][["title"]] <- "clustered WT-HSCs"
  WT_plot <- WT_plot + theme(plot.title = element_text(size=12), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
  WT_plot[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  
  WT_plot2 <- DimPlot(object = WT, reduction = "umap", label = FALSE, group.by = "cell_cycle", combine = FALSE, pt.size = 0.1)
  WT_plot2 <- WT_plot2[[1]] + labs(color = "cell cycle")
  WT_plot2[["labels"]][["title"]] <- "cell cycle classification \n in WT-HSCs"
  WT_plot2 <- WT_plot2 + theme(plot.title = element_text(size=12), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
  WT_plot2[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  
  KO_plot <- DimPlot(object = KO, reduction = "ref.umap", label = FALSE, group.by = "predicted.celltype", combine = FALSE, pt.size = 0.1)
  KO_plot <- KO_plot[[1]] + labs(color = "cell cycle")
  KO_plot[["labels"]][["title"]] <- "KO-HSCs with \n cell cycle annotation \n transfered from WT-HSCs"
  KO_plot <- KO_plot + theme(plot.title = element_text(size=12), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
  KO_plot[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  
  combined_plot <- WT_plot + WT_plot2 + KO_plot + plot_layout(guides = 'collect')
  
  ggsave(
    path_for_plot,
    plot = combined_plot,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 20,
    height = 7.5,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
}
plot_umaps_cell_cycle_updated <- function(WT, KO, path_for_plot){
  WT_plot <- DimPlot(object = WT, reduction = "umap_cycle", label = FALSE, group.by = "new_Phase", combine = FALSE, pt.size = 0.1)
  WT_plot <- WT_plot[[1]] + labs(color = "new_Phase")
  WT_plot[["labels"]][["title"]] <- "WT LT-HSC"
  WT_plot <- WT_plot + theme(plot.title = element_text(size=12), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
  WT_plot[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  
  KO_plot <- DimPlot(object = KO, reduction = "umap_cycle", label = FALSE, group.by = "new_Phase", combine = FALSE, pt.size = 0.1)
  KO_plot <- KO_plot[[1]] + labs(color = "new_Phase")
  KO_plot[["labels"]][["title"]] <- "KO LT-HSC"
  KO_plot <- KO_plot + theme(plot.title = element_text(size=12), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
  KO_plot[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  
  combined_plot <- WT_plot + KO_plot + plot_layout(guides = 'collect')
  
  ggsave(
    paste(path_for_plot, "LTs_cell_cycle_umap.pdf", sep="/"),
    plot = combined_plot,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 15,
    height = 7.5,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
}
# function that quantifies cells in a cluster and compares cell frequencies across conditions
quantify_cell_frequencies_by_cluster <- function(object_WT, object_KO_with_ref_labels, clusters_with_cycling_cells, path_for_plot){
  object_WT@meta.data[["Cluster"]] <- object_WT$cell_cycle
  object_KO_with_ref_labels@meta.data[["Cluster"]] <- object_KO_with_ref_labels$predicted.celltype
  
  object_WT$Cluster[which(object_WT$cell_cycle %in% clusters_with_cycling_cells)] <- "cycling"
  object_WT$Cluster[which(!(object_WT$cell_cycle %in% clusters_with_cycling_cells))] <- "non-cycling"
  
  object_KO_with_ref_labels$Cluster[which(object_KO_with_ref_labels$predicted.celltype %in% clusters_with_cycling_cells)] <- "cycling"
  object_KO_with_ref_labels$Cluster[which(!(object_KO_with_ref_labels$predicted.celltype %in% clusters_with_cycling_cells))] <- "non-cycling"
  
  percentage_non_cycling_WT <- length(which(object_WT$Cluster=="non-cycling"))/length(object_WT$Cluster)
  percentage_non_cycling_KO <- length(which(object_KO_with_ref_labels$Cluster=="non-cycling"))/length(object_KO_with_ref_labels$Cluster)
  df_plotting <- data.frame(condition=c("WT", "KO"),
                            percentage=c(percentage_non_cycling_WT, percentage_non_cycling_KO))
  df_plotting$condition <- factor(df_plotting$condition, levels = c("WT", "KO"))
  df_plotting$percentage <- df_plotting$percentage*100
  bar_plot <-ggplot(data=df_plotting, aes(x=condition, y=percentage, color=condition)) +
    geom_bar(stat="identity", fill="white") + theme_classic() + scale_color_manual(values=c("black", "red"))
  ggsave(
    path_for_plot,
    plot = bar_plot,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 5,
    height = 8,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
}
# compute log2(fold-changes) of cell type frequencies between conditions with 95% confidence intervalls; write pdf
quantify_compositional_changes_of_cell_types <- function(TBM_combined, path_for_plot="Figure_plot.pdf", pseudocount=0.01){
  # compute log2FCs for all populations with more than 100 cells in at least one condition
  vector_with_populations <- as.character(unique(TBM_combined$SingleR_labels))
  vector_with_populations <- vector_with_populations[!is.na(vector_with_populations)]
  selection_vector <- c()
  for (index in (1:length(vector_with_populations))){
    current_population <- vector_with_populations[index]
    if ((length(which((TBM_combined$SingleR_labels == current_population) & (TBM_combined$identity == "WT")))>=100) |
        (length(which((TBM_combined$SingleR_labels == current_population) & (TBM_combined$identity == "KO")))>=100)){
      selection_vector <- append(selection_vector, TRUE)
    } else {
      selection_vector <- append(selection_vector, FALSE)
    }
  }
  vector_with_large_populations <- vector_with_populations[selection_vector]
  
  log2FC <- c()
  for (index in (1:length(vector_with_large_populations))){
    current_population <- vector_with_large_populations[index]
    WT <- length(which((TBM_combined$SingleR_labels == current_population) & (TBM_combined$identity == "WT")))
    WT_freq <- WT/length(which(TBM_combined$identity == "WT"))
    KO <- length(which((TBM_combined$SingleR_labels == current_population) & (TBM_combined$identity == "KO")))
    KO_freq <- KO/length(which(TBM_combined$identity == "KO"))
    log2FC <- append(log2FC, log2((KO_freq+pseudocount)/(WT_freq+pseudocount)))
  }
  
  # compute confidence intervalls for foldchanges by bootstrapping
  num_repeats <- 10000
  foldchanges_after_sampling <- as.data.frame(matrix(0, ncol = length(vector_with_large_populations),
                                                     nrow = num_repeats))
  colnames(foldchanges_after_sampling) <- vector_with_large_populations
  for (index in (1:num_repeats)){
    WT_labels <- TBM_combined$SingleR_labels[which(TBM_combined$identity=="WT")]
    KO_labels <- TBM_combined$SingleR_labels[which(TBM_combined$identity=="KO")]
    WT_samples <- sample(WT_labels, size=length(WT_labels), replace = TRUE)
    KO_samples <- sample(KO_labels, size=length(KO_labels), replace = TRUE)
    for (index2 in (1:length(vector_with_large_populations))){
      current_label <- vector_with_large_populations[index2]
      perc_WT <- length(which(WT_samples == current_label)) / length(WT_samples)
      perc_KO <- length(which(KO_samples == current_label)) / length(KO_samples)
      foldchanges_after_sampling[index, current_label] <- log2((perc_KO+pseudocount)/(perc_WT+pseudocount))
    }
  }
  
  confidence_vector_lower <- c()
  confidence_vector_upper <- c()
  for (index in (1:length(vector_with_large_populations))){
    confidence_vector_lower <- append(confidence_vector_lower, quantile(as.numeric(foldchanges_after_sampling[, index]), c(.975, .025), na.rm=TRUE)[2])
    confidence_vector_upper <- append(confidence_vector_upper, quantile(as.numeric(foldchanges_after_sampling[, index]), c(.975, .025), na.rm=TRUE)[1])
  }
  
  # plot fold changes
  plotting_log2fc <- data.frame(Cell_Type=vector_with_large_populations, Log2FC=log2FC, lower=confidence_vector_lower, upper=confidence_vector_upper)
  plotting_log2fc <- plotting_log2fc[order(plotting_log2fc$Log2FC, decreasing = FALSE),]
  plotting_log2fc$Cell_Type <- factor(plotting_log2fc$Cell_Type, levels=plotting_log2fc$Cell_Type)
  
  bar_plot <- ggplot(data=plotting_log2fc, aes(x=Cell_Type, y=Log2FC)) +
    geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                  position=position_dodge(.9)) +
    theme(axis.text.x = element_text(angle = 45))
  
  ggsave(
    path_for_plot,
    plot = bar_plot,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 10,
    height = 5,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
}
# prepare LSK and Linminus cells for diffusion modeling: remove unnecessary cell populations
prepare_Linminus_for_diffusion_modeling <- function(Linminus_combined){
  # only retain cells in clusters of early progenitors
  Linminus_combined <- Seurat::FindNeighbors(Linminus_combined, reduction="pca", dims = 1:20, k.param = 50)
  Linminus_combined <- Seurat::FindClusters(Linminus_combined, resolution = 0.1)
  Linminus_combined_core <- subset(Linminus_combined, subset = seurat_clusters %in% c("0", "1", "5"))
  # re-compute SCT for Linminus
  Linminus_combined_core@assays[["SCT"]] <- NULL
  DefaultAssay(Linminus_combined_core) <- "RNA"
  Linminus_combined_core <- Seurat::SCTransform(Linminus_combined_core, verbose = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)
  DefaultAssay(Linminus_combined_core) <- "SCT"
  Linminus_combined_core <- Seurat::RunPCA(Linminus_combined_core, features = Seurat::VariableFeatures(object = Linminus_combined_core))
  return(Linminus_combined_core)
}

prepare_LSK_for_diffusion_modeling <- function(LSK_combined){
  LSK_combined <- subset(LSK_combined, subset = SingleR_labels %in% c("LTHSC", "STHSC", "MPP", "CMP", "GMP", "DC","Baso", "MEP", "Mono", "PreCFUE", "CLP", "proB", "preT"))
  return(LSK_combined)
}

merge_LSK_Linminus_using_Linminus_features <- function(LSK_combined, Linminus_combined_core){
  DefaultAssay(LSK_combined) <- "RNA"
  DefaultAssay(Linminus_combined_core) <- "RNA"
  Linminus_core_LSK <- merge(LSK_combined, y = c(Linminus_combined_core), add.cell.ids = c("LSK", "Linminus"), project = "Usp22", merge.data=FALSE)
  identity_vector <- c(rep("LSK", ncol(LSK_combined)), rep("Linminus", ncol(Linminus_combined_core)))
  Linminus_core_LSK@meta.data[["Sample"]] <- identity_vector
  Linminus_core_LSK@assays[["SCT"]] <- NULL
  # note that Linminus and LSK together are re-normalized based on a normalization model trained on Linminus only (this is done to avoid dominance of features related to cell cycle)
  DefaultAssay(Linminus_core_LSK) <- "RNA"
  DefaultAssay(Linminus_combined_core) <- "SCT"
  Linminus_core_LSK <- Seurat::SCTransform(Linminus_core_LSK, verbose = FALSE, reference.SCT.model = Linminus_combined_core@assays[["SCT"]]@SCTModel.list[["model1"]], residual.features = VariableFeatures(Linminus_combined_core))
  DefaultAssay(Linminus_core_LSK) <- "SCT"
  Linminus_core_LSK <- Seurat::RunPCA(Linminus_core_LSK, features = Seurat::VariableFeatures(object = Linminus_core_LSK))
  return(Linminus_core_LSK)
}

scatter_3d <- function(sce, dim.red.name = 'DiffusionMap', color.column, marker_size = 20, scene = ''){
  require(plotly)
  df = as.tibble(reducedDim(sce,dim.red.name)[,c(1:3)])
  colnames(df) = c('dr1','dr2','dr3')
  if(color.column %in% colnames(colData(sce))){
    df$color = colData(sce)[,color.column]
  }else{
    df$color = logcounts(sce)[color.column,]    
  }
  p <- plot_ly(df, x = ~dr1, y = ~dr2, z = ~dr3, color = ~color, scene = scene) %>%
    add_markers(size = marker_size) %>%
    layout(scene = list(xaxis = list(title = 'x'),
                        yaxis = list(title = 'y'),
                        zaxis = list(title = 'z')))
  return(p)
}

run_diffusion_analysis <- function(Linminus_core_LSK){
  # compute diffusion map ####
  Linminus_core_LSK.sce <- as.SingleCellExperiment(Linminus_core_LSK, assay = "SCT")
  dmap <- destiny::DiffusionMap(reducedDim(Linminus_core_LSK.sce,'PCA')[,1:30])
  reducedDim(Linminus_core_LSK.sce, 'DiffusionMap') <- dmap@eigenvectors
  
  # visualize intitial diffusion map
  initial_diff_map_plot_cell_types <- scatter_3d(Linminus_core_LSK.sce, dim.red.name = 'DiffusionMap', color.column="SingleR_labels", marker_size = 20, scene = '')
  initial_diff_map_plot_identity <- scatter_3d(Linminus_core_LSK.sce, dim.red.name = 'DiffusionMap', color.column="identity", marker_size = 20, scene = '')
  
  # remove outliers and cell types that are not of interest from diffusion map ####
  Linminus_core_LSK.sce$y_value <- dmap@eigenvectors[colnames(Linminus_core_LSK.sce), 2]
  Linminus_core_LSK.sce$outlier <- (Linminus_core_LSK.sce$y_value < (-0.01))
  outlier_plot <- scatter_3d(Linminus_core_LSK.sce, dim.red.name = 'DiffusionMap', color.column="outlier", marker_size = 20, scene = '')
  Linminus_core_LSK.sce <- Linminus_core_LSK.sce[, !Linminus_core_LSK.sce$outlier]
  
  # re-compute diffusion_map
  dmap <- destiny::DiffusionMap(reducedDim(Linminus_core_LSK.sce,'PCA')[,1:30])
  reducedDim(Linminus_core_LSK.sce, 'DiffusionMap') <- dmap@eigenvectors
  diff_map_after_outlier_removal <- scatter_3d(Linminus_core_LSK.sce, dim.red.name = 'DiffusionMap', color.column="SingleR_labels", marker_size = 20, scene = '')
  
  # keep only cell types of interest for both subgroups
  Interesting_cells <- Linminus_core_LSK.sce[, Linminus_core_LSK.sce$SingleR_labels %in% c("GMP", "CMP", "PreCFUE", "CLP", "MPP", "LTHSC", "STHSC", "CFUE", "MEP")]
  final_diffusion_map <- scatter_3d(Interesting_cells, dim.red.name = 'DiffusionMap', color.column="SingleR_labels", marker_size = 20, scene = '')
  return(Interesting_cells)
}

visualize_samples_and_cell_types_on_diffusion_map <- function(Interesting_cells, path_for_plot="current_Fig.pdf"){
  Interesting_cells.seurat <- as.Seurat(Interesting_cells)
  Fig1 <- DimPlot(Interesting_cells.seurat, reduction = "DiffusionMap", dims = c(1, 2), label = FALSE, group.by = "SingleR_labels", repel = TRUE, combine = FALSE)
  Fig1 <- Fig1[[1]]
  Fig1[["labels"]][["title"]] <- "Annotation by SingleR"
  Fig1[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  Fig2 <- DimPlot(Interesting_cells.seurat, reduction = "DiffusionMap", dims = c(1, 2), label = FALSE, group.by = "Sample", repel = TRUE, cols = c("grey", "red"),  combine = FALSE)
  Fig2 <- Fig2[[1]]
  Fig2[["labels"]][["title"]] <- "Sample Types"
  Fig2[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  
  # subset to Linminus and LSK
  Linminus_diff <- Interesting_cells.seurat[, Interesting_cells.seurat$Sample=="Linminus"]
  LSK_diff <- Interesting_cells.seurat[, Interesting_cells.seurat$Sample=="LSK"]
  
  # subsample WT cells to allow comparable visualization
  num_of_samples <- length(which(LSK_diff$identity=="KO"))
  WT_indices <- which(LSK_diff$identity=="WT")
  WT_sampled <- sample(WT_indices, size=num_of_samples, replace = FALSE)
  LSK_diff_corrected <- LSK_diff[, c(which(LSK_diff$identity=="KO"), WT_sampled)]
  Fig3 <- DimPlot(LSK_diff_corrected, reduction = "DiffusionMap", dims = c(1, 2), label = FALSE, group.by = "identity", repel = TRUE, combine = FALSE, cols = c('KO' = 'red', 'WT' = 'black'))
  Fig3 <- Fig3[[1]]
  Fig3[["labels"]][["title"]] <- "Conditions"
  Fig3[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  
  presentation_figure1 <- gridExtra::arrangeGrob(
    grobs = list(Fig1, Fig2, Fig3),
    widths = c(1),
    layout_matrix = rbind(c(1), c(2), c(3))
  )
  
  ggsave(
    path_for_plot,
    plot = presentation_figure1,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 15,
    height = 30,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
}

define_cluster_identities_on_dmap_for_slingshot_analysis <- function(Interesting_cells){
  Interesting_cells.seurat <- as.Seurat(Interesting_cells)
  
  # cluster cells to check cluster identity of stem cell clusters and define root cluster for slingshot analysis
  Interesting_cells.seurat <- cluster_cells(Interesting_cells.seurat, reduction = "DiffusionMap", ndims=3, knn=50, create_neighbor_object=FALSE, resolution=0.2)
  Interesting_cells@colData@listData[["Cluster"]] <- Interesting_cells.seurat$seurat_clusters
  dmap_plot_clusters <- scatter_3d(Interesting_cells, dim.red.name = 'DiffusionMap', color.column="Cluster", marker_size = 20, scene = '')
  dmap_plot_cell_types <- scatter_3d(Interesting_cells, dim.red.name = 'DiffusionMap', color.column="SingleR_labels", marker_size = 20, scene = '')
  dmap_plot_samples <- scatter_3d(Interesting_cells, dim.red.name = 'DiffusionMap', color.column="Sample", marker_size = 20, scene = '')
  # check markers for the two different stem cell clusters
  clusters_for_DEA <- subset(Interesting_cells.seurat, subset= seurat_clusters %in% c(3, 5))
  Idents(clusters_for_DEA) <- clusters_for_DEA$seurat_clusters
  markers <- FindAllMarkers(clusters_for_DEA)
  # cluster 3 is used as the root (Ifitm1, Hlf, Ly6a, Cd74); cluster 5 consists of more proliferating pf4 positive stem cells
  root_cells <- which(as.character(Interesting_cells$Cluster)=="3")
  
  # do second clustering for defining CLP and GMP branches and the cluster in which they separate
  Interesting_cells.seurat <- cluster_cells(Interesting_cells.seurat, reduction = "DiffusionMap", ndims=3, knn=50, create_neighbor_object=FALSE, resolution=0.1)
  Interesting_cells@colData@listData[["Cluster2"]] <- as.character(Interesting_cells.seurat$seurat_clusters)
  dmap_plot_branch_clusters <- scatter_3d(Interesting_cells, dim.red.name = 'DiffusionMap', color.column="Cluster2", marker_size = 20, scene = '')
  # subcluster cluster3 and visualize
  Interesting_cells.seurat_3 <- subset(Interesting_cells.seurat, subset=seurat_clusters %in% c("3"))
  Interesting_cells.seurat_3 <- cluster_cells(Interesting_cells.seurat_3, reduction = "DiffusionMap", ndims=3, knn=50, create_neighbor_object=FALSE, resolution=0.1)
  Interesting_cells_3 <- Interesting_cells[,which(Interesting_cells$Cluster2=="3")]
  Interesting_cells_3@colData@listData[["Cluster2"]] <- Interesting_cells.seurat_3$seurat_clusters
  Interesting_cells_3$Cluster2 <- as.character(Interesting_cells_3$Cluster2)
  dmap_plot_subclustering_cluster3 <- scatter_3d(Interesting_cells_3, dim.red.name = 'DiffusionMap', color.column="Cluster2", marker_size = 20, scene = '')
  Interesting_cells_3$Cluster2[which(Interesting_cells_3$Cluster2 %in% c("0", "1", "3"))] <- "GMP_branch"
  Interesting_cells_3$Cluster2[which(Interesting_cells_3$Cluster2 %in% c("2"))] <- "CLP_branch"
  Interesting_cells$Cluster2[which(Interesting_cells$Cluster2=="3")] <- Interesting_cells_3$Cluster2
  Interesting_cells$Cluster2[which(Interesting_cells$Cluster2 %in% c("0", "1", "2"))] <- "split"
  
  # do third clustering to identify endpoints for slingshot
  Interesting_cells.seurat <- Seurat::FindNeighbors(Interesting_cells.seurat, reduction="DiffusionMap", dims = 1:3, k.param = 20)
  Interesting_cells.seurat <- Seurat::FindClusters(Interesting_cells.seurat, resolution = 0.8)
  Interesting_cells@colData@listData[["Cluster3"]] <- Interesting_cells.seurat$seurat_clusters
  dmap_plot_clusters3 <- scatter_3d(Interesting_cells, dim.red.name = 'DiffusionMap', color.column="Cluster3", marker_size = 20, scene = '')
  CFUEs <- which(as.character(Interesting_cells$Cluster3)=="22")
  GMPs <- which(as.character(Interesting_cells$Cluster3) %in% c("25", "7", "24"))
  CLPs <- which(as.character(Interesting_cells$Cluster3) %in% c("27"))
  Baso <- which(as.character(Interesting_cells$Cluster3) %in% c("26"))
  
  # annotate subclusters of cluster3 and merge annotations into complete object
  Interesting_cells@colData@listData[["location"]] <- Interesting_cells$Cluster2
  Interesting_cells$location[root_cells] <- "root_cells"
  Interesting_cells$location[CFUEs] <- "CFUEs"
  Interesting_cells$location[GMPs] <- "GMPs"
  Interesting_cells$location[CLPs] <- "CLPs"
  Interesting_cells$location[Baso] <- "Baso"
  
  dmap_plot_with_clusters_for_slingshot <- scatter_3d(Interesting_cells, dim.red.name = 'DiffusionMap', color.column="location", marker_size = 20, scene = '')
  return(list(object_with_endpoints=Interesting_cells, 
              dmap_plot_with_clusters_for_slingshot=dmap_plot_with_clusters_for_slingshot,
              dmap_plot_clusters1=dmap_plot_clusters,
              dmap_plot_clusters2=dmap_plot_branch_clusters,
              dmap_plot_clusters3=dmap_plot_clusters3,
              dmap_plot_cell_types=dmap_plot_cell_types,
              dmap_plot_samples=dmap_plot_samples,
              markers_stem_cells=markers,
              dmap_plot_subclustering_of_2nd_clustering=dmap_plot_subclustering_cluster3))
}

annotate_cells_by_assigned_slingshot_trajectories_and_rename_pseudotimes <- function(Interesting_cells){
  trajectory_assignment <- slingBranchID(Interesting_cells)
  trajectory_assignment <- as.character(trajectory_assignment)
  trajectory_assignment[which(trajectory_assignment=="1")] <- "CLP_lineage"
  trajectory_assignment[which(trajectory_assignment=="2")] <- "GMP_lineage"
  trajectory_assignment[which(trajectory_assignment=="3")] <- "CFUE_lineage"
  trajectory_assignment[which(trajectory_assignment=="4")] <- "Baso_lineage"
  trajectory_assignment[which(trajectory_assignment=="1,2")] <- "CLP_GMP"
  trajectory_assignment[which(trajectory_assignment=="1,2,3,4")] <- "shared"
  trajectory_assignment[which(trajectory_assignment=="3,4")] <- "CFUE_Baso"
  Interesting_cells@colData@listData[["Lineage_assignment"]] <- trajectory_assignment
  Interesting_cells@colData@listData[["pseudotime_CLP"]] <- Interesting_cells$slingPseudotime_1
  Interesting_cells@colData@listData[["pseudotime_GMP"]] <- Interesting_cells$slingPseudotime_2
  Interesting_cells@colData@listData[["pseudotime_CFUE"]] <- Interesting_cells$slingPseudotime_3
  Interesting_cells@colData@listData[["pseudotime_Baso"]] <- Interesting_cells$slingPseudotime_4
  return(Interesting_cells)
}

plot_dmap_with_slingshot_trajectories <- function(Interesting_cells, color_cells_by="Lineage_assignment", path_for_plot="presentation_figure1.pdf"){
  Seurat_slingshot <- as.Seurat(Interesting_cells)
  Fig1 <- DimPlot(Seurat_slingshot, reduction = "DiffusionMap", dims = c(1, 2), label = FALSE, group.by = color_cells_by, repel = TRUE, combine = FALSE)
  Fig1 <- Fig1[[1]]
  # add trajectories computed by slingshot
  order_of_projections1 <- Interesting_cells@int_metadata[["slingshot"]]@curves[["curve1"]][["ord"]]
  order_of_projections2 <- Interesting_cells@int_metadata[["slingshot"]]@curves[["curve2"]][["ord"]]
  order_of_projections3 <- Interesting_cells@int_metadata[["slingshot"]]@curves[["curve3"]][["ord"]]
  order_of_projections4 <- Interesting_cells@int_metadata[["slingshot"]]@curves[["curve4"]][["ord"]]
  Fig1 <- Fig1 + geom_path(x=Interesting_cells@int_metadata[["slingshot"]]@curves[["curve1"]][["s"]][order_of_projections1,1],
                           y=Interesting_cells@int_metadata[["slingshot"]]@curves[["curve1"]][["s"]][order_of_projections1,2]) +
    geom_path(x=Interesting_cells@int_metadata[["slingshot"]]@curves[["curve2"]][["s"]][order_of_projections2,1],
              y=Interesting_cells@int_metadata[["slingshot"]]@curves[["curve2"]][["s"]][order_of_projections2,2])+
    geom_path(x=Interesting_cells@int_metadata[["slingshot"]]@curves[["curve3"]][["s"]][order_of_projections3,1],
              y=Interesting_cells@int_metadata[["slingshot"]]@curves[["curve3"]][["s"]][order_of_projections3,2])+
    geom_path(x=Interesting_cells@int_metadata[["slingshot"]]@curves[["curve4"]][["s"]][order_of_projections4,1],
              y=Interesting_cells@int_metadata[["slingshot"]]@curves[["curve4"]][["s"]][order_of_projections4,2])
  
  Fig1[["labels"]][["title"]] <- "Trajectories inferred by Slingshot"
  Fig1[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  
  ggsave(
    path_for_plot,
    plot = Fig1,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 15,
    height = 10,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
  
}

plot_dmap_with_slingshot_trajectories_LSK <- function(Interesting_cells, color_cells_by="identity", path_for_plot="presentation_figure1.pdf"){
  Seurat_slingshot <- as.Seurat(Interesting_cells)
  Seurat_slingshot <- subset(Seurat_slingshot, subset=Sample=="LSK")
  SCE_slingshot <- Interesting_cells[,which(Interesting_cells$Sample=="LSK")]
  Fig1 <- DimPlot(Seurat_slingshot, reduction = "DiffusionMap", dims = c(1, 2), label = FALSE, group.by = "identity", repel = TRUE, combine = FALSE, shuffle=TRUE, cols = c('KO' = 'red', 'WT' = 'black'))
  Fig1 <- Fig1[[1]]
  # add trajectories computed by slingshot
  order_of_projections2 <- SCE_slingshot@int_metadata[["slingshot"]]@curves[["curve2"]][["ord"]]
  order_of_projections3 <- SCE_slingshot@int_metadata[["slingshot"]]@curves[["curve3"]][["ord"]]
  projections_curve2_x <- SCE_slingshot@int_metadata[["slingshot"]]@curves[["curve2"]][["s"]][order_of_projections2,1]
  projections_curve2_y <- SCE_slingshot@int_metadata[["slingshot"]]@curves[["curve2"]][["s"]][order_of_projections2,2]
  projections_curve3_x <- SCE_slingshot@int_metadata[["slingshot"]]@curves[["curve3"]][["s"]][order_of_projections3,1]
  projections_curve3_y <- SCE_slingshot@int_metadata[["slingshot"]]@curves[["curve3"]][["s"]][order_of_projections3,2]
  
  projections_curve2_x <- projections_curve2_x[which(names(projections_curve2_x) %in% colnames(SCE_slingshot))]
  projections_curve2_y <- projections_curve2_y[which(names(projections_curve2_y) %in% colnames(SCE_slingshot))]
  projections_curve3_x <- projections_curve3_x[which(names(projections_curve3_x) %in% colnames(SCE_slingshot))]
  projections_curve3_y <- projections_curve3_y[which(names(projections_curve3_y) %in% colnames(SCE_slingshot))]
  
  
  Fig1 <- Fig1 + geom_path(x=projections_curve2_x,
                           y=projections_curve2_y)+
    geom_path(x=projections_curve3_x,
              y=projections_curve3_y)
  Fig1 <- Fig1 + annotate("text", x=0.001, y=0.008, label= "Myeloid Traj.") +
    annotate("text", x=0.008, y=0.003, label= "MEP Traj.") +
    annotate("text", x=0.0087, y=-0.0023, label= "HSCs")
  
  Fig1[["labels"]][["title"]] <- "Trajectories inferred by Slingshot"
  Fig1[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  
  ggsave(
    path_for_plot,
    plot = Fig1,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 15,
    height = 10,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
  
}


run_condition_imbalance_analysis <- function(Interesting_cells, sample_to_use="LSK", pseudo_time_trajectories_to_use=c("pseudotime_GMP", "pseudotime_CFUE"), 
                                             size_of_local_neighborhoods=50, num_bootstrap_iter=2, path_to_helper_functions="condition_imbalance_analysis_helper_functions.R", 
                                             size_moving_window=0.05){
  source(path_to_helper_functions, local = TRUE)
  
  # 1. before sampling the neighbors of each cell are identified and stored in a matrix called "nn_lists_in_rows";
  # element (i,j) contains the index of the j-th-nearest neighbor of cell i
  Interesting_cells.Seurat <- as.Seurat(Interesting_cells)
  seurat_sample <- subset(Interesting_cells.Seurat, subset=Sample==sample_to_use)
  # The factor used to define the k.param of the FindNeighbors function might need to be increased for small neighborhoods
  seurat_sample <- Seurat::FindNeighbors(seurat_sample, reduction="DiffusionMap_subset", dims = 1:3, k.param = size_of_local_neighborhoods*2, return.neighbor = TRUE)
  nn_lists_in_rows <- seurat_sample@neighbors[["RNA.nn"]]@nn.idx
  Indices_WT <- which(seurat_sample$identity=="WT")
  Indices_KO <- which(seurat_sample$identity=="KO")
  
  # 2. Bootstrap-samples from WT and KO cells are created in multiple iterations; In each iteration the same number of cells are sampled from both conditions (using the minimum size of both conditions);
  # At the end of each iteration, a vector is created: the j-th element describes how often the cell with index j was sampled;
  # The obtained vectors are used as rows of a matrix so that each row captures the outcome of one iteration
  matrix_samples <- .create_matrix_with_number_of_observations_per_cell___bootstrap_iterations_correspond_to_rows(seurat_sample, Indices_WT, Indices_KO, sample_size=min(c(length(Indices_WT), length(Indices_KO))), num_bootstrap_iter)
  
  # 3. A matrix with the same dimensions as in (2.) is created: Element(i,j) gives the proportion of KO cells in the neighborhood of cell j in the i-th iteration after sampling;
  # For finding the nearest neighbors of cell j in sampling iteration i, it is considered that each neighbor in the cell's non-sampled neighborhood can after sampling be either missing or occur once or multiple times;
  # For example: nearest neighbors before sampling: cellA, cellB, cellC; after sampling: cellA, cellB, cellB
  # Therefor, the number of "truly" distinct cells contributing to a neighborhood can be smaller than the defined size of the neighborhood 
  # but the sum of the observations of all neighbors must equal the size of the neighborhood.
  matrix_with_scores <- .improved___compute_matrix_with_proportion_of_WT_cells_in_neighborhood___bootstrap_iterations_correspond_to_rows(seurat_sample, matrix_samples, num_bootstrap_iter, nn_lists_in_rows, size_of_local_neighborhoods)
  
  # 4. Use the matrices from the two previous steps to construct a list containing one dataframe per bootstrap iteration;
  # Each dataframe (with rows corresponding to cells) contains a column with scores (proportion of KO cells in the neighborhood) and another one with the corresponding pseudo-time along a specified trajectory. 
  # As some cells were sampled multiple times, the dataframes contain duplicated rows.
  # In each dataframe, the rows (cells) are ordered by increasing pseudo-time values and rows with missing pseudo-time values are removed. 
  # Those rows correspond to cells which were not assigned to the trajectory of interest.
  # For each row a moving average of the score is computed. For calculating the moving average of the score at pseudo-time X, all scores within a defined range of pseudo-time values are used.
  # After calculating the moving average, the rows are made unique and scores at defined pseudo-time values are interpolated.
  # The interpolated values from the different iteration (dataframes in the list) are written into columns of a dataframe (columns corresponding to iterations)
  # and for each row (pseudo-time value) the median and 95% confidence intervall is computed.
  # The whole procedure needs to be done for each trajectory separately as each cell can be assigned to one or multiple trajectories.
  list_with_trajectories_and_their_medians_and_CIs_ordered_by_pseudotime <- list()
  for (index in (1:length(pseudo_time_trajectories_to_use))){
    pseudo_time_trajectory_to_use <- pseudo_time_trajectories_to_use[index]
    
    list_sampling_results <- .create_list_with_dataframes___pseudotime_and_score_columns___each_dataframe_corresponds_to_a_bootstrap_iteration(seurat_sample=seurat_sample, 
                                                                                                                                               matrix_with_scores=matrix_with_scores, 
                                                                                                                                               matrix_samples=matrix_samples, 
                                                                                                                                               pseudo_time_trajectory_to_use=pseudo_time_trajectory_to_use)
    
    list_with_interpolated_values <- .interpolate_values(list_sampling_results=list_sampling_results, 
                                                         seurat_sample=seurat_sample,
                                                         pseudo_time_trajectory_to_use=pseudo_time_trajectory_to_use,
                                                         size_moving_window=size_moving_window)
    
    
    df_confidence <- .compute_median_and_ci_of_interpolated_values(list_with_interpolated_values=list_with_interpolated_values)
    list_with_trajectories_and_their_medians_and_CIs_ordered_by_pseudotime <- list.append(list_with_trajectories_and_their_medians_and_CIs_ordered_by_pseudotime, df_confidence)
  }
  names(list_with_trajectories_and_their_medians_and_CIs_ordered_by_pseudotime) <- pseudo_time_trajectories_to_use
  return(list_with_trajectories_and_their_medians_and_CIs_ordered_by_pseudotime)
}

plot_condition_imbalance_analysis <- function(df_confidence, x_label="pseudo-time along myeloid trajectory", y_label="|KO| / (|KO|+|WT|)", path_for_plot="presentation_figure1.pdf"){
  plot <- ggplot(df_confidence, aes(x = pseudo_time, y = median)) +
    #geom_line() +
    geom_line(aes(y = median, color="Median Moving Average", linetype="Median Moving Average")) + 
    geom_ribbon(aes(ymin = lower, ymax = upper, color="95% Confidence Intervall", linetype="95% Confidence Intervall"), alpha = 0.5) +
    geom_line(aes(y = 0.5, color="|WT| = |KO|", linetype="|WT| = |KO|")) + ylab(y_label)+ xlab(x_label) + ylim(0,1) +
    labs(linetype = "statistics",
         color = "statistics") +
    scale_color_manual(name = "statistics", values = c("Median Moving Average" = "black", "95% Confidence Intervall" = "grey", "|WT| = |KO|" = "blue")) +
    scale_linetype_manual(name = "statistics",
                          #breaks = c("Median Moving Average", "|WT| = |KO|"),
                          values = c("Median Moving Average" = "dotted", "95% Confidence Intervall" = "solid", "|WT| = |KO|" = "solid"))
  
  ggsave(
    path_for_plot,
    plot = plot,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 12,
    height = 6.5,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
}

plot_Fig1C <- function(TBM_combined, path_for_plot=folder_with_plots){
  DefaultAssay(TBM_combined) <- "SCT"
  TBM_combined_plotting_KO <- subset(TBM_combined, subset = identity == "KO")
  TBM_combined_plotting_WT <- subset(TBM_combined, subset = identity == "WT")
  test1 <- FeaturePlot(TBM_combined_plotting_WT, features = c("Usp22"), slot = "data", combine = FALSE, cols = c("lightgrey", "#CC0000"))
  test1 <- test1[[1]]
  test1[["labels"]][["title"]] <- "Usp22 mRNA in WT"
  test1[["layers"]][[1]][["aes_params"]][["alpha"]] <- 1.0
  limits <- test1[["scales"]][["scales"]][[3]][["limits"]]
  test2 <- FeaturePlot(TBM_combined_plotting_KO, features = c("Usp22"), slot = "data", combine = FALSE, cols = c("lightgrey", "#CC0000"))
  test2 <- test2[[1]]
  test2[["scales"]][["scales"]][[3]][["limits"]] <- limits
  test2[["labels"]][["title"]] <- "Usp22 mRNA in KO"
  test2[["layers"]][[1]][["aes_params"]][["alpha"]] <- 1.0
  annotation_plot <- DimPlot(object = TBM_combined, reduction = "umap", label = TRUE, group.by = "SingleR_labels", repel = TRUE, combine = FALSE)
  test3 <- annotation_plot[[1]] + NoLegend()
  test3[["labels"]][["title"]] <- "Cell Types"
  test3[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  condition_plot <- DimPlot(object = TBM_combined, reduction = "umap", label = FALSE, group.by = "identity", combine = FALSE, cols = c('KO' = 'red', 'WT' = 'black'), shuffle = TRUE)
  test4 <- condition_plot[[1]]
  test4[["labels"]][["title"]] <- "Condition"
  test4[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.2
  #gridExtra::grid.arrange(test1, test2, nrow = 1, top = "Normalized expression of Usp22 mRNA")
  
  TBM_paper_figure <- gridExtra::arrangeGrob(
    grobs = list(test1, test2, test3, test4),
    widths = c(1, 1),
    layout_matrix = rbind(c(3, 4),
                          c(1, 2))
  )
  
  combined_plot_upper <- test3 + test4 + plot_layout(guides = 'collect')
  combined_plot_lower <- test1 + test2 + plot_layout(guides = 'collect')
  
  ggsave(
    paste(path_for_plot, "Fig1C_panel.pdf", sep="/"),
    plot = TBM_paper_figure,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 20,
    height = 20,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
  
  ggsave(
    paste(path_for_plot, "TBM_annotation_condition.pdf", sep="/"),
    plot = combined_plot_upper,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 20,
    height = 10,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
  
  ggsave(
    paste(path_for_plot, "TBM_usp22_expression.pdf", sep="/"),
    plot = combined_plot_lower,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 20,
    height = 10,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
  
}

create_usp22_panel_for_all_stem_cell_populations <- function(HSC_WT, ST_WT, MPP_WT, HSC_KO, ST_KO, MPP_KO, path_for_plot){
  DefaultAssay(HSC_WT) <- "SCT"
  DefaultAssay(HSC_KO) <- "SCT"
  HSC1 <- FeaturePlot(HSC_WT, features=c("Usp22"), slot = "data", combine = FALSE, pt.size = 0.3)
  HSC1 <- HSC1[[1]] + xlab("umap1") + ylab("umap2") + theme(plot.title = element_blank(), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
  HSC2 <- FeaturePlot(HSC_KO, features=c("Usp22"), slot = "data", combine = FALSE, reduction = "ref.umap", pt.size = 0.3)
  HSC2 <- HSC2[[1]] + xlab("umap1") + ylab("umap2") + theme(plot.title = element_blank(), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
  HSC1[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  HSC2[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  limits <- HSC1[["scales"]][["scales"]][[3]][["limits"]]
  HSC2[["scales"]][["scales"]][[3]][["limits"]] <- limits
  
  
  DefaultAssay(ST_WT) <- "SCT"
  DefaultAssay(ST_KO) <- "SCT"
  ST1 <- FeaturePlot(ST_WT, features=c("Usp22"), slot = "data", combine = FALSE, pt.size = 0.3)
  ST1 <- ST1[[1]] + xlab("umap1") + ylab("umap2") + theme(plot.title = element_blank(), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
  ST2 <- FeaturePlot(ST_KO, features=c("Usp22"), slot = "data", combine = FALSE, reduction = "ref.umap", pt.size = 0.3)
  ST2 <- ST2[[1]] + xlab("umap1") + ylab("umap2") + theme(plot.title = element_blank(), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
  ST1[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  ST2[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  limits <- ST1[["scales"]][["scales"]][[3]][["limits"]]
  ST2[["scales"]][["scales"]][[3]][["limits"]] <- limits
  
  
  DefaultAssay(MPP_WT) <- "SCT"
  DefaultAssay(MPP_KO) <- "SCT"
  MPP1 <- FeaturePlot(MPP_WT, features=c("Usp22"), slot = "data", combine = FALSE, pt.size = 0.3)
  MPP1 <- MPP1[[1]] + xlab("umap1") + ylab("umap2") + theme(plot.title = element_blank(), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
  MPP2 <- FeaturePlot(MPP_KO, features=c("Usp22"), slot = "data", combine = FALSE, reduction = "ref.umap", pt.size = 0.3)
  MPP2 <- MPP2[[1]] + xlab("umap1") + ylab("umap2") + theme(plot.title = element_blank(), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
  MPP1[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  MPP2[["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.5
  limits <- MPP1[["scales"]][["scales"]][[3]][["limits"]]
  MPP2[["scales"]][["scales"]][[3]][["limits"]] <- limits
  
  
  TBM_paper_figure <- gridExtra::arrangeGrob(
    grobs = list(gridExtra::arrangeGrob(HSC1, top="Usp22 mRNA in WT", left="LT-HSC"), gridExtra::arrangeGrob(HSC2, top="Usp22 mRNA in KO"), 
                 gridExtra::arrangeGrob(ST1, left="ST-HSC"), ST2, 
                 gridExtra::arrangeGrob(MPP1, left="MPP"), MPP2),
    widths = c(1, 1),
    layout_matrix = rbind(c(1,2),
                          c(3,4),
                          c(5,6))
  )
  
  ggsave(
    path_for_plot,
    plot = TBM_paper_figure,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 15,
    height = 15,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
  
}

# # the function below was used to translate human symbols of cell cycle genes to mouse symbols; to save time and to avoid connection issues with biomart, those genes were saved in files
# convertHumanGeneList <- function(x){
#   require("biomaRt")
#   human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#   mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#   
#   genes = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
#   
#   humanx <- unique(genes[, 2])
#   
#   # Print the first 6 genes found to the screen
#   print(head(humanx))
#   return(humanx)
# }
# g2m.genes <- convertHumanGeneList(cc.genes$g2m.genes)
# s.genes <- convertHumanGeneList(cc.genes$s.genes)
g2m.genes <- readRDS(file = "g2m.genes")
s.genes <- readRDS(file = "s.genes")

run_GSEA <- function(seurat_object, min.pct = 0.1, ont_category="BP", minGSSize = 10, maxGSSize = 800, pvalueCutoff = 0.05, showCategory=50){
  genes_for_GO_analysis <- FindMarkers(seurat_object, group.by = "condition", ident.1 = "KO", ident.2 = "WT", assay = "SCT",
                                       slot = "data", logfc.threshold = 0.0, test.use = "wilcox", min.pct = min.pct)
  named_gene_vector <- genes_for_GO_analysis$avg_log2FC
  names(named_gene_vector) <- rownames(genes_for_GO_analysis)
  named_gene_vector = sort(named_gene_vector, decreasing = TRUE)
  
  library("org.Mm.eg.db", character.only = TRUE)
  geneID <- mapIds(org.Mm.eg.db, keys=names(named_gene_vector),
                   keytype="SYMBOL", column="ENTREZID")
  indices_to_keep <- which(!(is.na(geneID)))
  readable_foldchanges <- named_gene_vector[indices_to_keep]
  names(named_gene_vector) <- geneID
  named_gene_vector <- named_gene_vector[indices_to_keep]
  named_gene_vector <- named_gene_vector[which(!(duplicated(names(named_gene_vector))))] # needs to be done because some mapping produces duplicated names
  
  library(clusterProfiler)
  gse <- gseGO(geneList=named_gene_vector, 
               ont =ont_category, 
               minGSSize = minGSSize, 
               maxGSSize = maxGSSize, 
               pvalueCutoff = pvalueCutoff, 
               verbose = TRUE, 
               OrgDb = org.Mm.eg.db, 
               pAdjustMethod = "BH")
  number_of_significant_GO_terms <- length(which(gse@result$p.adjust < 0.05))
  
  require(DOSE)
  library(enrichplot)
  results <- pairwise_termsim(gse, method = "JC", showCategory=number_of_significant_GO_terms)
  plot <- emapplot(results, showCategory = showCategory,          
                   color = "p.adjust",
                   layout = "nicely",
                   cex_category = 0.2,
                   cex_line = 0.2,
                   min_edge = 0.2,
                   cex_label_category = 0.4)
  return(list(foldchanges=named_gene_vector, readable_foldchanges=readable_foldchanges, results=results, plot=plot))
}

get_GO_term_groups <- function(enrichment_object, number_of_top_GO_terms_to_consider=65, GO_terms_to_consider=NULL, edge_threshold_Jaccard=0.2){
  results_readable <- setReadable(enrichment_object, 'org.Mm.eg.db', 'ENTREZID')
  if (is.null(GO_terms_to_consider)){
    top_GO_terms <- results_readable@result[["Description"]][1:number_of_top_GO_terms_to_consider]
    Jaccard_matrix <- results_readable@termsim[top_GO_terms, top_GO_terms]
  } else {
    Jaccard_matrix <- results_readable@termsim[GO_terms_to_consider, GO_terms_to_consider]
  }
  Jaccard_matrix_threshold_applied <- Jaccard_matrix
  for (index in (1:nrow(Jaccard_matrix_threshold_applied))){
    current_row <- Jaccard_matrix_threshold_applied[index,]
    Jaccard_matrix_threshold_applied[index, which((is.na(current_row)) | (current_row < edge_threshold_Jaccard))] <- 0
  }
  
  list_with_grouped_GO_term <- list()
  for (index in (1:ncol(Jaccard_matrix_threshold_applied))){
    current_colname <- colnames(Jaccard_matrix_threshold_applied)[index]
    direct_neighbors <- rownames(Jaccard_matrix_threshold_applied)[which(Jaccard_matrix_threshold_applied[index,]>0)]
    column_group <- c(current_colname, direct_neighbors)
    list_indices_to_pool <- c()
    if (length(list_with_grouped_GO_term) > 0){
      for (index2 in (1:length(list_with_grouped_GO_term))){
        if (length(intersect(column_group, list_with_grouped_GO_term[[index2]])) > 0){
          list_indices_to_pool <- append(list_indices_to_pool, index2)
        }
      }
    }
    if (length(list_indices_to_pool) > 0){
      new_group <- unique(c(column_group, unlist(list_with_grouped_GO_term[list_indices_to_pool])))
      list_with_grouped_GO_term[list_indices_to_pool] <- NULL
      list_with_grouped_GO_term <- list.append(list_with_grouped_GO_term, new_group)
    } else {
      list_with_grouped_GO_term <- list.append(list_with_grouped_GO_term, column_group)
    }
  }
  return(list_with_grouped_GO_term)
}

get_core_genes_for_each_GO_term <- function(results_readable, GO_terms_of_interest=results_readable@result[["Description"]][1:65]){
  core_genes <- results_readable@result[["core_enrichment"]][which(results_readable@result[["Description"]] %in% GO_terms_of_interest)]
  core_genes_list <- list()
  for (index in (1:length(core_genes))){
    core_genes_list <- list.append(core_genes_list, strsplit(core_genes[index], split="/"))
  }
  names(core_genes_list) <- GO_terms_of_interest
  return(core_genes_list)
}

get_number_of_core_genes_per_GO_term_group <- function(GO_term_groups, core_genes_for_each_GO_term){
  number_of_unique_genes_per_group <- c()
  for (index in (1:length(GO_term_groups))){
    current_GO_terms <- GO_term_groups[[index]]
    number_of_unique_genes_per_group <- append(number_of_unique_genes_per_group, length(unique(unlist(core_genes_for_each_GO_term[current_GO_terms]))))
  }
  return(number_of_unique_genes_per_group)
}

subcluster_GO_group <- function(results_readable, mother_cluster, GO_term_groups_cut_global){
  results_readable_immune <- results_readable
  results_readable_immune@result <- results_readable_immune@result[which(results_readable_immune@result$Description %in% GO_term_groups_cut_global[[mother_cluster]]),]
  JC_sim_immune_terms <- results_readable_immune@termsim[results_readable_immune@result$Description, results_readable_immune@result$Description]
  JC_dist_immune_terms <- (JC_sim_immune_terms*(-1))+1
  diag(JC_dist_immune_terms) <- 0
  gdata::lowerTriangle(JC_dist_immune_terms) <- gdata::upperTriangle(JC_dist_immune_terms, byrow=TRUE)
  JC_dist_immune_terms_dist <- as.dist(JC_dist_immune_terms)
  hclust_ward <- hclust(JC_dist_immune_terms_dist, method = 'ward.D2')
  return(list(hclust_ward=hclust_ward, JC_dist=JC_dist_immune_terms_dist))
}

cut_tree_and_get_list_of_clustered_GO_terms <- function(hclust_ward, number_of_clusters){
  cut_ward <- cutree(hclust_ward, k = number_of_clusters)
  list_of_clustered_immune_GO_terms <- list()
  for (index in (1:number_of_clusters)){
    list_of_clustered_immune_GO_terms <- list.append(list_of_clustered_immune_GO_terms, names(cut_ward)[which(cut_ward==index)])
  }
  return(list_of_clustered_immune_GO_terms)
}

plot_GSEA_heatmap_global <- function(processed_GSEA_results=processed_GSEA_results, log2FC_threshold=1.2, path_to_file="plots_20220216", sample_name = "LTs", figure_height=9, rowlabel_size=6, collabel_size=8, replicate_ID=1){
  
  basic_file_name <- paste(path_to_file, sample_name, sep="/")
  
  core_genes_for_each_GO_term <- processed_GSEA_results[["core_genes_for_each_GO_term"]]
  GO_term_groups_cut_global <- processed_GSEA_results[["GO_term_groups_cut_global"]]
  named_foldchanges_sorted <- processed_GSEA_results[["GSEA_results_sampled"]]$readable_foldchanges
  number_of_clusters <- length(GO_term_groups_cut_global)
  
  all_filtered_GO_terms <- unique(unlist(GO_term_groups_cut_global))
  all_core_genes <- unique(unlist(core_genes_for_each_GO_term[all_filtered_GO_terms]))
  all_core_genes_sorted_FCs <- named_foldchanges_sorted[which(names(named_foldchanges_sorted) %in% all_core_genes)]
  all_core_genes_sorted_FCs_filtered <- all_core_genes_sorted_FCs[which((all_core_genes_sorted_FCs > log2(log2FC_threshold)) | (all_core_genes_sorted_FCs < (-log2(log2FC_threshold))))]
  global_matrix <- matrix(0, ncol = number_of_clusters, nrow = length(all_core_genes_sorted_FCs_filtered))
  for (index in (1:number_of_clusters)){
    global_matrix[,index] <- all_core_genes_sorted_FCs_filtered
  }
  rownames(global_matrix) <- names(all_core_genes_sorted_FCs_filtered)
  colnames(global_matrix) <- names(GO_term_groups_cut_global)
  # set NAs
  for (index in (1:ncol(global_matrix))){
    current_GO_terms <- GO_term_groups_cut_global[[index]]
    current_core_genes <- unique(unlist(core_genes_for_each_GO_term[current_GO_terms]))
    global_matrix[which(!(rownames(global_matrix) %in% current_core_genes)), index] <- NA
  }
  
  mean_values_of_columns <- apply(global_matrix, 2, FUN=function(x)mean(x, na.rm=TRUE))
  genes_per_column <- apply(global_matrix, 2, FUN=function(x)(length(which(!is.na(x)))))
  columns_with_negative_values <- which(mean_values_of_columns<0)
  negative_columns_gene_numbers <- genes_per_column[columns_with_negative_values]
  negative_sorted <- sort(negative_columns_gene_numbers, decreasing = FALSE)
  
  columns_with_positive_values <- which(mean_values_of_columns>0)
  positive_columns_gene_numbers <- genes_per_column[columns_with_positive_values]
  positive_sorted <- sort(positive_columns_gene_numbers, decreasing = TRUE)
  
  global_matrix <- as.matrix(global_matrix[, names(c(positive_sorted, negative_sorted))])
  colnames(global_matrix) <- names(c(positive_sorted, negative_sorted))
  
  color <- colorRampPalette(c("blue", "white", "red"))(100)
  myBreaks <- c(seq(min(global_matrix, na.rm = TRUE), 0, length.out=ceiling(100/2) + 1), 
                seq(max(global_matrix, na.rm = TRUE)/100, max(global_matrix, na.rm = TRUE), length.out=floor(100/2)))
  pheatmap::pheatmap(global_matrix, na_col = "grey", color = color, breaks = myBreaks, cluster_cols = FALSE, cluster_rows = FALSE, 
                     angle_col=270, fontsize_row=rowlabel_size, fontsize_col=collabel_size, fontsize=8, 
                     main = "", 
                     filename = paste(basic_file_name, paste0(paste("global_clusters_core_genes_heat", replicate_ID, sep = "_"), ".pdf"), sep="_"), width = 3, height = figure_height)
}

prepare_MPPs_for_FateID <- function(MPP_WT, MPP_KO, g2m.genes=g2m.genes, s.genes=s.genes){
  set.seed(1)
  print(sample(5))
  DefaultAssay(MPP_WT) <- "SCT"
  MPP_WT <- CellCycleScoring(
    object = MPP_WT,
    g2m.features = g2m.genes,
    s.features = s.genes,
    assay="SCT",
    seed=1
  )
  
  set.seed(1)
  print(sample(5))
  DefaultAssay(MPP_KO) <- "SCT"
  MPP_KO <- CellCycleScoring(
    object = MPP_KO,
    g2m.features = g2m.genes,
    s.features = s.genes,
    assay="SCT",
    seed=1
  )
  
  set.seed(1)
  print(sample(5))
  DefaultAssay(MPP_WT) <- "RNA"
  MPP_WT <- Seurat::SCTransform(MPP_WT, verbose = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3, vars.to.regress = c("S.Score", "G2M.Score"))
  
  set.seed(1)
  print(sample(5))
  DefaultAssay(MPP_KO) <- "RNA"
  MPP_KO <- Seurat::SCTransform(MPP_KO, verbose = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3, vars.to.regress = c("S.Score", "G2M.Score"))
  
  # from both conditions remove cell types which cannot be used to identify lineage markers via FateID ("NK", "T", "preT", "proB", "Baso", "Eo", "DC", "Mono")
  MPP_WT <- MPP_WT[,which(!(MPP_WT$SingleR_labels %in% c("NK", "T", "preT", "B", "proB", "Baso", "Eo", "DC", "Mono")))]
  MPP_KO <- MPP_KO[,which(!(MPP_KO$SingleR_labels %in% c("NK", "T", "preT", "B", "proB", "Baso", "Eo", "DC", "Mono")))]
  
  return(list(MPP_WT=MPP_WT, MPP_KO=MPP_KO))
}

select_lineage_specific_genes_by_FateID <- function(MPP_WT, MPP_KO){
  # select WT features by FateID
  x <- MPP_WT@assays[["SCT"]]@scale.data
  y <- MPP_WT$SingleR_labels
  y[which(y %in% c("GMP"))] <- "1"
  y[which(y %in% c("CLP"))] <- "2"
  y[which(y %in% c("PreCFUE", "MEP"))] <- "3"
  y[which(!(y %in% c("1", "2", "3")))] <- "0"
  y <- as.numeric(y)
  tar <- c(1, 2, 3)
  fb_WT  <- fateBias(x, y, tar, z=NULL, minnr=5, minnrh=11, adapt=TRUE, confidence=0.75, nbfactor=5, use.dist=FALSE, seed=12345, nbtree=NULL)
  # only use important features used in the first 25% of all iterations because classification is less accurate in later iterations
  fb_WT$rfl <- fb_WT$rfl[1:(length(fb_WT$rfl)/4)]
  k1 <- impGenes(fb_WT,"t1",ithr=.02,zthr=4)
  k2 <- impGenes(fb_WT,"t2",ithr=.02,zthr=4)
  k3 <- impGenes(fb_WT,"t3",ithr=.02,zthr=4)
  feature_for_embedding_WT <- unique(c(rownames(k1[["d"]]), rownames(k2[["d"]]), rownames(k3[["d"]])))
  
  # select KO features by FateID
  x <- MPP_KO@assays[["SCT"]]@scale.data
  y <- MPP_KO$SingleR_labels
  y[which(y %in% c("GMP"))] <- "1"
  y[which(y %in% c("CLP"))] <- "2"
  y[which(y %in% c("PreCFUE", "MEP"))] <- "3"
  y[which(!(y %in% c("1", "2", "3")))] <- "0"
  y <- as.numeric(y)
  tar <- c(1, 2, 3)
  fb_KO  <- fateBias(x, y, tar, z=NULL, minnr=5, minnrh=20, adapt=TRUE, confidence=0.75, nbfactor=5, use.dist=FALSE, seed=12345, nbtree=NULL)
  fb_KO$rfl <- fb_KO$rfl[1:(length(fb_KO$rfl)/4)]
  k1 <- impGenes(fb_KO,"t1",ithr=.02,zthr=4)
  k2 <- impGenes(fb_KO,"t2",ithr=.02,zthr=4)
  k3 <- impGenes(fb_KO,"t3",ithr=.02,zthr=4)
  feature_for_embedding_KO <- unique(c(rownames(k1[["d"]]), rownames(k2[["d"]]), rownames(k3[["d"]])))
  
  # select marker genes both groups have in common to generate a common embedding
  common_features <- feature_for_embedding_WT[which(feature_for_embedding_WT %in% feature_for_embedding_KO)]
  return(common_features)
}

compute_diffusion_map <- function(MPP_combined, common_features, replicate_ID=1){
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
  
  # compute diffusion map based on selected features
  dmap <- destiny::DiffusionMap(reducedDim(MPP_combined.sce,'FILTERED'))
  reducedDim(MPP_combined.sce, 'DiffusionMap') <- dmap@eigenvectors
  
  # visualize initial diffusion map
  dmap_plot_unfiltered <- scatter_3d(MPP_combined.sce, dim.red.name = 'DiffusionMap', color.column="SingleR_labels", marker_size = 20, scene = '')
  
  # some extreme outliers in replicate 1 need to be removed to increase the resolution of lineage specification
  if (is.null(replicate_ID) == FALSE){
    if (replicate_ID==1){
      # remove extreme cells (some MEPs and PreCFUEs) and re-compute diffusion map
      cells_to_keep <- rownames(reducedDim(MPP_combined.sce, 'DiffusionMap'))[which(reducedDim(MPP_combined.sce, 'DiffusionMap')[,3]>(-0.1))]
      MPP_combined.sce <- MPP_combined.sce[,cells_to_keep]
      MPP_combined <- MPP_combined[,cells_to_keep]
      dmap <- destiny::DiffusionMap(reducedDim(MPP_combined.sce,'FILTERED'))
      reducedDim(MPP_combined.sce, 'DiffusionMap') <- dmap@eigenvectors
      dmap_plot_filtered <- scatter_3d(MPP_combined.sce, dim.red.name = 'DiffusionMap', color.column="SingleR_labels", marker_size = 20, scene = '')
    } else {
      dmap_plot_filtered <- dmap_plot_unfiltered
    }
  } else {
    dmap_plot_filtered <- dmap_plot_unfiltered
  }
  return(list(MPP_combined=MPP_combined, MPP_combined.sce=MPP_combined.sce, dmap_plot_unfiltered=dmap_plot_unfiltered, dmap_plot_filtered=dmap_plot_filtered))
}

assign_cells_to_branches_using_slingshot <- function(MPP_combined, MPP_combined.sce, replicate_ID=1){
  # create clusters (parameters for clustering were selected in a way to ensure meaningful curve fitting by the slingshot software)
  MPP_combined.seurat <- as.Seurat(MPP_combined.sce)
  if (replicate_ID==2){
    MPP_combined.seurat <- cluster_cells(MPP_combined.seurat, reduction = "DiffusionMap", ndims=3, knn=50, create_neighbor_object=FALSE, resolution=0.3)
  } else {
    MPP_combined.seurat <- cluster_cells(MPP_combined.seurat, reduction = "DiffusionMap", ndims=3, knn=50, create_neighbor_object=FALSE, resolution=0.1)
  }
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
  
  if (replicate_ID==2){
    # summarize all clusters that are not endpoint clusters
    endpoint_clusters <- c(clusters_corresponding_to_SingleR_labels["PreCFUE"], clusters_corresponding_to_SingleR_labels["CLP"], clusters_corresponding_to_SingleR_labels["GMP"])
    all_clusters <- as.character(unique(MPP_combined.seurat$seurat_clusters))
    central_clusters <- all_clusters[which(!(all_clusters %in% endpoint_clusters))]
    MPP_combined.sce$Cluster1[which(MPP_combined.sce$Cluster1 %in% central_clusters)] <- "central"
  }
  
  # run slingshot
  MPP_combined.sce@int_colData@listData[["reducedDims"]]@listData[["DiffusionMap_subset"]] <- MPP_combined.sce@int_colData@listData[["reducedDims"]]@listData[["DiffusionMap"]][,1:3]
  MPP_combined.sce <- slingshot(MPP_combined.sce, clusterLabels = 'Cluster1', reducedDim = 'DiffusionMap_subset', start.clus=c(clusters_corresponding_to_SingleR_labels["PreCFUE"]), end.clus=c(clusters_corresponding_to_SingleR_labels["CLP"], clusters_corresponding_to_SingleR_labels["GMP"]), allow.breaks=TRUE)
  
  # assign curves to lineages
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
  return(list(MPP_combined=MPP_combined, MPP_combined.sce=MPP_combined.sce, clusters_plot=clusters_plot))
}

compute_per_cell_scores_of_GMP_LMPP_priming_and_plot_distribution <- function(object_for_signature=MPP_combined, 
                                                                              objects_for_scoring=list(MPP_combined.corrected_cell_numbers=MPP_combined.corrected_cell_numbers, 
                                                                                                       HSC_combined=HSC_combined),
                                                                              path_for_plot= "plots_20220216/",
                                                                              replicate_ID=NULL){
  if (replicate_ID==1){
    if ((file.exists("GMP_signature")==FALSE) | (file.exists("CLP_signature")==FALSE)){
      MPP_combined <- object_for_signature
      MPP_combined.WT <- MPP_combined[,which(MPP_combined$condition=="WT")]
      WT_GMP_markers1 <- FindMarkers(MPP_combined.WT, group.by = "branches", assay = "SCT", slot = "data",ident.1 = "GMP", ident.2 = "CLP", test.use = "wilcox", min.pct = 0.05, logfc.threshold = log2(1.25))
      WT_GMP_markers2 <- FindMarkers(MPP_combined.WT, group.by = "branches", assay = "SCT", slot = "data", ident.1 = "GMP", ident.2 = "undefined", test.use = "wilcox", min.pct = 0.05, logfc.threshold = log2(1.25))
      WT_GMP_markers3 <- FindMarkers(MPP_combined.WT, group.by = "branches", assay = "SCT", slot = "data", ident.1 = "GMP", ident.2 = "MEP", test.use = "wilcox", min.pct = 0.05, logfc.threshold = log2(1.25))
      WT_GMP_markers1_positve <- rownames(WT_GMP_markers1[which((WT_GMP_markers1$avg_log2FC>log2(1.25)) & (WT_GMP_markers1$p_val_adj<0.05)),])
      WT_GMP_markers2_positve <- rownames(WT_GMP_markers2[which((WT_GMP_markers2$avg_log2FC>log2(1.25)) & (WT_GMP_markers2$p_val_adj<0.05)),])
      WT_GMP_markers3_positve <- rownames(WT_GMP_markers3[which((WT_GMP_markers3$avg_log2FC>log2(1.25)) & (WT_GMP_markers3$p_val_adj<0.05)),])
      GMP_signature <- intersect(intersect(WT_GMP_markers1_positve, WT_GMP_markers2_positve), WT_GMP_markers3_positve)
      saveRDS(object = GMP_signature, file = "GMP_signature")
      
      WT_CLP_markers1 <- FindMarkers(MPP_combined.WT, group.by = "branches", assay = "SCT", slot = "data", ident.1 = "CLP", ident.2 = "GMP", test.use = "wilcox", min.pct = 0.05, logfc.threshold = log2(1.25))
      WT_CLP_markers2 <- FindMarkers(MPP_combined.WT, group.by = "branches", assay = "SCT", slot = "data", ident.1 = "CLP", ident.2 = "undefined", test.use = "wilcox", min.pct = 0.05, logfc.threshold = log2(1.25))
      WT_CLP_markers3 <- FindMarkers(MPP_combined.WT, group.by = "branches", assay = "SCT", slot = "data", ident.1 = "CLP", ident.2 = "MEP", test.use = "wilcox", min.pct = 0.05, logfc.threshold = log2(1.25))
      WT_CLP_markers1_positve <- rownames(WT_CLP_markers1[which((WT_CLP_markers1$avg_log2FC>log2(1.25)) & (WT_CLP_markers1$p_val_adj<0.05)),])
      WT_CLP_markers2_positve <- rownames(WT_CLP_markers2[which((WT_CLP_markers2$avg_log2FC>log2(1.25)) & (WT_CLP_markers2$p_val_adj<0.05)),])
      WT_CLP_markers3_positve <- rownames(WT_CLP_markers3[which((WT_CLP_markers3$avg_log2FC>log2(1.25)) & (WT_CLP_markers3$p_val_adj<0.05)),])
      CLP_signature <- intersect(intersect(WT_CLP_markers1_positve, WT_CLP_markers2_positve), WT_CLP_markers3_positve)
      saveRDS(object = CLP_signature, file = "CLP_signature")
    } else {
      GMP_signature <- readRDS(file = "GMP_signature")
      CLP_signature <- readRDS(file = "CLP_signature")
    }
  } else {
    GMP_signature <- readRDS(file = "GMP_signature")
    CLP_signature <- readRDS(file = "CLP_signature")
  }
  
  MPP_combined.corrected_cell_numbers <- objects_for_scoring[["MPP_combined.corrected_cell_numbers"]]
  
  HSC_combined <- objects_for_scoring[["HSC_combined"]]
  
  
  # compute signature scores for MPPs and rename braches and set levels of categorical variables for plotting
  MPP_combined.corrected_cell_numbers <- AddModuleScore(
    object=MPP_combined.corrected_cell_numbers,
    features=list(GMP_signature, CLP_signature),
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    k = FALSE,
    assay = "SCT",
    name = "module_score",
    seed = 1,
    search = FALSE
  )
  MPP_combined.corrected_cell_numbers@meta.data[["condition_violin_plot"]] <- MPP_combined.corrected_cell_numbers$condition
  MPP_combined.corrected_cell_numbers$condition_violin_plot <- factor(MPP_combined.corrected_cell_numbers$condition_violin_plot, levels=c("WT", "KO"))
  MPP_combined.corrected_cell_numbers@meta.data[["branches_violin_plot"]] <- MPP_combined.corrected_cell_numbers$branches
  MPP_combined.corrected_cell_numbers$branches_violin_plot[which(MPP_combined.corrected_cell_numbers$branches_violin_plot=="CLP")] <- "LMPPs"
  MPP_combined.corrected_cell_numbers$branches_violin_plot[which(MPP_combined.corrected_cell_numbers$branches_violin_plot=="GMP")] <- "GMP-biased MPPs"
  MPP_combined.corrected_cell_numbers$branches_violin_plot[which(MPP_combined.corrected_cell_numbers$branches_violin_plot=="MEP")] <- "MEP-biased MPPs"
  MPP_combined.corrected_cell_numbers$branches_violin_plot[which(MPP_combined.corrected_cell_numbers$branches_violin_plot=="undefined")] <- "unbiased MPPs"
  
  
  # compute signature scores for HSCs and rename branches and set levels of categorical variables for plotting
  HSC_combined <- AddModuleScore(
    object=HSC_combined,
    features=list(GMP_signature, CLP_signature),
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    k = FALSE,
    assay = "SCT",
    name = "module_score",
    seed = 1,
    search = FALSE
  )
  HSC_combined$condition <- factor(HSC_combined$condition, levels = c("WT", "KO"))
  HSC_combined$sample[which(HSC_combined$sample=="LTHSC")] <- "LT-HSCs"
  HSC_combined$sample[which(HSC_combined$sample=="STHSC")] <- "ST-HSCs"
  HSC_combined$sample <- factor(HSC_combined$sample, levels = c("LT-HSCs", "ST-HSCs"))
  
  
  # create violin plots
  Vln_Plot_MPPs_GMPsignature <- VlnPlot(
    object=MPP_combined.corrected_cell_numbers,
    features="module_score1",
    cols = c("grey", "red"),
    pt.size = NULL,
    idents = NULL,
    sort = FALSE,
    assay = NULL,
    group.by = "branches_violin_plot",
    split.by = "condition_violin_plot",
    adjust = 1,
    y.max = NULL,
    same.y.lims = FALSE,
    log = FALSE,
    ncol = NULL,
    slot = "data",
    split.plot = FALSE,
    stack = FALSE,
    combine = FALSE,
    fill.by = "feature",
    flip = FALSE
  )[[1]]
  Vln_Plot_MPPs_GMPsignature_prep <- Vln_Plot_MPPs_GMPsignature + ggtitle("") + theme(axis.text.x = element_text(angle = 270, size=8), axis.title.x = element_blank(), plot.title = element_text(size=10)) +
    stat_compare_means(aes(group = split), method="wilcox.test", paired=FALSE, label = "p.signif", hide.ns=FALSE)
  
  
  Vln_Plot_MPPs_LMPPsignature <- VlnPlot(
    object=MPP_combined.corrected_cell_numbers,
    features="module_score2",
    cols = c("grey", "red"),
    pt.size = NULL,
    idents = NULL,
    sort = FALSE,
    assay = NULL,
    group.by = "branches_violin_plot",
    split.by = "condition_violin_plot",
    adjust = 1,
    y.max = NULL,
    same.y.lims = FALSE,
    log = FALSE,
    ncol = NULL,
    slot = "data",
    split.plot = FALSE,
    stack = FALSE,
    combine = TRUE,
    fill.by = "feature",
    flip = FALSE
  )[[1]]
  
  Vln_Plot_MPPs_LMPPsignature_prep <- Vln_Plot_MPPs_LMPPsignature + ggtitle("") + theme(axis.text.x = element_text(angle = 270, size=8), axis.title.x = element_blank(), plot.title = element_text(size=10)) +
    stat_compare_means(aes(group = split), method="wilcox.test", paired=FALSE, label = "p.signif", hide.ns=FALSE)
  
  
  Vln_Plot_HSCs_GMPsignature <- VlnPlot(
    object=HSC_combined,
    features="module_score1",
    cols = c("grey", "red"),
    pt.size = NULL,
    idents = NULL,
    sort = FALSE,
    assay = NULL,
    group.by = "sample",
    split.by = "condition",
    adjust = 1,
    y.max = NULL,
    same.y.lims = FALSE,
    log = FALSE,
    ncol = NULL,
    slot = "data",
    split.plot = FALSE,
    stack = FALSE,
    combine = TRUE,
    fill.by = "feature",
    flip = FALSE
  )[[1]]
  
  Vln_Plot_HSCs_GMPsignature_prep <- Vln_Plot_HSCs_GMPsignature + ggtitle("") + theme(axis.text.x = element_text(angle = 270, size=8), axis.title.x = element_blank(), plot.title = element_text(size=10)) + 
    stat_compare_means(aes(group = split), method="wilcox.test", paired=FALSE, label = "p.signif", hide.ns=FALSE)
  
  
  Vln_Plot_HSCs_LMPPsignature <- VlnPlot(
    object=HSC_combined,
    features="module_score2",
    cols = c("grey", "red"),
    pt.size = NULL,
    idents = NULL,
    sort = FALSE,
    assay = NULL,
    group.by = "sample",
    split.by = "condition",
    adjust = 1,
    y.max = NULL,
    same.y.lims = FALSE,
    log = FALSE,
    ncol = NULL,
    slot = "data",
    split.plot = FALSE,
    stack = FALSE,
    combine = TRUE,
    fill.by = "feature",
    flip = FALSE
  )[[1]]
  
  Vln_Plot_HSCs_LMPPsignature_prep <- Vln_Plot_HSCs_LMPPsignature + ggtitle("") + theme(axis.text.x = element_text(angle = 270, size=8), axis.title.x = element_blank(), plot.title = element_text(size=10)) +
    stat_compare_means(aes(group = split), method="wilcox.test", paired=FALSE, label = "p.signif", hide.ns=FALSE)
  
  
  
  # group plots by signature Vln_Plot_MPPs_GMPsignature, Vln_Plot_MPPs_LMPPsignature, Vln_Plot_HSCs_GMPsignature, Vln_Plot_HSCs_LMPPsignature
  combined_plot_GMP <- Vln_Plot_MPPs_GMPsignature_prep + Vln_Plot_HSCs_GMPsignature_prep + plot_layout(guides = 'collect', ncol = 1, nrow = 2)
  combined_plot_LMPP <- Vln_Plot_MPPs_LMPPsignature_prep + Vln_Plot_HSCs_LMPPsignature_prep + plot_layout(guides = 'collect', ncol = 1, nrow = 2)
  
  ggsave(
    paste(path_for_plot, "Vln_Plot_panel_GMP.pdf", sep=""),
    plot = combined_plot_GMP,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 10,
    height = 20,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
  
  ggsave(
    paste(path_for_plot, "Vln_Plot_panel_LMPP.pdf", sep=""),
    plot = combined_plot_LMPP,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 10,
    height = 20,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
  
  return(list(GMP_signature=GMP_signature, CLP_signature=CLP_signature))
}

assign_cells_to_cell_cycle_phases_based_on_gene_lists <- function(object_to_annotate=MPP_combined,
                                                                  G2M_thres_in_neighborhood=1,
                                                                  neighborhood_size=200,
                                                                  g2m.genes=g2m.genes,
                                                                  s.genes=s.genes){
  if ((is.null(g2m.genes)) | (is.null(s.genes))){
    stop("Lists with cell cycle genes could not be loaded!")
  }
  if (!("SCT" %in% names(object_to_annotate@assays))){
    stop("SCT assay required for cell cycle classification!")
  }
  DefaultAssay(object_to_annotate) <- "SCT"
  object_to_annotate <- CellCycleScoring(
    object = object_to_annotate,
    g2m.features = g2m.genes,
    s.features = s.genes)
  object_to_annotate <- Seurat::RunUMAP(object_to_annotate, features=c(s.genes, g2m.genes), assay="SCT", slot="data", return.model = TRUE, n.neighbors = 30, min.dist = 0.7, metric="cosine", reduction.name="umap_cycle")#25 0.7
  DimPlot(object_to_annotate, reduction = "umap_cycle", group.by = "Phase")
  object_to_annotate <- FindNeighbors(object_to_annotate, reduction = "umap_cycle", dims=1:2, k.param = neighborhood_size, return.neighbor = TRUE)
  neighbor_matrix <- object_to_annotate@neighbors[["SCT.nn"]]@nn.idx
  G1S_indices <- which(object_to_annotate$Phase %in% c("G1", "S"))
  G1_indices <- which(object_to_annotate$Phase %in% c("G1"))
  S_indices <- which(object_to_annotate$Phase %in% c("S"))
  G2M_indices <- which(object_to_annotate$Phase %in% c("G2M"))
  neighbor_matrix_G1S <- neighbor_matrix[G1S_indices,]
  new_G1S_members <- c()
  for (index in (1:nrow(neighbor_matrix_G1S))){
    current_row <- neighbor_matrix_G1S[index,]
    if (length(which(current_row %in% G2M_indices)) < G2M_thres_in_neighborhood){
      number_of_cells_in_G1 <- length(which(current_row %in% G1_indices))
      number_of_cells_in_S <- length(which(current_row %in% S_indices))
      proportion_of_smaller_fraction <- min(c(number_of_cells_in_G1, number_of_cells_in_S)) / length(current_row)
      if (proportion_of_smaller_fraction > 0.05){
        new_G1S_members <- append(new_G1S_members, current_row[1])
      }
    }
  }
  object_to_annotate@meta.data[["new_Phase"]] <- object_to_annotate$Phase
  object_to_annotate$new_Phase[new_G1S_members] <- "G1S"
  umap_cell_cycle <- DimPlot(object_to_annotate, reduction = "umap_cycle", group.by = "new_Phase")
  return(list(umap_cell_cycle=umap_cell_cycle, annotated_object=object_to_annotate))
}

plot_cell_cycle_changes_as_heatmap <- function(annotated_object=MPP_combined, path_to_plots="plots_20220216"){
  matrix_diff_cycle <- matrix(0, nrow=4, ncol = length(unique(annotated_object$branches)))
  colnames(matrix_diff_cycle) <- c("undefined", "CLP", "GMP", "MEP")
  rownames(matrix_diff_cycle) <- c("G1", "G1S", "S", "G2M")
  for (index in (1:length(unique(annotated_object$branches)))){
    current_cluster <- unique(annotated_object$branches)[index]
    number_of_WT_cells_in_current_cluster <- length(which((annotated_object$condition=="WT") & (annotated_object$branches==current_cluster)))
    number_of_KO_cells_in_current_cluster <- length(which((annotated_object$condition=="KO") & (annotated_object$branches==current_cluster)))
    if (min(c(number_of_WT_cells_in_current_cluster, number_of_KO_cells_in_current_cluster)) > 100){
      proportion_G1_WT <- length(which((annotated_object$condition=="WT") & (annotated_object$branches==current_cluster) & (annotated_object$new_Phase=="G1")))/number_of_WT_cells_in_current_cluster
      proportion_G1_KO <- length(which((annotated_object$condition=="KO") & (annotated_object$branches==current_cluster) & (annotated_object$new_Phase=="G1")))/number_of_KO_cells_in_current_cluster
      matrix_diff_cycle["G1", current_cluster] <- proportion_G1_KO - proportion_G1_WT
      
      proportion_G1S_WT <- length(which((annotated_object$condition=="WT") & (annotated_object$branches==current_cluster) & (annotated_object$new_Phase=="G1S")))/number_of_WT_cells_in_current_cluster
      proportion_G1S_KO <- length(which((annotated_object$condition=="KO") & (annotated_object$branches==current_cluster) & (annotated_object$new_Phase=="G1S")))/number_of_KO_cells_in_current_cluster
      matrix_diff_cycle["G1S", current_cluster] <- proportion_G1S_KO - proportion_G1S_WT
      
      proportion_S_WT <- length(which((annotated_object$condition=="WT") & (annotated_object$branches==current_cluster) & (annotated_object$new_Phase=="S")))/number_of_WT_cells_in_current_cluster
      proportion_S_KO <- length(which((annotated_object$condition=="KO") & (annotated_object$branches==current_cluster) & (annotated_object$new_Phase=="S")))/number_of_KO_cells_in_current_cluster
      matrix_diff_cycle["S", current_cluster] <- proportion_S_KO - proportion_S_WT
      
      proportion_G2M_WT <- length(which((annotated_object$condition=="WT") & (annotated_object$branches==current_cluster) & (annotated_object$new_Phase=="G2M")))/number_of_WT_cells_in_current_cluster
      proportion_G2M_KO <- length(which((annotated_object$condition=="KO") & (annotated_object$branches==current_cluster) & (annotated_object$new_Phase=="G2M")))/number_of_KO_cells_in_current_cluster
      matrix_diff_cycle["G2M", current_cluster] <- proportion_G2M_KO - proportion_G2M_WT
    }
  }
  colnames(matrix_diff_cycle) <- c("unbiased MPPs", "LMPPs", "GMP-biased MPPs", "MEP-biased MPPs")
  
  color <- colorRampPalette(c("blue", "white", "red"))(100)
  myBreaks <- c(seq(min(matrix_diff_cycle), 0, length.out=ceiling(100/2) + 1), 
                seq(max(matrix_diff_cycle)/100, max(matrix_diff_cycle), length.out=floor(100/2)))
  pheatmap(matrix_diff_cycle, color = color, breaks = myBreaks, cluster_rows=FALSE, cluster_cols=FALSE, angle_col=270,
           fontsize_row=8, fontsize_col=8, fontsize=8, 
           main = "Differences of\ncell cycle\nphase proportions", 
           filename = paste(path_to_plots, "cell_cycle_heatmap_unsampled.pdf", sep = "/"), width = 2.5, height = 3)
  dev.off()
}

plot_cell_cycle_changes_in_barplot_MPPs <- function(annotated_object=MPP_combined, path_to_plots=folder_with_plots, width = 5, height = 8){
  matrix_diff_cycle <- data.frame(branch=c(rep("undefined", 8),
                                           rep("CLP", 8),
                                           rep("GMP", 8),
                                           rep("MEP", 8)),
                                  condition=c(rep(c(rep("WT", 4), rep("KO", 4)), 4)),
                                  phase=c(rep(c("G1", "G1S", "S", "G2M"), 8)),
                                  percentage=rep(0, 4*2*4))
  
  for (index in (1:nrow(matrix_diff_cycle))){
    branch <- matrix_diff_cycle[index, "branch"]
    condition <- matrix_diff_cycle[index, "condition"]
    phase <- matrix_diff_cycle[index, "phase"]
    number_of_cells_in_current_condition_and_cluster <- length(which((annotated_object$condition==condition) & (annotated_object$branches==branch)))
    matrix_diff_cycle[index, "percentage"] <- (length(which((annotated_object$condition==condition) & (annotated_object$branches==branch) & (annotated_object$new_Phase==phase)))/number_of_cells_in_current_condition_and_cluster)*100
  }
  
  matrix_diff_cycle[which(matrix_diff_cycle[,"branch"]=="undefined"), "branch"] <- "unbiased MPPs"
  matrix_diff_cycle[which(matrix_diff_cycle[,"branch"]=="CLP"), "branch"] <- "LMPPs"
  matrix_diff_cycle[which(matrix_diff_cycle[,"branch"]=="GMP"), "branch"] <- "GMP-biased MPPs"
  matrix_diff_cycle[which(matrix_diff_cycle[,"branch"]=="MEP"), "branch"] <- "MEP-biased MPPs"
  
  matrix_diff_cycle[,"branch"] <- factor(matrix_diff_cycle[,"branch"], levels = c("unbiased MPPs", "LMPPs", "GMP-biased MPPs", "MEP-biased MPPs"))
  matrix_diff_cycle[,"condition"] <- factor(matrix_diff_cycle[,"condition"], levels = c("WT", "KO"))
  matrix_diff_cycle[,"phase"] <- factor(matrix_diff_cycle[,"phase"], levels = c("G1", "G1S", "S", "G2M"))
  
  
  bar_plot <- ggplot(matrix_diff_cycle,
                     aes(x = condition,
                         y = percentage,
                         fill = phase)) + 
    geom_bar(stat = "identity",
             position = "stack") +
    facet_grid(~ branch, switch=NULL) + 
    ylab("%") +
    theme(strip.text.x.top = element_text(angle=270), axis.title.y = element_text(angle=0), panel.grid = element_blank())
  
  ggsave(
    paste(path_to_plots, "MPPs_cell_cycle_barplot.pdf", sep="/"),
    plot = bar_plot,
    device = NULL,
    path = NULL,
    scale = 1,
    width = width,
    height = height,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
}

plot_cell_cycle_changes_in_barplot_LTs <- function(annotated_object=LTHSC_combined, path_to_plots=folder_with_plots, width = 5, height = 8){
  matrix_diff_cycle <- data.frame(condition=c(c(rep("WT", 4), rep("KO", 4))),
                                  phase=c(rep(c("G1", "G1S", "S", "G2M"), 2)),
                                  percentage=rep(0, 4*2))
  
  for (index in (1:nrow(matrix_diff_cycle))){
    condition <- matrix_diff_cycle[index, "condition"]
    phase <- matrix_diff_cycle[index, "phase"]
    number_of_cells_in_current_condition_and_cluster <- length(which(annotated_object$condition==condition))
    matrix_diff_cycle[index, "percentage"] <- (length(which((annotated_object$condition==condition) & (annotated_object$new_Phase==phase)))/number_of_cells_in_current_condition_and_cluster)*100
  }
  
  matrix_diff_cycle[,"condition"] <- factor(matrix_diff_cycle[,"condition"], levels = c("WT", "KO"))
  matrix_diff_cycle[,"phase"] <- factor(matrix_diff_cycle[,"phase"], levels = c("G1", "G1S", "S", "G2M"))
  
  bar_plot <- ggplot(matrix_diff_cycle,
                     aes(x = condition,
                         y = percentage,
                         fill = phase)) + 
    geom_bar(stat = "identity",
             position = "stack") +
    ylab("%") +
    theme(axis.title.y = element_text(angle=0), panel.grid = element_blank())
  
  ggsave(
    paste(path_to_plots, "LTs_cell_cycle_barplot.pdf", sep="/"),
    plot = bar_plot,
    device = NULL,
    path = NULL,
    scale = 1,
    width = width,
    height = height,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
}

create_barplot_with_percentage_of_cycling_LTs <- function(LTHSC_combined=LTHSC_combined, path_to_plots="plots_20220216"){
  percentage_non_cycling_WT <- length(which((LTHSC_combined$condition == "WT") & (LTHSC_combined$new_Phase %in% c("G1", "G1S")))) / length(which(LTHSC_combined$condition == "WT"))
  percentage_non_cycling_KO <- length(which((LTHSC_combined$condition == "KO") & (LTHSC_combined$new_Phase %in% c("G1", "G1S")))) / length(which(LTHSC_combined$condition == "KO"))
  df_plotting <- data.frame(condition=c("WT", "KO"),
                            percentage=c(percentage_non_cycling_WT, percentage_non_cycling_KO))
  df_plotting$condition <- factor(df_plotting$condition, levels = c("WT", "KO"))
  df_plotting$percentage <- df_plotting$percentage*100
  bar_plot <- ggplot(data=df_plotting, aes(x=condition, y=percentage, color=condition)) +
    geom_bar(stat="identity", fill="white") + theme_classic() + scale_color_manual(values=c("black", "red"))
  ggsave(
    paste(path_to_plots, "cycling_LTs_barplot.pdf", sep="/"),
    plot = bar_plot,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 5,
    height = 8,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
}

correct_MPP_embedding_for_differential_ditribution_of_conditions <- function(MPP_combined=MPP_combined, MPP_combined.sce=MPP_combined.sce, numbers_of_direct_neighbors_to_consider=1){
  # add diffusion map coordinates to seurat object
  temp_seurat <- as.Seurat(MPP_combined.sce)
  MPP_combined@reductions[["DiffusionMap_subset"]] <- temp_seurat@reductions[["DiffusionMap_subset"]]
  rm(temp_seurat)
  # compute matrix with nearest neighbor indices based on first three components of diffusion map
  MPP_combined <- FindNeighbors(MPP_combined, reduction = "DiffusionMap_subset", dims=1:3, k.param = 50, return.neighbor = TRUE)
  neighbor_matrix <- MPP_combined@neighbors[["RNA.nn"]]@nn.idx
  # create table with WT cells whose first neighbor is a KO cell; one WT-KO pair per row; each KO cell can be included multiple times
  KO_indices <- which(MPP_combined$condition=="KO")
  nearest_KO_neighbors <- c()
  for (index in (1:nrow(neighbor_matrix))){
    current_neighbors <- neighbor_matrix[index, (2 : (numbers_of_direct_neighbors_to_consider + 1))]
    nearest_KO_neighbors <- append(nearest_KO_neighbors, current_neighbors[which(current_neighbors %in% KO_indices)][1])
  }
  table_for_cell_selection <- data.frame(central_cell=c(1:ncol(MPP_combined)), neighbor=nearest_KO_neighbors)
  rownames(table_for_cell_selection) <- colnames(MPP_combined)
  table_for_cell_selection <- table_for_cell_selection[which(!(is.na(table_for_cell_selection[,"neighbor"]))),]
  table_for_cell_selection <- table_for_cell_selection[table_for_cell_selection[,"central_cell"] %in% (which(MPP_combined$condition=="WT")),]
  # subset MPP embedding to obtain a new embedding including all WT-KO pairs identified above; some KO are sampled multiple times during subsetting
  MPP_combined.corrected_cell_numbers <- MPP_combined[, c(table_for_cell_selection[,"central_cell"], table_for_cell_selection[,"neighbor"])]
  MPP_combined.sce.corrected_cell_numbers <- MPP_combined.sce[, c(table_for_cell_selection[,"central_cell"], table_for_cell_selection[,"neighbor"])]
  # rename cells that were sampled multiple times to make cell names unique again
  copy_numbers <- rep(0, length(unique(colnames(MPP_combined.sce.corrected_cell_numbers))))
  names(copy_numbers) <- unique(colnames(MPP_combined.sce.corrected_cell_numbers))
  cell_specifications <- rep(0, ncol(MPP_combined.sce.corrected_cell_numbers))
  for (index in (1:ncol(MPP_combined.sce.corrected_cell_numbers))){
    current_barcode <- colnames(MPP_combined.sce.corrected_cell_numbers)[index]
    copy_numbers[current_barcode] <- copy_numbers[current_barcode]+1
    cell_specifications[index] <- copy_numbers[current_barcode]
  }
  colnames(MPP_combined.sce.corrected_cell_numbers) <- paste(colnames(MPP_combined.sce.corrected_cell_numbers), as.character(cell_specifications), sep="_")
  MPP_combined.corrected_cell_numbers <- RenameCells(MPP_combined.corrected_cell_numbers, new.names = paste(colnames(MPP_combined.corrected_cell_numbers), as.character(cell_specifications), sep="_"))
  scatterplot_sampled_object <- scatter_3d(MPP_combined.sce.corrected_cell_numbers, dim.red.name = 'DiffusionMap', color.column="condition", marker_size = 20, scene = '')
  
  # calculate how often individual KO cells were sampled:
  max_sampling_number <- max(table(table_for_cell_selection[,"neighbor"]))
  percentage_of_KO_cells_sampled_once <- length(which(table(table_for_cell_selection[,"neighbor"])==1)) / length(table(table_for_cell_selection[,"neighbor"]))
  percentage_of_KO_cells_sampled_twice <- length(which(table(table_for_cell_selection[,"neighbor"])==2)) / length(table(table_for_cell_selection[,"neighbor"]))
  
  return(list(scatterplot_sampled_object=scatterplot_sampled_object, 
              MPP_combined.corrected_cell_numbers=MPP_combined.corrected_cell_numbers,
              MPP_combined.sce.corrected_cell_numbers=MPP_combined.sce.corrected_cell_numbers,
              max_sampling_number=max_sampling_number,
              percentage_of_KO_cells_sampled_once=percentage_of_KO_cells_sampled_once,
              percentage_of_KO_cells_sampled_twice=percentage_of_KO_cells_sampled_twice))
}

create_plots_showing_MPP_embedding <- function(MPP_combined.sce, path_to_plot_parameters="3d_scatter_parameters.rds", path_to_plots="plots_20220216"){
  
  # sample WT and KO embedding to equal cell numbers and read plotting parameters
  number_of_cells_to_plot_per_condition <- min(c(length(which(MPP_combined.sce$condition=="WT")), length(which(MPP_combined.sce$condition=="KO"))))
  WT_cells_to_plot <- sample(which(MPP_combined.sce$condition=="WT"), size = number_of_cells_to_plot_per_condition, replace = FALSE)
  KO_cells_to_plot <- sample(which(MPP_combined.sce$condition=="KO"), size = number_of_cells_to_plot_per_condition, replace = FALSE)
  MPP_combined.sce_sampled <- MPP_combined.sce[, c(WT_cells_to_plot, KO_cells_to_plot)]
  #sds <- SlingshotDataSet(MPP_combined.sce_sampled)
  save <- readRDS(file=path_to_plot_parameters)
  
  # plot WT cells
  plotcol <- MPP_combined.sce_sampled[,which(MPP_combined.sce_sampled$condition=="WT")]$branches
  plotcol[which(plotcol=="CLP")] <- "orange"
  plotcol[which(plotcol=="GMP")] <- "blue"
  plotcol[which(plotcol=="MEP")] <- "red"
  plotcol[which(plotcol=="undefined")] <- "grey"
  plot3d(reducedDims(MPP_combined.sce_sampled[,which(MPP_combined.sce_sampled$condition=="WT")])$DiffusionMap_subset, col = plotcol, size = 4, alpha=0.5)
  #plot3d.SlingshotDataSet(sds, lwd = 2, add = TRUE)
  par3d(save)
  snapshot3d(filename = paste(path_to_plots, "WT_diffusion_map.png", sep="/"), 
             fmt = "png", width = 5000, height = 5000)
  img <- image_read(paste(path_to_plots, "WT_diffusion_map.png", sep="/"))
  image_write(img, path = paste(path_to_plots, "WT_diffusion_map.pdf", sep="/"), format = "pdf")
  rgl.close()
  
  # plot KO cells
  plotcol <- MPP_combined.sce_sampled[,which(MPP_combined.sce_sampled$condition=="KO")]$branches
  plotcol[which(plotcol=="CLP")] <- "orange"
  plotcol[which(plotcol=="GMP")] <- "blue"
  plotcol[which(plotcol=="MEP")] <- "red"
  plotcol[which(plotcol=="undefined")] <- "grey"
  plot3d(reducedDims(MPP_combined.sce_sampled[,which(MPP_combined.sce_sampled$condition=="KO")])$DiffusionMap_subset, col = plotcol, size = 4, alpha=0.5)
  #plot3d.SlingshotDataSet(sds, lwd = 2, add = TRUE)
  par3d(save)
  snapshot3d(filename = paste(path_to_plots, "KO_diffusion_map.png", sep="/"), 
             fmt = "png", width = 5000, height = 5000)
  img <- image_read(paste(path_to_plots, "KO_diffusion_map.png", sep="/"))
  image_write(img, path = paste(path_to_plots, "KO_diffusion_map.pdf", sep="/"), format = "pdf")
  rgl.close()
}

create_plots_showing_MPP_embedding_with_curves <- function(MPP_combined.sce, path_to_plot_parameters="3d_scatter_parameters.rds", path_to_plots="plots_20220216"){
  
  # sample WT and KO embedding to equal cell numbers and read plotting parameters
  number_of_cells_to_plot_per_condition <- min(c(length(which(MPP_combined.sce$condition=="WT")), length(which(MPP_combined.sce$condition=="KO"))))
  WT_cells_to_plot <- sample(which(MPP_combined.sce$condition=="WT"), size = number_of_cells_to_plot_per_condition, replace = FALSE)
  KO_cells_to_plot <- sample(which(MPP_combined.sce$condition=="KO"), size = number_of_cells_to_plot_per_condition, replace = FALSE)
  MPP_combined.sce_sampled <- MPP_combined.sce[, c(WT_cells_to_plot, KO_cells_to_plot)]
  sds <- SlingshotDataSet(MPP_combined.sce)
  save <- readRDS(file=path_to_plot_parameters)
  
  # plot WT cells
  plotcol <- MPP_combined.sce_sampled[,which(MPP_combined.sce_sampled$condition=="WT")]$branches
  plotcol[which(plotcol=="CLP")] <- "orange"
  plotcol[which(plotcol=="GMP")] <- "blue"
  plotcol[which(plotcol=="MEP")] <- "red"
  plotcol[which(plotcol=="undefined")] <- "grey"
  plot3d(reducedDims(MPP_combined.sce_sampled[,which(MPP_combined.sce_sampled$condition=="WT")])$DiffusionMap_subset, col = plotcol, size = 4, alpha=0.5)
  plot3d.SlingshotDataSet(sds, lwd = 2, add = TRUE)
  par3d(save)
  snapshot3d(filename = paste(path_to_plots, "WT_diffusion_map.png", sep="/"), 
             fmt = "png", width = 5000, height = 5000)
  img <- image_read(paste(path_to_plots, "WT_diffusion_map.png", sep="/"))
  image_write(img, path = paste(path_to_plots, "WT_diffusion_map.pdf", sep="/"), format = "pdf")
  rgl.close()
  
  # plot KO cells
  plotcol <- MPP_combined.sce_sampled[,which(MPP_combined.sce_sampled$condition=="KO")]$branches
  plotcol[which(plotcol=="CLP")] <- "orange"
  plotcol[which(plotcol=="GMP")] <- "blue"
  plotcol[which(plotcol=="MEP")] <- "red"
  plotcol[which(plotcol=="undefined")] <- "grey"
  plot3d(reducedDims(MPP_combined.sce_sampled[,which(MPP_combined.sce_sampled$condition=="KO")])$DiffusionMap_subset, col = plotcol, size = 4, alpha=0.5)
  plot3d.SlingshotDataSet(sds, lwd = 2, add = TRUE)
  par3d(save)
  snapshot3d(filename = paste(path_to_plots, "KO_diffusion_map.png", sep="/"), 
             fmt = "png", width = 5000, height = 5000)
  img <- image_read(paste(path_to_plots, "KO_diffusion_map.png", sep="/"))
  image_write(img, path = paste(path_to_plots, "KO_diffusion_map.pdf", sep="/"), format = "pdf")
  rgl.close()
}


plot_differential_expression_of_signatures <- function(signature_genes=GMP_signature, MPP_combined.corrected_cell_numbers=MPP_combined.corrected_cell_numbers,
                                                       HSC_combined=HSC_combined,
                                                       plot_title="Differential expression\nof genes specific\nfor GMP-biased MPPs",
                                                       full_path="plots_20220216/GMP_signature_heatmap.pdf"){
  matrix_deg_signature_MPPs <- matrix(0, nrow=length(signature_genes), ncol = length(unique(MPP_combined.corrected_cell_numbers$branches)))
  colnames(matrix_deg_signature_MPPs) <- c("undefined", "CLP", "GMP", "MEP")
  rownames(matrix_deg_signature_MPPs) <- signature_genes
  for (index in (1:length(unique(MPP_combined.corrected_cell_numbers$branches)))){
    current_cluster <- unique(MPP_combined.corrected_cell_numbers$branches)[index]
    values_WT <- MPP_combined.corrected_cell_numbers@assays[["SCT"]]@counts[, which((MPP_combined.corrected_cell_numbers$branches==current_cluster) & (MPP_combined.corrected_cell_numbers$condition=="WT"))]
    values_KO <- MPP_combined.corrected_cell_numbers@assays[["SCT"]]@counts[, which((MPP_combined.corrected_cell_numbers$branches==current_cluster) & (MPP_combined.corrected_cell_numbers$condition=="KO"))]
    for (index2 in (1:nrow(matrix_deg_signature_MPPs))){
      current_gene <- rownames(matrix_deg_signature_MPPs)[index2]
      current_WT_values <- values_WT[current_gene,]
      current_KO_values <- values_KO[current_gene,]
      p_value <- wilcox.test(current_WT_values, current_KO_values)$p.value
      if ((is.numeric(p_value)==TRUE) & (!is.nan(p_value))){
        if (p_value<0.05){
          matrix_deg_signature_MPPs[current_gene, current_cluster] <- log2((mean(current_KO_values)+1)/(mean(current_WT_values)+1))
        } else {
          matrix_deg_signature_MPPs[current_gene, current_cluster] <- NA
        }
      } else {
        matrix_deg_signature_MPPs[current_gene, current_cluster] <- NA
      }
    }
  }
  colnames(matrix_deg_signature_MPPs) <- c("unbiased MPPs", "LMPPs", "GMP-biased MPPs", "MEP-biased MPPs")
  
  
  matrix_deg_signature_HSCs <- matrix(0, nrow=length(signature_genes), ncol = length(unique(HSC_combined$sample)))
  colnames(matrix_deg_signature_HSCs) <- c("LT-HSCs", "ST-HSCs")
  rownames(matrix_deg_signature_HSCs) <- signature_genes
  for (index in (1:length(unique(HSC_combined$sample)))){
    current_cluster <- unique(HSC_combined$sample)[index]
    values_WT <- HSC_combined@assays[["SCT"]]@counts[, which((HSC_combined$sample==current_cluster) & (HSC_combined$condition=="WT"))]
    values_KO <- HSC_combined@assays[["SCT"]]@counts[, which((HSC_combined$sample==current_cluster) & (HSC_combined$condition=="KO"))]
    for (index2 in (1:nrow(matrix_deg_signature_HSCs))){
      current_gene <- rownames(matrix_deg_signature_HSCs)[index2]
      if (current_gene %in% rownames(HSC_combined@assays[["SCT"]]@counts)){
        current_WT_values <- values_WT[current_gene,]
        current_KO_values <- values_KO[current_gene,]
        p_value <- wilcox.test(current_WT_values, current_KO_values)$p.value
        if ((is.numeric(p_value)==TRUE) & (!is.nan(p_value))){
          if (p_value<0.05){
            matrix_deg_signature_HSCs[current_gene, current_cluster] <- log2((mean(current_KO_values)+1)/(mean(current_WT_values)+1))
          } else {
            matrix_deg_signature_HSCs[current_gene, current_cluster] <- NA
          }
        } else {
          matrix_deg_signature_HSCs[current_gene, current_cluster] <- NA
        }
      } else {
        matrix_deg_signature_HSCs[current_gene, current_cluster] <- NA
      }
    }
  }
  
  
  
  combined_matrix <- cbind(matrix_deg_signature_HSCs, matrix_deg_signature_MPPs)
  
  
  
  color <- colorRampPalette(c("blue", "white", "red"))(100)
  myBreaks <- c(seq(min(combined_matrix, na.rm = TRUE), 0, length.out=ceiling(100/2) + 1), 
                seq(max(combined_matrix, na.rm = TRUE)/100, max(combined_matrix, na.rm = TRUE), length.out=floor(100/2)))
  pheatmap::pheatmap(combined_matrix, na_col = "grey", color = color, breaks = myBreaks, cluster_rows=FALSE, cluster_cols=FALSE,
                     angle_col=270, fontsize_row=8, fontsize_col=8, fontsize=8,
                     main = plot_title,
                     filename = full_path, width = 3, height = 7)
  
}

run_GSEA_and_group_significant_GO_terms_by_overlap_of_gene_sets <- function(seurat_object=MPP_combined.corrected_cell_numbers, significance_threshold=0.05, edge_threshold_Jaccard=0.3, min_number_of_leading_edge_genes=100, replicate_ID=NULL, cell_type=NULL, reuse_min_number_of_leading_edge_genes=NULL){
  GSEA_results_sampled <- run_GSEA(seurat_object=seurat_object, min.pct = 0.1, ont_category="BP", minGSSize = 10, maxGSSize = 800, pvalueCutoff = 0.05, showCategory=50)
  
  results_readable <- setReadable(GSEA_results_sampled[["results"]], 'org.Mm.eg.db', 'ENTREZID')
  
  GO_terms_to_consider <- results_readable@result$Description[which(results_readable@result$p.adjust < significance_threshold)]
  
  network_of_all_significant_terms <- emapplot(results_readable, showCategory = length(GO_terms_to_consider),          
                                               color = "p.adjust",
                                               layout = "nicely",
                                               cex_category = 0.2,
                                               cex_line = 0.2,
                                               min_edge = edge_threshold_Jaccard,
                                               cex_label_category = 0.4)
  
  GO_term_groups <- get_GO_term_groups(enrichment_object=results_readable, number_of_top_GO_terms_to_consider=NULL, GO_terms_to_consider=GO_terms_to_consider, edge_threshold_Jaccard=edge_threshold_Jaccard)
  
  core_genes_for_each_GO_term <- get_core_genes_for_each_GO_term(results_readable, GO_terms_of_interest=GO_terms_to_consider)
  
  number_of_unique_genes_per_group <- get_number_of_core_genes_per_GO_term_group(GO_term_groups=GO_term_groups, core_genes_for_each_GO_term=core_genes_for_each_GO_term)
  
  size_of_GO_term_groups <- unlist(lapply(GO_term_groups, length))
  number_of_leading_edge_genes_per_group <- number_of_unique_genes_per_group
  
  if (replicate_ID == 1){
    # save number of GO terms in largest group among those groups that will be excluded in the next step
    maximum_size_of_excluded_group_rep1 <- max(size_of_GO_term_groups[which(!(number_of_leading_edge_genes_per_group > min_number_of_leading_edge_genes))])
    saveRDS(object = maximum_size_of_excluded_group_rep1, file = paste0(paste("maximum_size_of_excluded_group_rep1", cell_type, sep = "_"), ".rds"))
  }
  
  if (replicate_ID != 1){
    if (reuse_min_number_of_leading_edge_genes == TRUE){
      min_number_of_leading_edge_genes <- readRDS(file=paste0(paste("min_number_of_leading_edge_genes", replicate_ID, sep = "_"), ".rds"))
    } else {
      # find a value for min_number_of_leading_edge_genes that causes retention of GO groups of similar size compared to replicate1
      maximum_size_of_excluded_group_rep1 <- readRDS(file = paste0(paste("maximum_size_of_excluded_group_rep1", cell_type, sep = "_"), ".rds"))
      threshold_values <- c(20:150)
      max_excluded_group_size <- c()
      for (index in (1:length(threshold_values))){
        current_threshold <- threshold_values[index]
        if (length(which(!(number_of_leading_edge_genes_per_group > current_threshold)))>0){
          max_excluded_group_size <- append(max_excluded_group_size, max(size_of_GO_term_groups[which(!(number_of_leading_edge_genes_per_group > current_threshold))]))
        } else {
          max_excluded_group_size <- append(max_excluded_group_size, 0)
        }
      }
      min_number_of_leading_edge_genes <- max(threshold_values[which(max_excluded_group_size <= maximum_size_of_excluded_group_rep1)])
      saveRDS(object = min_number_of_leading_edge_genes, file = paste0(paste("min_number_of_leading_edge_genes", replicate_ID, sep = "_"), ".rds"))
    }
  }
  
  GO_term_groups_cut_global <- GO_term_groups[which(number_of_unique_genes_per_group > min_number_of_leading_edge_genes)]
  
  return(list(results_readable=results_readable, GO_terms_to_consider=GO_terms_to_consider,
              GO_term_groups=GO_term_groups, core_genes_for_each_GO_term=core_genes_for_each_GO_term,
              number_of_unique_genes_per_group=number_of_unique_genes_per_group,
              GO_term_groups_cut_global=GO_term_groups_cut_global,
              network_of_all_significant_terms=network_of_all_significant_terms,
              GSEA_results_sampled=GSEA_results_sampled))
}

plot_filtered_GO_network <- function(processed_GSEA_results=processed_GSEA_results, min_edge=0.3, path_to_plots="plots_20220216", sample_name="MPPs", use_old_coordinates=FALSE, save_new_coordinates=FALSE, file_name_new_coordinates=NULL, file_name_old_coordinates=NULL, replicate_ID=1){
  basic_file_name <- paste(path_to_plots, sample_name, sep="/")
  # create filtered gsea-results object
  retained_GO_terms <- unique(unlist(processed_GSEA_results[["GO_term_groups_cut_global"]]))
  results_global_filtered <- processed_GSEA_results[["results_readable"]]
  results_global_filtered@result <- results_global_filtered@result[which(results_global_filtered@result$Description %in% retained_GO_terms),]
  # plot global filtered GO-term groups
  grouped_GO_terms_plot_labeled <- emapplot(results_global_filtered, showCategory = length(retained_GO_terms),          
                                            color = "p.adjust",
                                            layout = "nicely",
                                            cex_category = 0.2,
                                            cex_line = 0.2,
                                            min_edge = min_edge,
                                            cex_label_category = 0.2)
  
  if(use_old_coordinates==TRUE){
    old_coordinates <- readRDS(file=file_name_old_coordinates)
    grouped_GO_terms_plot_labeled[["data"]][["x"]] <- old_coordinates[,"x"]
    grouped_GO_terms_plot_labeled[["data"]][["y"]] <- old_coordinates[,"y"]
  }
  
  ggsave(
    paste(basic_file_name, paste0(paste("GO_network_labeled", replicate_ID, sep = "_"), ".pdf"), sep="_"),
    plot = grouped_GO_terms_plot_labeled,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 20,
    height = 20,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
  
  coordinates_of_labeled_plot <- data.frame(x=grouped_GO_terms_plot_labeled[["data"]][["x"]], y=grouped_GO_terms_plot_labeled[["data"]][["y"]])
  
  if(save_new_coordinates==TRUE){
    saveRDS(object=coordinates_of_labeled_plot, file=file_name_new_coordinates)
  }
  
  grouped_GO_terms_plot <- emapplot(results_global_filtered, showCategory = length(retained_GO_terms),          
                                    color = "p.adjust",
                                    layout = "nicely",
                                    cex_category = 0.2,
                                    cex_line = 0.2,
                                    min_edge = min_edge,
                                    cex_label_category = 0.0)
  
  grouped_GO_terms_plot[["data"]][["x"]] <- coordinates_of_labeled_plot[,"x"]
  grouped_GO_terms_plot[["data"]][["y"]] <- coordinates_of_labeled_plot[,"y"]
  
  ggsave(
    paste(basic_file_name, paste0(paste("GO_network_unlabeled", replicate_ID, sep = "_"), ".pdf"), sep="_"),
    plot = grouped_GO_terms_plot,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 20,
    height = 20,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
  
}

visualize_subclustering_results <- function(number_of_clusters=suggested_number_of_clusters, subclustering_results=subclustering_results, processed_GSEA_results=processed_GSEA_results, min_edge=0.3){
  list_of_clustered_immune_GO_terms <- cut_tree_and_get_list_of_clustered_GO_terms(hclust_ward=subclustering_results[["hclust_ward"]], number_of_clusters=number_of_clusters)
  results_readable_clust <- processed_GSEA_results[["results_readable"]]
  results_readable_clust@result <- results_readable_clust@result[which(results_readable_clust@result$Description %in% unlist(list_of_clustered_immune_GO_terms)),]
  
  GO_term_vector <- c()
  Cluster_index_vector <- c()
  for (index in (1:length(list_of_clustered_immune_GO_terms))){
    GO_term_vector <- append(GO_term_vector, list_of_clustered_immune_GO_terms[[index]])
    Cluster_index_vector <- append(Cluster_index_vector, rep(index, length(list_of_clustered_immune_GO_terms[[index]])))
  }
  names(Cluster_index_vector) <- GO_term_vector
  results_readable_clust@result[["GO_subcluster"]] <- Cluster_index_vector[results_readable_clust@result$Description]
  
  
  grouped_GO_terms_plot <- emapplot(results_readable_clust, showCategory = nrow(results_readable_clust@result),
                                    color = "GO_subcluster",
                                    layout = "nicely",
                                    cex_category = 0.2,
                                    cex_line = 0.2,
                                    min_edge = min_edge,
                                    cex_label_category = 0.4)
  return(list(grouped_GO_terms_plot=grouped_GO_terms_plot, list_of_clustered_immune_GO_terms=list_of_clustered_immune_GO_terms))
}

plot_GSEA_heatmap_subcluster <- function(processed_GSEA_results=processed_GSEA_results,
                                         list_of_clustered_immune_GO_terms=subclusters_and_visualization[["list_of_clustered_immune_GO_terms"]], 
                                         names_of_subclusters=names_of_subclusters, mother_cluster="immune processes", log2FC_threshold=1.2, 
                                         number_of_clusters=suggested_number_of_clusters,
                                         path_to_plot="plots_20220216"){
  
  core_genes_for_each_GO_term <- processed_GSEA_results[["core_genes_for_each_GO_term"]]
  GO_term_groups_cut_global <- processed_GSEA_results[["GO_term_groups_cut_global"]]
  named_foldchanges_sorted=processed_GSEA_results[["GSEA_results_sampled"]]$readable_foldchanges
  
  names(list_of_clustered_immune_GO_terms) <- names_of_subclusters
  all_immune_core_genes <- unique(unlist(core_genes_for_each_GO_term[GO_term_groups_cut_global[[mother_cluster]]]))
  all_immune_core_genes_sorted_FCs <- named_foldchanges_sorted[which(names(named_foldchanges_sorted) %in% all_immune_core_genes)]
  all_immune_core_genes_sorted_FCs_filtered <- all_immune_core_genes_sorted_FCs[which((all_immune_core_genes_sorted_FCs > log2(log2FC_threshold)) | (all_immune_core_genes_sorted_FCs < (-log2(log2FC_threshold))))]
  immune_matrix <- matrix(0, ncol = number_of_clusters, nrow = length(all_immune_core_genes_sorted_FCs_filtered))
  for (index in (1:number_of_clusters)){
    immune_matrix[,index] <- all_immune_core_genes_sorted_FCs_filtered
  }
  rownames(immune_matrix) <- names(all_immune_core_genes_sorted_FCs_filtered)
  colnames(immune_matrix) <- names(list_of_clustered_immune_GO_terms)
  # set NAs
  for (index in (1:ncol(immune_matrix))){
    current_GO_terms <- list_of_clustered_immune_GO_terms[[index]]
    current_core_genes <- unique(unlist(core_genes_for_each_GO_term[current_GO_terms]))
    immune_matrix[which(!(rownames(immune_matrix) %in% current_core_genes)), index] <- NA
  }
  
  #  mean_values_of_columns <- apply(immune_matrix, 2, FUN=function(x)mean(x, na.rm=TRUE))
  genes_per_column <- apply(immune_matrix, 2, FUN=function(x)(length(which(!is.na(x)))))
  #  columns_with_negative_values <- which(mean_values_of_columns<0)
  #  negative_columns_gene_numbers <- genes_per_column[columns_with_negative_values]
  #  negative_sorted <- sort(negative_columns_gene_numbers, decreasing = FALSE)
  
  #  columns_with_positive_values <- which(mean_values_of_columns>0)
  #  positive_columns_gene_numbers <- genes_per_column[columns_with_positive_values]
  genes_per_column_sorted <- sort(genes_per_column, decreasing = TRUE)
  
  immune_matrix <- immune_matrix[, names(genes_per_column_sorted)]
  
  color <- colorRampPalette(c("white", "red"))(100)
  myBreaks <- seq(0, max(immune_matrix, na.rm = TRUE), length.out=100)
  pheatmap::pheatmap(immune_matrix, na_col = "grey", color = color, breaks = myBreaks, cluster_cols = FALSE, cluster_rows = FALSE, 
                     angle_col=270, fontsize_row=6, fontsize_col=8, fontsize=8, 
                     main = "", 
                     filename = paste(path_to_plot, "GO_heatmap_immune.pdf", sep="/"), width = 2, height = 7)
  
}

# replicate-specific parameters: 
load_replicate_parameters <- function(replicate=1){
  if (replicate==1){
    parameters <- list(replicate=1,
                       minimum_number_of_genes_per_cell_MPP=2500,
                       minimum_number_of_genes_per_cell_LT=2500,
                       minimum_number_of_genes_per_cell_ST=2500,
                       file_name_GO_network_coordinates_MPP="MPPs_emapplot_coordinates_20220316",
                       file_name_GO_network_coordinates_LT="LTs_emapplot_coordinates_20220316",
                       parameter_file_for_3D_plot_MPP="3d_scatter_parameters_20220312.rds")
  }
  if (replicate==2){
    parameters <- list(replicate=2,
                       minimum_number_of_genes_per_cell_MPP=1500,
                       minimum_number_of_genes_per_cell_LT=1500,
                       minimum_number_of_genes_per_cell_ST=1500,
                       file_name_GO_network_coordinates_MPP="MPPs_emapplot_coordinates_2",
                       file_name_GO_network_coordinates_LT="LTs_emapplot_coordinates_2",
                       parameter_file_for_3D_plot_MPP="3d_scatter_parameters_2.rds")
  }
  return(parameters)
}

# make gene length function reproducible
getGeneLengthAndGCContent_reproducible <- function (id, org, mode = c("biomart", "org.db")) 
{
  id.type <- .autoDetectGeneIdType(id[1])
  if (is.na(id.type)) 
    stop("Only ENTREZ or ENSEMBL gene IDs are supported.")
  mode <- match.arg(mode)
  inp.id <- id
  if (mode == "org.db") {
    txdb.pkg <- .org2pkg(org, type = "TxDb")
    .isAvailable(txdb.pkg, type = "TxDb")
    bsgen.pkg <- .org2pkg(org, type = "BSgenome")
    .isAvailable(bsgen.pkg, type = "BSgenome")
    txdb.spl <- unlist(strsplit(txdb.pkg, "\\."))
    txdb.id.type <- txdb.spl[length(txdb.spl)]
    if (txdb.id.type == "ensGene") {
      txdb.id.type <- "ensembl"
    }
    else if (txdb.id.type == "knownGene") {
      txdb.id.type <- "entrez"
    }
    else if (txdb.id.type == "sgdGene") {
      txdb.id.type <- "sgd"
    }
    else {
      stop(paste("TxDb does not use ENSEMBL or ENTREZ gene IDs"))
    }
    if (id.type != txdb.id.type) {
      orgdb.pkg <- .org2pkg(org)
      .isAvailable(orgdb.pkg)
      orgdb.pkg <- get(orgdb.pkg)
      id.map <- mapIds(orgdb.pkg, keys = id, column = ifelse(id.type == 
                                                               "entrez", "ENSEMBL", "ENTREZID"), keytype = ifelse(id.type == 
                                                                                                                    "entrez", "ENTREZID", "ENSEMBL"))
      id <- id.map[!is.na(id.map)]
    }
    txdb.pkg <- get(txdb.pkg)
    coords <- exonsBy(txdb.pkg, by = "gene")
    id <- id[id %in% names(coords)]
    coords <- reduce(coords[id])
    len <- sum(width(coords))
    bsgen.pkg <- get(bsgen.pkg)
    seqs <- getSeq(bsgen.pkg, coords)
    af <- alphabetFrequency(unlist(seqs, use.names = FALSE), 
                            baseOnly = TRUE, as.prob = TRUE)
    gc.cont <- mean(relist(rowSums(af[, c("C", "G")]), seqs))
  }
  else {
    id.type <- paste0(id.type, ifelse(id.type == "entrez", 
                                      "gene", "_gene_id"))
    message("Connecting to BioMart ...")
    ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "https://dec2021.archive.ensembl.org")
    ds <- listDatasets(ensembl)[, "dataset"]
    ds <- grep(paste0("^", org), ds, value = TRUE)
    if (length(ds) == 0) 
      stop(paste("Mart not found for:", org))
    else if (length(ds) > 1) {
      message("Found several marts")
      sapply(ds, function(d) message(paste(which(ds == 
                                                   d), d, sep = ": ")))
      n <- readline(paste0("Choose mart (1-", length(ds), 
                           ") : "))
      ds <- ds[as.integer(n)]
    }
    ensembl <- useDataset(ds, mart = ensembl)
    message(paste0("Downloading sequence", ifelse(length(id) > 
                                                    1, "s", ""), " ..."))
    if (length(id) > 100) 
      message("This may take a few minutes ...")
    attrs <- c(id.type, "ensembl_exon_id", "chromosome_name", 
               "exon_chrom_start", "exon_chrom_end")
    coords <- getBM(filters = id.type, attributes = attrs, 
                    values = id, mart = ensembl)
    id <- unique(coords[, id.type])
    coords <- GRangesList(sapply(id, function(i) {
      i.coords <- coords[coords[, 1] == i, 3:5]
      g <- GRanges(i.coords[, 1], IRanges(i.coords[, 2], 
                                          i.coords[, 3]))
      return(g)
    }), compress = FALSE)
    coords <- reduce(coords)
    len <- sum(width(coords))
    sel <- c(id.type, "start_position", "end_position")
    gene.pos <- getBM(attributes = sel, filters = id.type, 
                      values = id, mart = ensembl)
    gene.seqs <- getSequence(id = id, type = id.type, seqType = "gene_exon_intron", 
                             mart = ensembl)
    gc.cont <- sapply(id, function(i) {
      ecoords <- coords[[i]]
      gpos <- gene.pos[gene.pos[, id.type] == i, c("start_position", 
                                                   "end_position")]
      gseq <- DNAString(gene.seqs[gene.seqs[, id.type] == 
                                    i, "gene_exon_intron"])
      start <- start(ranges(ecoords)) - gpos[1, 1] + 1
      end <- end(ranges(ecoords)) - gpos[1, 1] + 1
      eseq <- gseq[IRanges(start, end)]
      gc.cont <- sum(alphabetFrequency(eseq, as.prob = TRUE)[c("C", 
                                                               "G")])
      return(gc.cont)
    })
  }
  res <- cbind(len, gc.cont)
  colnames(res) <- c("length", "gc")
  rownames(res) <- id
  if (mode == "org.db") 
    if (id.type != txdb.id.type) 
      rownames(res) <- names(id)
  not.found <- !(inp.id %in% rownames(res))
  na.col <- rep(NA, sum(not.found))
  rn <- c(rownames(res), inp.id[not.found])
  res <- rbind(res, cbind(na.col, na.col))
  rownames(res) <- rn
  res <- res[inp.id, ]
  return(res)
}
environment(getGeneLengthAndGCContent_reproducible) <- asNamespace('EDASeq')

# construct paths for saving and reading intermediate results
get_paths_for_intermediate_results <- function(replicate_ID = NULL, paths_to_results=NULL, prefixMPPs="MPPs_after_annotation", prefixTBM="TBM_after_annotation", prefixHSC="LT_ST_after_merge"){
  MPP_WT_path <- paste0(paste(paste0(paths_to_results, replicate_ID), prefixMPPs, sep="/"), "_WT.rds")
  MPP_KO_path <- paste0(paste(paste0(paths_to_results, replicate_ID), prefixMPPs, sep="/"), "_KO.rds")
  MPP_combined_after_slingshot <- paste0(paste0(paths_to_results, replicate_ID), "/MPP_combined_after_slingshot.rds")
  MPP_combined.sce_after_slingshot <- paste0(paste0(paths_to_results, replicate_ID), "/MPP_combined.sce_after_slingshot.rds")
  MPP_combined_after_sampling <- paste0(paste0(paths_to_results, replicate_ID), "/MPP_combined_after_sampling.rds")
  MPP_combined.sce_after_sampling <- paste0(paste0(paths_to_results, replicate_ID), "/MPP_combined.sce_after_sampling.rds")
  TBM_combined_path <- paste0(paste(paste0(paths_to_results, replicate_ID), prefixTBM, sep="/"), ".rds")
  HSC_combined_path <- paste0(paste(paste0(paths_to_results, replicate_ID), prefixHSC, sep="/"), ".rds")
  return(list(MPP_WT_path=MPP_WT_path,
              MPP_KO_path=MPP_KO_path,
              TBM_combined_path=TBM_combined_path,
              MPP_combined_after_slingshot=MPP_combined_after_slingshot,
              MPP_combined.sce_after_slingshot=MPP_combined.sce_after_slingshot,
              MPP_combined_after_sampling=MPP_combined_after_sampling,
              MPP_combined.sce_after_sampling=MPP_combined.sce_after_sampling,
              HSC_combined_path=HSC_combined_path))
}

create_seurat_object_with_gene_detection_rates_in_MPP_subsets <- function(MPP=NULL){
  MPP@meta.data[["sample"]] <- paste(paste(MPP$condition, MPP$experiment, sep="_"), MPP$branches, sep = "_")
  sample_names <- as.character(unique(MPP$sample))
  detection_rate_matrix <- matrix(0, nrow = nrow(MPP@assays$RNA@counts), ncol=length(sample_names))
  rownames(detection_rate_matrix) <- rownames(MPP@assays$RNA@counts)
  colnames(detection_rate_matrix) <- sample_names
  for (index in (1:length(sample_names))){
    current_sample <- sample_names[index]
    
    number_of_cells_in_current_sample <- length(which(MPP$sample == current_sample))
    numbers_of_zeros <- proxyC::rowZeros(MPP@assays$RNA@counts[, which(MPP$sample == current_sample)])
    numbers_of_detections <- number_of_cells_in_current_sample - numbers_of_zeros
    detection_rate_vector <- numbers_of_detections / number_of_cells_in_current_sample
    detection_rate_matrix[, current_sample] <- detection_rate_vector
  }
  
  Seurat_detection_MPP <- CreateSeuratObject(counts=detection_rate_matrix)
  list_with_substrings = strsplit(colnames(Seurat_detection_MPP), '_')
  condition_str <- unlist(lapply(list_with_substrings, FUN=function(x){x[1]}))
  experiment_str <- unlist(lapply(list_with_substrings, FUN=function(x){x[2]}))
  branch_str <- unlist(lapply(list_with_substrings, FUN=function(x){x[3]}))
  Seurat_detection_MPP@meta.data[["condition"]] <- condition_str
  Seurat_detection_MPP@meta.data[["experiment"]] <- experiment_str
  Seurat_detection_MPP@meta.data[["branches"]] <- branch_str
  Seurat_detection_MPP$condition <- factor(Seurat_detection_MPP$condition, levels = c("WT", "KO"))
  
  # rename MPP branches
  Seurat_detection_MPP$branches[which(Seurat_detection_MPP$branches=="CLP")] <- "LMPP"
  Seurat_detection_MPP$branches[which(Seurat_detection_MPP$branches=="GMP")] <- "GMP-biased MPP"
  Seurat_detection_MPP$branches[which(Seurat_detection_MPP$branches=="MEP")] <- "MEP-biased MPP"
  Seurat_detection_MPP$branches[which(Seurat_detection_MPP$branches=="undefined")] <- "unbiased MPP"
  
  return(Seurat_detection_MPP)
}

create_seurat_object_with_average_expression_per_gene_in_MPP_subsets <- function(MPP=NULL, merge_sample_subsets=FALSE){
  if (merge_sample_subsets==FALSE){
    MPP@meta.data[["sample"]] <- paste(paste(MPP$condition, MPP$experiment, sep="_"), MPP$branches, sep = "_")
  } else {
    MPP@meta.data[["sample"]] <- paste(paste(MPP$condition, MPP$experiment, sep="_"), MPP$sample, sep = "_")
  }
  sample_names <- as.character(unique(MPP$sample))
  average_expression_matrix <- matrix(0, nrow = nrow(MPP@assays$SCT@counts), ncol=length(sample_names))
  rownames(average_expression_matrix) <- rownames(MPP@assays$SCT@counts)
  colnames(average_expression_matrix) <- sample_names
  for (index in (1:length(sample_names))){
    current_sample <- sample_names[index]
    sum_of_all_expression_values <- rowSums(MPP@assays$SCT@counts[, which(MPP$sample == current_sample)])
    average_expression_matrix[, current_sample] <- sum_of_all_expression_values/length(which(MPP$sample == current_sample))
  }
  
  Seurat_detection_MPP <- CreateSeuratObject(counts=average_expression_matrix)
  list_with_substrings = strsplit(colnames(Seurat_detection_MPP), '_')
  condition_str <- unlist(lapply(list_with_substrings, FUN=function(x){x[1]}))
  experiment_str <- unlist(lapply(list_with_substrings, FUN=function(x){x[2]}))
  branch_str <- unlist(lapply(list_with_substrings, FUN=function(x){x[3]}))
  Seurat_detection_MPP@meta.data[["condition"]] <- condition_str
  Seurat_detection_MPP@meta.data[["experiment"]] <- experiment_str
  Seurat_detection_MPP@meta.data[["branches"]] <- branch_str
  Seurat_detection_MPP$condition <- factor(Seurat_detection_MPP$condition, levels = c("WT", "KO"))
  
  # rename MPP branches
  if (merge_sample_subsets==FALSE){
    Seurat_detection_MPP$branches[which(Seurat_detection_MPP$branches=="CLP")] <- "LMPP"
    Seurat_detection_MPP$branches[which(Seurat_detection_MPP$branches=="GMP")] <- "GMP-biased MPP"
    Seurat_detection_MPP$branches[which(Seurat_detection_MPP$branches=="MEP")] <- "MEP-biased MPP"
    Seurat_detection_MPP$branches[which(Seurat_detection_MPP$branches=="undefined")] <- "unbiased MPP"
  } else {
    Seurat_detection_MPP$branches[which(Seurat_detection_MPP$branches %in% c("CLP", "GMP", "MEP", "undefined"))] <- "MPP"
  }
  return(Seurat_detection_MPP)
}



compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice <- function(seurat_object=NULL, gene_set=NULL, folder_with_plots=NULL, 
                                                                               plot_width=25, plot_hight=15, 
                                                                               y_axis_label="% of cells with detectable mRNA", 
                                                                               file_name="gene_plots_patch_NK_TFs.pdf",
                                                                               mode="detection_frequency",
                                                                               colnumber_panel=NULL,
                                                                               rownumber_panel=NULL,
                                                                               pt_size=1.2,
                                                                               line_thickness=0.4){
  genes_to_visualize <- gene_set
  list_with_gene_plots <- list()
  for (index in (1:length(genes_to_visualize))){
    current_gene <- genes_to_visualize[index]
    if ((current_gene %in% rownames(seurat_object@assays$RNA@counts)) == FALSE){
      print(paste(current_gene, "not in data", sep = " "))
      next
    }
    data <- seurat_object@meta.data[, c("branches", "condition", "experiment")]
    data[, "gene"] <- seurat_object@assays$RNA@counts[current_gene,]
    data[, "pair"] <- substr(rownames(data), start = 4, stop = nchar(rownames(data)))
    if (mode=="detection_frequency"){
      data$gene <- data$gene*100
    }
    data$branches <- as.character(data$branches)
    data$branches[which(data$branches=="GMP-biased MPP")] <- "GMP-\nbiased\nMPP"
    data$branches[which(data$branches=="MEP-biased MPP")] <- "MEP-\nbiased\nMPP"
    data$branches[which(data$branches=="unbiased MPP")] <- "unbiased\nMPP"
    data$branches[which(data$branches=="LT-HSCs")] <- "LT"
    data$branches[which(data$branches=="ST-HSCs")] <- "ST"
    data$branches[which(data$branches=="Granulocytes")] <- "Granulo-\ncytes"
    data$branches[which(data$branches=="Monocytes")] <- "Mono-\ncytes"
    data$branches <- factor(data$branches, levels = c("LT", "ST", "MPP", "GMP-\nbiased\nMPP", "LMPP", "MEP-\nbiased\nMPP", "unbiased\nMPP",
                                                      "Granulo-\ncytes", "Mono-\ncytes", "DC", "B-cells", "T", "NK", "EryBl", "Retic"))
    
    p <- ggpaired(data, x = "condition", y = "gene",
                  color = "condition", palette = c("black", "red"), 
                  line.color = "gray", line.size = line_thickness,
                  facet.by = "branches", short.panel.labs = FALSE, id="pair")+ 
      ylab(y_axis_label)+ggtitle(current_gene)+
      geom_point(aes(shape=experiment, color=condition), size=pt_size)+scale_shape_manual(values=c(17, 16))+#scale_size_manual(values = c(pt_size, pt_size))+
      theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 8), axis.text.x = element_text(size=8), 
            axis.text.y = element_text(size=8), plot.title = element_text(hjust = 0.5), legend.position="right")
    p <- facet(p, facet.by = "branches", nrow=1, ncol=length(unique(data$branches)), panel.labs.font = list(face = NULL, color = NULL, size = 8, angle = 0))
    p[["layers"]][[1]] <- NULL
    final_plot <- p # + stat_compare_means(method="wilcox.test", method.args = list(alternative = "two.sided"), label = "p.format", paired = TRUE, size=2) + stat_summary(fun=mean, geom="point", shape=23, size=4)
    list_with_gene_plots <- list.append(list_with_gene_plots, final_plot)
  }
  gene_plots_patch <- wrap_plots(list_with_gene_plots, guides = 'collect', ncol=colnumber_panel, nrow=rownumber_panel)
  ggsave(
    paste(folder_with_plots, file_name, sep="/"),
    plot = gene_plots_patch,
    device = NULL,
    path = NULL,
    scale = 1,
    width = plot_width,
    height = plot_hight,
    units = "cm",
    dpi = 800,
    limitsize = TRUE,
    bg = NULL
  )
}




