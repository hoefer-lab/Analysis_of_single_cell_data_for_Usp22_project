source("script2.R", local=FALSE)
path_to_working_directory <- getwd()
folder_with_plots <- "plots_of_individual_genes_20220817"
dir.create(paste(path_to_working_directory, folder_with_plots, sep="/"), showWarnings = FALSE)



# create gene sets ####
# to analyse in all LSK cells (including MPP subsets)
genes_mobility_in_LSK <- c("Cxcr4", "Kit", "Flt3", "Cd44", "Itga4", "Itga5", "Itgb1", "Efnb2")
myeloid_TFs <- c("Cebpa", "Cebpb", "Cebpe", "Egr1", "Egr2", "Gfi1", "Hoxa10", "Irf8", "Spi1")
NK_cell_lineage_committment_markers <- c("Gzmb", "Klrb1b", "Klrb1c", "Klre1", "Klrk1", "Il2rb")
NK_cell_maturation_TFs <- c("Tbx21", "Eomes", "Tox", "Tox2", "Prdm1", "Zeb2", "Gata3", "Smad4", "Smad3", "Foxo1")
NK_cell_lineage_committment_TFs <- c("Nfil3", "Tcf7", "Ets1", "Id2", "Stat5", "Id3", "Spi1")
DEGs <- c("Txnip", "Mif", "Cd74")
TFs_to_validate_branch_identity <- c("Spi1", "Tcf3", "Tal1", "Klf1", 
                                     "Ikzf1", "Gata1", "Cebpa", "Cebpb", 
                                     "Cebpe", "Egr1", "Egr2", "Gfi1", 
                                     "Hoxa10", "Irf8", "Ebf1", "Bcl11a", 
                                     "Foxo1", "Notch1", "Tcf12","Tcf7", 
                                     "Gata3", "Nfil3", "Ets1", "Id2")
interleukins <- c("Il1a", "Il1b", "Il6", "Ifng")

# create plots for LSKs ####
HSC1 <- readRDS(file = "/Users/metz/Nextcloud/backup_Usp22_data_and_Robjects/R_objects_Usp22/objects_with_intermediate_results/1/LT_ST_after_merge.rds")
HSC2 <- readRDS(file = "/Users/metz/Nextcloud/backup_Usp22_data_and_Robjects/R_objects_Usp22/objects_with_intermediate_results/2/LT_ST_after_merge.rds")
HSC <- merge(HSC1, y = c(HSC2), add.cell.ids = NULL, project = "U22")
HSC@meta.data[["experiment"]] <- c(rep("1", ncol(HSC1)), rep("2", ncol(HSC2)))
rm(HSC1)
rm(HSC2)
DefaultAssay(HSC) <- "RNA"
HSC@assays$SCT <- NULL
saveRDS(object=HSC, file="HSC_all_experiments")
rm(HSC)

MPP1 <- readRDS(file = "/Users/metz/Nextcloud/backup_Usp22_data_and_Robjects/R_objects_Usp22/objects_with_intermediate_results/1/MPP_combined_after_slingshot.rds")
MPP2 <- readRDS(file = "/Users/metz/Nextcloud/backup_Usp22_data_and_Robjects/R_objects_Usp22/objects_with_intermediate_results/2/MPP_combined_after_slingshot.rds")
MPP <- merge(MPP1, y = c(MPP2), add.cell.ids = NULL, project = "U22")
MPP@meta.data[["experiment"]] <- c(rep("1", ncol(MPP1)), rep("2", ncol(MPP2)))
rm(MPP1)
rm(MPP2)
DefaultAssay(MPP) <- "RNA"
MPP@assays$SCT <- NULL
saveRDS(object=MPP, file="MPP_all_experiments")

HSC <- readRDS(file="HSC_all_experiments")
LSK <- merge(HSC, y = c(MPP), add.cell.ids = NULL, project = "U22")
rm(HSC)
rm(MPP)
LSK$branches[which(is.na(LSK$branches))] <- LSK$sample[which(is.na(LSK$branches))]

# The SCTransform function is used in the conserve.memory mode and the vst function which is internally called in line 114 returns the warning that no corrected counts are returned
# This warning can be neglected because corrected counts are computed in line 172 as the variable residual.type was not overwritten and is still set to pearson
DefaultAssay(LSK) <- "RNA"
LSK <- Seurat::SCTransform(LSK, verbose = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3, conserve.memory = TRUE, do.correct.umi = TRUE)
saveRDS(object=LSK, file="LSK_all_experiments_normalized_together.rds")

# 1. for TFs regulating NK cell lineage commitment
Seurat_detection_LSK <- create_seurat_object_with_gene_detection_rates_in_MPP_subsets(MPP=LSK)
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_detection_LSK, gene_set=NK_cell_lineage_committment_TFs, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=30, plot_hight=15, 
                                                                   y_axis_label="% of cells with detectable mRNA", 
                                                                   file_name="gene_plots_patch_NK_TFs_LSK.pdf",
                                                                   mode="detection_frequency",
                                                                   pt_size=3.0,
                                                                   line_thickness=0.6)


Seurat_average_expression_LSK <- create_seurat_object_with_average_expression_per_gene_in_MPP_subsets(MPP=LSK)
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_average_expression_LSK, gene_set=NK_cell_lineage_committment_TFs, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=30, plot_hight=15, 
                                                                   y_axis_label="average of normalized mRNA counts", 
                                                                   file_name="gene_plots_patch_NK_TFs_average_corrected_counts_LSK.pdf",
                                                                   mode="gene_expression",
                                                                   pt_size=3.0,
                                                                   line_thickness=0.6)

# 2. for NK_cell_maturation_TFs
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_detection_LSK, gene_set=NK_cell_maturation_TFs, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=40, plot_hight=25, 
                                                                   y_axis_label="% of cells with detectable mRNA", 
                                                                   file_name="gene_plots_patch_NK_maturation_TFs_LSK.pdf",
                                                                   mode="detection_frequency",
                                                                   pt_size=3.0,
                                                                   line_thickness=0.6)
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_average_expression_LSK, gene_set=NK_cell_maturation_TFs, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=40, plot_hight=25, 
                                                                   y_axis_label="average of normalized mRNA counts", 
                                                                   file_name="gene_plots_patch_NK_maturation_TFs_average_corrected_counts_LSK.pdf",
                                                                   mode="gene_expression",
                                                                   pt_size=3.0,
                                                                   line_thickness=0.6)

# 3. for NK_cell_lineage_committment_markers
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_detection_LSK, gene_set=NK_cell_lineage_committment_markers, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=30, plot_hight=15, 
                                                                   y_axis_label="% of cells with detectable mRNA", 
                                                                   file_name="gene_plots_patch_NK_markers_LSK.pdf",
                                                                   mode="detection_frequency",
                                                                   pt_size=3.0,
                                                                   line_thickness=0.6)
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_average_expression_LSK, gene_set=NK_cell_lineage_committment_markers, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=30, plot_hight=15, 
                                                                   y_axis_label="average of normalized mRNA counts", 
                                                                   file_name="gene_plots_patch_NK_markers_average_corrected_counts_LSK.pdf",
                                                                   mode="gene_expression",
                                                                   pt_size=3.0,
                                                                   line_thickness=0.6)

# 4. for myeloid_TFs
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_detection_LSK, gene_set=myeloid_TFs, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=30, plot_hight=22.5, 
                                                                   y_axis_label="% of cells with detectable mRNA", 
                                                                   file_name="gene_plots_patch_myeloid_TFs_LSK.pdf",
                                                                   mode="detection_frequency",
                                                                   pt_size=3.0,
                                                                   line_thickness=0.6)
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_average_expression_LSK, gene_set=myeloid_TFs, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=30, plot_hight=22.5, 
                                                                   y_axis_label="average of normalized mRNA counts", 
                                                                   file_name="gene_plots_patch_myeloid_TFs_average_corrected_counts_LSK.pdf",
                                                                   mode="gene_expression",
                                                                   pt_size=3.0,
                                                                   line_thickness=0.6)

# 5. for genes_mobility_in_LSK
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_detection_LSK, gene_set=genes_mobility_in_LSK, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=30, plot_hight=22.5, 
                                                                   y_axis_label="% of cells with detectable mRNA", 
                                                                   file_name="gene_plots_patch_genes_mobility_LSK.pdf",
                                                                   mode="detection_frequency",
                                                                   pt_size=3.0,
                                                                   line_thickness=0.6)
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_average_expression_LSK, gene_set=genes_mobility_in_LSK, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=30, plot_hight=22.5, 
                                                                   y_axis_label="average of normalized mRNA counts", 
                                                                   file_name="gene_plots_patch_genes_mobility_average_corrected_counts_LSK.pdf",
                                                                   mode="gene_expression",
                                                                   pt_size=3.0,
                                                                   line_thickness=0.6)

# plot Important DEGs
DEGs <- c("Txnip", "Mif", "Cd74")
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_average_expression_LSK, gene_set=DEGs, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=30, plot_hight=7.5, 
                                                                   y_axis_label="average of normalized mRNA counts", 
                                                                   file_name="DEGs_LSK.pdf",
                                                                   mode="gene_expression",
                                                                   pt_size=3.0,
                                                                   line_thickness=0.6)
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_detection_LSK, gene_set=DEGs, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=30, plot_hight=7.5, 
                                                                   y_axis_label="% of cells with detectable mRNA", 
                                                                   file_name="DEGs_LSK_detection.pdf",
                                                                   mode="detection_frequency",
                                                                   pt_size=3.0,
                                                                   line_thickness=0.6)

# show all branch specific TFs
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_average_expression_LSK, gene_set=TFs_to_validate_branch_identity, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=45, plot_hight=45, 
                                                                   y_axis_label="average of normalized mRNA counts", 
                                                                   file_name="TFs_all_lineages_LSK.pdf",
                                                                   mode="gene_expression",
                                                                   colnumber_panel=4,
                                                                   rownumber_panel=6,
                                                                   pt_size=3.0,
                                                                   line_thickness=0.6)

# plot interleukin expression
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_average_expression_LSK, gene_set=interleukins, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=30, plot_hight=15, 
                                                                   y_axis_label="average of normalized mRNA counts", 
                                                                   file_name="interleukin_average_expression_in_LSK.pdf",
                                                                   mode="gene_expression",
                                                                   pt_size=3.0,
                                                                   line_thickness=0.6)
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_detection_LSK, gene_set=interleukins, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=30, plot_hight=15, 
                                                                   y_axis_label="% of cells with detectable mRNA", 
                                                                   file_name="interleukin_detection_frequency_in_LSK.pdf",
                                                                   mode="detection_frequency",
                                                                   pt_size=3.0,
                                                                   line_thickness=0.6)


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


# 6. analyse mobility factors in TBM ####
TBM1 <- readRDS(file = "/Users/metz/Nextcloud/backup_Usp22_data_and_Robjects/R_objects_Usp22/objects_with_intermediate_results/1/TBM_after_annotation.rds")
TBM2 <- readRDS(file = "/Users/metz/Nextcloud/backup_Usp22_data_and_Robjects/R_objects_Usp22/objects_with_intermediate_results/2/TBM_after_annotation.rds")
TBM <- merge(TBM1, y = c(TBM2), add.cell.ids = NULL, project = "U22")
TBM@meta.data[["experiment"]] <- c(rep("1", ncol(TBM1)), rep("2", ncol(TBM2)))
rm(TBM1)
rm(TBM2)
DefaultAssay(TBM) <- "RNA"
TBM@assays$SCT <- NULL
TBM <- Seurat::SCTransform(TBM, verbose = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3, conserve.memory = TRUE, do.correct.umi = TRUE)
TBM@meta.data[["branches"]] <- TBM$SingleR_labels
TBM@meta.data[["condition"]] <- TBM$identity
TBM <- TBM[, which(TBM$SingleR_labels %in% c("B-cells", "Granulocytes", "Retic", "Monocytes", "EryBl", "T", "NK", "DC"))]
TBM$branches <- factor(TBM$branches, levels = c("Granulocytes", "Monocytes", "DC", "B-cells", "T", "NK", "EryBl", "Retic"))
saveRDS(object=TBM, file="TBM_main_populations_20220804.rds")

Seurat_detection_TBM <- create_seurat_object_with_gene_detection_rates_in_MPP_subsets(MPP=TBM)
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_detection_TBM, gene_set=genes_mobility_in_TBM, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=30, plot_hight=15, 
                                                                   y_axis_label="% of cells with detectable mRNA", 
                                                                   file_name="gene_plots_patch_genes_mobility_TBM.pdf",
                                                                   mode="detection_frequency")
Seurat_average_expression_TBM <- create_seurat_object_with_average_expression_per_gene_in_MPP_subsets(MPP=TBM)
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_average_expression_TBM, gene_set=genes_mobility_in_TBM, 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=30, plot_hight=15, 
                                                                   y_axis_label="average of normalized mRNA counts", 
                                                                   file_name="gene_plots_patch_genes_mobility_average_corrected_counts_TBM.pdf",
                                                                   mode="gene_expression")

# plot Usp22 expression
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_average_expression_TBM, gene_set=c("Usp22"), 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=17, plot_hight=7.5, 
                                                                   y_axis_label="average of normalized mRNA counts", 
                                                                   file_name="Usp22_TBM.pdf",
                                                                   mode="gene_expression")


# show Usp22 expression profile across all populations
TBM <- readRDS(file="TBM_main_populations_20220804.rds")
LSK <- readRDS(file="LSK_all_experiments_normalized_together.rds")
DefaultAssay(TBM) <- "RNA"
DefaultAssay(LSK) <- "RNA"
LSK@assays$SCT <- NULL
all_cells <- merge(TBM, y = c(LSK), add.cell.ids = NULL, project = "U22")
rm(LSK)
rm(TBM)
saveRDS(object=all_cells, file="all_populations_all_replicates_20220804.rds")
DefaultAssay(all_cells) <- "RNA"
all_cells <- Seurat::SCTransform(all_cells, verbose = FALSE, residual.features = c("Usp22"), conserve.memory = TRUE, do.correct.umi = TRUE)
all_cells$branches[which(is.na(all_cells$branches))] <- all_cells$sample[which(is.na(all_cells$branches))]

Seurat_average_expression_all_cells <- create_seurat_object_with_average_expression_per_gene_in_MPP_subsets(MPP=all_cells)
compare_expression_of_selected_genes_in_MPP_between_WT_and_KO_mice(seurat_object=Seurat_average_expression_all_cells, gene_set=c("Usp22"), 
                                                                   folder_with_plots=folder_with_plots, 
                                                                   plot_width=25, plot_hight=7, 
                                                                   y_axis_label="average of normalized mRNA counts", 
                                                                   file_name="Usp22_all_cells.pdf",
                                                                   mode="gene_expression")









