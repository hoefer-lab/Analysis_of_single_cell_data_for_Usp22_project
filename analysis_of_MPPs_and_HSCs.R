# load MPPs, annotate them by SingleR and prepare them for FateID ####
set.seed(1)
sample(5)
MPP_WT <- load_process_data(path=path_MPP_WT, nFeature_lower_threshold=2500, nFeature_upper_threshold=10000, percent.mt_threshold=6, including_normalization=TRUE)
MPP_KO <- load_process_data(path=path_MPP_KO, nFeature_lower_threshold=2500, nFeature_upper_threshold=10000, percent.mt_threshold=6, including_normalization=TRUE)

# annotate MPPs by SingleR
MPP_WT <- annotate_by_SingleR(object_to_annotate=MPP_WT, path_helper_functions="annotation_helper_functions.R")
MPP_KO <- annotate_by_SingleR(object_to_annotate=MPP_KO, path_helper_functions="annotation_helper_functions.R")
# In case that the biomart server cannot be reached, annotated objects can be loaded:
# saveRDS(object=MPP_WT, file="/Volumes/addition_storage/R_objects_Usp22/MPP_WT_annotated_20220311_2.rds")
# saveRDS(object=MPP_KO, file="/Volumes/addition_storage/R_objects_Usp22/MPP_KO_annotated_20220311_2.rds")


# assign cell-cycle phases to cells of both conditions independently and regress out cell-cycle signal from data
prepared_MPPs <- prepare_MPPs_for_FateID(MPP_WT=MPP_WT, MPP_KO=MPP_KO, g2m.genes=g2m.genes, s.genes=s.genes)
MPP_WT <- prepared_MPPs[["MPP_WT"]]
MPP_KO <- prepared_MPPs[["MPP_KO"]]
rm(prepared_MPPs)

# feature selection for MPPs by FateID from both conditions independently ####
common_features <- select_lineage_specific_genes_by_FateID(MPP_WT, MPP_KO)
# the select_lineage_specific_genes_by_FateID function is very slow; to save time, features can also be loaded from file: features_from_FateID_20220311

# merge WT and KO and redo normalization on combined objects ####
MPP_combined <- merge(MPP_WT, y = MPP_KO, add.cell.ids = NULL, project = "U22")
MPP_combined@meta.data[["condition"]] <- c(rep("WT", ncol(MPP_WT)), rep("KO", ncol(MPP_KO)))
DefaultAssay(MPP_combined) <- "RNA"
MPP_combined <- Seurat::SCTransform(MPP_combined, verbose = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)
rm(MPP_WT)
rm(MPP_KO)

# compute diffusion map of MPPs based on normalized expression values of marker genes identified by FateID ####
dmap_results <- compute_diffusion_map(MPP_combined, common_features)
MPP_combined <- dmap_results[["MPP_combined"]]
MPP_combined.sce <- dmap_results[["MPP_combined.sce"]]
rm(dmap_results)

# assign MPPs to branches using Slingshot trajectories ####
objects_with_classified_cells <- assign_cells_to_branches_using_slingshot(MPP_combined, MPP_combined.sce)
MPP_combined <- objects_with_classified_cells[["MPP_combined"]]
MPP_combined.sce <- objects_with_classified_cells[["MPP_combined.sce"]]
rm(objects_with_classified_cells)
# compute percentage of cells in each MPP subset for both conditions ####
LMPP_WT <- length(which((MPP_combined$branches=="CLP") & (MPP_combined$condition=="WT"))) / length(which(MPP_combined$condition=="WT"))
GMP_biased_MPP_WT <- length(which((MPP_combined$branches=="GMP") & (MPP_combined$condition=="WT"))) / length(which(MPP_combined$condition=="WT"))
MEP_biased_MPP_WT <- length(which((MPP_combined$branches=="MEP") & (MPP_combined$condition=="WT"))) / length(which(MPP_combined$condition=="WT"))
unbiased_MPP_WT <- length(which((MPP_combined$branches=="undefined") & (MPP_combined$condition=="WT"))) / length(which(MPP_combined$condition=="WT"))

LMPP_KO <- length(which((MPP_combined$branches=="CLP") & (MPP_combined$condition=="KO"))) / length(which(MPP_combined$condition=="KO"))
GMP_biased_MPP_KO <- length(which((MPP_combined$branches=="GMP") & (MPP_combined$condition=="KO"))) / length(which(MPP_combined$condition=="KO"))
MEP_biased_MPP_KO <- length(which((MPP_combined$branches=="MEP") & (MPP_combined$condition=="KO"))) / length(which(MPP_combined$condition=="KO"))
unbiased_MPP_KO <- length(which((MPP_combined$branches=="undefined") & (MPP_combined$condition=="KO"))) / length(which(MPP_combined$condition=="KO"))


# inspect Slingshot trajectories and plot embedding for manuscript ####
create_plots_showing_MPP_embedding_with_curves(MPP_combined.sce=MPP_combined.sce, path_to_plot_parameters="3d_scatter_parameters_20220312.rds", path_to_plots=folder_with_plots)
create_plots_showing_MPP_embedding(MPP_combined.sce=MPP_combined.sce, path_to_plot_parameters="3d_scatter_parameters_20220312.rds", path_to_plots=folder_with_plots)

# cell cycle analysis of MPPs ####
cell_cycle_results_MPPs <- assign_cells_to_cell_cycle_phases_based_on_gene_lists(object_to_annotate=MPP_combined,
                                                                                 G2M_thres_in_neighborhood=1,
                                                                                 neighborhood_size=200,
                                                                                 g2m.genes=g2m.genes,
                                                                                 s.genes=s.genes)
MPP_combined <- cell_cycle_results_MPPs[["annotated_object"]]
cell_cycle_results_MPPs[["umap_cell_cycle"]]
rm(cell_cycle_results_MPPs)

plot_cell_cycle_changes_as_heatmap(annotated_object=MPP_combined, path_to_plots=folder_with_plots)
plot_cell_cycle_changes_in_barplot_MPPs(annotated_object=MPP_combined, path_to_plots=folder_with_plots, width = 10, height = 8)

# cell-cycle analysis of LT-HSCs ####
LTHSC_WT <- load_process_data(path=path_HSC_WT, nFeature_lower_threshold=2500, nFeature_upper_threshold=10000, percent.mt_threshold=6, including_normalization=FALSE)
LTHSC_KO <- load_process_data(path=path_HSC_KO, nFeature_lower_threshold=2500, nFeature_upper_threshold=10000, percent.mt_threshold=6, including_normalization=FALSE)
LTHSC_combined <- merge(LTHSC_WT, y = c(LTHSC_KO), add.cell.ids = NULL, project = "U22")
LTHSC_combined@meta.data[["condition"]] <- c(rep("WT", ncol(LTHSC_WT)),
                                             rep("KO", ncol(LTHSC_KO)))
rm(LTHSC_WT)
rm(LTHSC_KO)
DefaultAssay(LTHSC_combined) <- "RNA"
LTHSC_combined <- Seurat::SCTransform(LTHSC_combined, verbose = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)
cell_cycle_results_LTs <- assign_cells_to_cell_cycle_phases_based_on_gene_lists(object_to_annotate=LTHSC_combined,
                                                                                G2M_thres_in_neighborhood=150,
                                                                                neighborhood_size=200,
                                                                                g2m.genes=g2m.genes,
                                                                                s.genes=s.genes)
LTHSC_combined <- cell_cycle_results_LTs[["annotated_object"]]
plot_umaps_cell_cycle_updated(WT=LTHSC_combined[,which(LTHSC_combined$condition=="WT")], 
                              KO=LTHSC_combined[,which(LTHSC_combined$condition=="KO")], 
                              path_for_plot=folder_with_plots)
create_barplot_with_percentage_of_cycling_LTs(LTHSC_combined=LTHSC_combined, path_to_plots=folder_with_plots)
plot_cell_cycle_changes_in_barplot_LTs(annotated_object=LTHSC_combined, path_to_plots=folder_with_plots, width = 7, height = 8)
visualize_cell_cycle_gene_expression_in_LTHSCs(LTHSC_combined=LTHSC_combined, folder_with_plots=folder_with_plots)
plot_cell_cycle_changes_in_barplot_LTs(annotated_object=LTHSC_combined, path_to_plots=folder_with_plots, width = 6.5, height = 8)
rm(cell_cycle_results_LTs)

# create sampled version of objects in which densities of WT and KO cells are comparable throughout the embedding ####
sampled_embeddings <- correct_MPP_embedding_for_differential_ditribution_of_conditions(MPP_combined=MPP_combined, MPP_combined.sce=MPP_combined.sce, numbers_of_direct_neighbors_to_consider=1)
MPP_combined.corrected_cell_numbers <- sampled_embeddings[["MPP_combined.corrected_cell_numbers"]]
MPP_combined.sce.corrected_cell_numbers <- sampled_embeddings[["MPP_combined.sce.corrected_cell_numbers"]]
sampled_embeddings[["scatterplot_sampled_object"]]
sampled_embeddings[["max_sampling_number"]]
sampled_embeddings[["percentage_of_KO_cells_sampled_once"]]
sampled_embeddings[["percentage_of_KO_cells_sampled_twice"]]
rm(sampled_embeddings)

# test if KO and WT cells are equally distributed along trajectories after sampling ####
MPP_combined.sce.corrected_cell_numbers@colData@listData[["Sample"]] <- rep("MPP", ncol(MPP_combined.sce.corrected_cell_numbers))
MPP_combined.sce.corrected_cell_numbers@colData@listData[["identity"]] <- MPP_combined.sce.corrected_cell_numbers$condition

list_with_trajectories_and_their_medians_and_CIs_ordered_by_pseudotime <- run_condition_imbalance_analysis(Interesting_cells=MPP_combined.sce.corrected_cell_numbers, sample_to_use="MPP",
                                                                                                           pseudo_time_trajectories_to_use=c("GMP_pseudotime_na", "CLP_pseudotime_na", "MEP_pseudotime_na"),
                                                                                                           size_of_local_neighborhoods=50, num_bootstrap_iter=10,
                                                                                                           path_to_helper_functions="condition_imbalance_analysis_helper_functions.R",
                                                                                                           size_moving_window=0.05)

df_confidence <- list_with_trajectories_and_their_medians_and_CIs_ordered_by_pseudotime[[1]]
plot_condition_imbalance_analysis(df_confidence, x_label="pseudo-time from HSCs to GMPs", y_label="KO/(KO+WT)", path_for_plot=paste(folder_with_plots, "condition_imbalance_analysis_myeloid_biased_MPPs_sampled.pdf", sep="/"))
df_confidence <- list_with_trajectories_and_their_medians_and_CIs_ordered_by_pseudotime[[2]]
plot_condition_imbalance_analysis(df_confidence, x_label="pseudo-time from HSCs to CLPs", y_label="KO/(KO+WT)", path_for_plot=paste(folder_with_plots, "condition_imbalance_analysis_LMPPs_sampled.pdf", sep="/"))
df_confidence <- list_with_trajectories_and_their_medians_and_CIs_ordered_by_pseudotime[[3]]
plot_condition_imbalance_analysis(df_confidence, x_label="pseudo-time from HSCs to MEPs", y_label="KO/(KO+WT)", path_for_plot=paste(folder_with_plots, "condition_imbalance_analysis_MEP_biased_MPPs_sampled.pdf", sep="/"))
rm(list_with_trajectories_and_their_medians_and_CIs_ordered_by_pseudotime)
rm(df_confidence)

# compute scores of signatures of GMP-biased MPPs and LMPPs for each cell in MPPs, LTs and STs; visualize as violin plot ####
LTHSC_WT <- load_process_data(path=path_HSC_WT, nFeature_lower_threshold=2500, nFeature_upper_threshold=10000, percent.mt_threshold=6, including_normalization=FALSE)
LTHSC_KO <- load_process_data(path=path_HSC_KO, nFeature_lower_threshold=2500, nFeature_upper_threshold=10000, percent.mt_threshold=6, including_normalization=FALSE)
STHSC_WT <- load_process_data(path=path_ST_WT, nFeature_lower_threshold=2500, nFeature_upper_threshold=10000, percent.mt_threshold=6, including_normalization=FALSE)
STHSC_KO <- load_process_data(path=path_ST_KO, nFeature_lower_threshold=2500, nFeature_upper_threshold=10000, percent.mt_threshold=6, including_normalization=FALSE)

HSC_combined <- merge(LTHSC_WT, y = c(LTHSC_KO, STHSC_WT, STHSC_KO), add.cell.ids = NULL, project = "U22")
HSC_combined@meta.data[["condition"]] <- c(rep("WT", ncol(LTHSC_WT)),
                                           rep("KO", ncol(LTHSC_KO)),
                                           rep("WT", ncol(STHSC_WT)),
                                           rep("KO", ncol(STHSC_KO)))
HSC_combined@meta.data[["sample"]] <- c(rep("LT-HSCs", (ncol(LTHSC_WT) + ncol(LTHSC_KO))),
                                        rep("ST-HSCs", (ncol(STHSC_WT) + ncol(STHSC_KO))))
rm(LTHSC_WT)
rm(LTHSC_KO)
rm(STHSC_WT)
rm(STHSC_KO)
DefaultAssay(HSC_combined) <- "RNA"
HSC_combined <- Seurat::SCTransform(HSC_combined, verbose = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)

signatures <- compute_per_cell_scores_of_GMP_LMPP_priming_and_plot_distribution(object_for_signature=MPP_combined,
                                                                  objects_for_scoring=list(MPP_combined.corrected_cell_numbers=MPP_combined.corrected_cell_numbers,
                                                                                           HSC_combined=HSC_combined),
                                                                  path_for_plot= paste(folder_with_plots, "/",sep=""))

# create heat maps showing differential expression of signature genes in all HSPC subsets ####
plot_differential_expression_of_signatures(signature_genes=signatures[["GMP_signature"]], MPP_combined.corrected_cell_numbers=MPP_combined.corrected_cell_numbers,
                                           HSC_combined=HSC_combined,
                                           plot_title="Differential expression\nof genes specific\nfor GMP-biased MPPs",
                                           full_path=paste(folder_with_plots, "GMP_signature_heatmap.pdf", sep="/"))

plot_differential_expression_of_signatures(signature_genes=signatures[["CLP_signature"]], MPP_combined.corrected_cell_numbers=MPP_combined.corrected_cell_numbers,
                                           HSC_combined=HSC_combined,
                                           plot_title="Differential expression\nof genes specific\nfor LMPPs",
                                           full_path=paste(folder_with_plots, "LMPP_signature_heatmap.pdf", sep="/"))
rm(HSC_combined)

# do GSEA for MPPs ####
processed_GSEA_results <- run_GSEA_and_group_significant_GO_terms_by_overlap_of_gene_sets(seurat_object=MPP_combined.corrected_cell_numbers, significance_threshold=0.05, edge_threshold_Jaccard=0.3, min_number_of_leading_edge_genes=70)
# after inspection, the groups of GO terms need to be annotated manually:
names(processed_GSEA_results[["GO_term_groups_cut_global"]]) <- c("cellular protein-containing complex assembly",
                                                                  "mitochondrion organization",
                                                                  "response to DNA damage stimulus",
                                                                  "response to growth factor, receptor signaling pathway",
                                                                  "exocytosis, endocytosis, vesicle-mediated transport",
                                                                  "nucleotide metabolic process",
                                                                  "aerobic respiration",
                                                                  "RNA processing",
                                                                  "DNA replication",
                                                                  "autophagy",
                                                                  "movement of cell or subcellular component",
                                                                  "immune processes",
                                                                  "tissue development"
                                                                  )
plot_filtered_GO_network(processed_GSEA_results=processed_GSEA_results, min_edge=0.3, path_to_plots=folder_with_plots, sample_name = "MPPs", use_old_coordinates=TRUE, save_new_coordinates=FALSE, file_name_new_coordinates=NULL, file_name_old_coordinates="MPPs_emapplot_coordinates_20220316")
# plot heatmap at level of global clusters
plot_GSEA_heatmap_global(processed_GSEA_results=processed_GSEA_results, log2FC_threshold=1.2, path_to_file=folder_with_plots, sample_name = "MPPs")

# heatmap immune subclusters
subclustering_results <- subcluster_GO_group(results_readable=processed_GSEA_results[["results_readable"]], mother_cluster="immune processes", GO_term_groups_cut_global=processed_GSEA_results[["GO_term_groups_cut_global"]])
kelley_results <- kgs(cluster=subclustering_results[["hclust_ward"]], diss=subclustering_results[["JC_dist"]], alpha=1, maxclust=NULL)
plot(names(kelley_results), kelley_results, xlab="# clusters", ylab="penalty")
suggested_number_of_clusters <- as.integer(names(kelley_results)[which(kelley_results==min(kelley_results))])
subclusters_and_visualization <- visualize_subclustering_results(number_of_clusters=suggested_number_of_clusters, subclustering_results=subclustering_results, processed_GSEA_results=processed_GSEA_results, min_edge=0.3)
subclusters_and_visualization[["grouped_GO_terms_plot"]]
subclusters_and_visualization[["list_of_clustered_immune_GO_terms"]]
# based on clustering results, decide which short names are most appropriate for each group
names_of_subclusters <- c("antigen processing and presentation",
                          "regulation of immune responses (mainly adaptive)",
                          "regulation of cytotoxicity; antigen processing and presentation",
                          "innate immune response; response to bacterium / other organism",
                          "leukocyte activation / proliferation; response to IFN-gamma; regulation of inflammation",
                          "production of molecular mediator / cytokine",
                          "survival; fatty acid/prostaglandin metabolism")
plot_GSEA_heatmap_subcluster(processed_GSEA_results=processed_GSEA_results,
                             list_of_clustered_immune_GO_terms=subclusters_and_visualization[["list_of_clustered_immune_GO_terms"]],
                             names_of_subclusters=names_of_subclusters, mother_cluster="immune processes", log2FC_threshold=1.2,
                             number_of_clusters=suggested_number_of_clusters,
                             path_to_plot=folder_with_plots)




# do GSEA for of LT-HSCs ####
processed_GSEA_results_LTHSC <- run_GSEA_and_group_significant_GO_terms_by_overlap_of_gene_sets(seurat_object=LTHSC_combined, significance_threshold=0.05, edge_threshold_Jaccard=0.3, min_number_of_leading_edge_genes=100)
# after inspection, the groups of GO terms need to be annotated manually:
names(processed_GSEA_results_LTHSC[["GO_term_groups_cut_global"]]) <- c("cellular protein-containing complex assembly",
                                                                        "DNA replication",
                                                                        "development",
                                                                        "response to growth factor",
                                                                        "aerobic respiration; nucleotide metabolic processes",
                                                                        "RNA processing",
                                                                        "cell division",
                                                                        "cellular response to DNA damage stimulus",
                                                                        "mitochondrion organization",
                                                                        "small molecule biosynthetic process",
                                                                        "receptor internalization; endocytosis; vesicle-mediated transport")
plot_filtered_GO_network(processed_GSEA_results=processed_GSEA_results_LTHSC, min_edge=0.3, path_to_plots=folder_with_plots, sample_name = "LTs")
# plot heatmap at level of global clusters
plot_GSEA_heatmap_global(processed_GSEA_results=processed_GSEA_results_LTHSC, log2FC_threshold=1.3, path_to_file=folder_with_plots, sample_name = "LTs", figure_height=12, rowlabel_size=6, collabel_size=6)



# clean work space
elements_to_remove <- setdiff(ls(), lsf.str())
elements_to_remove <- c(elements_to_remove, "elements_to_remove")
elements_to_remove <- elements_to_remove[which(!(elements_to_remove %in% c("path_HSC_KO", "path_HSC_WT", "path_Linminus_KO", 
                                                                           "path_Linminus_WT", "path_LSK_KO", "path_LSK_WT",
                                                                           "path_MPP_KO", "path_MPP_WT", "path_ST_KO", 
                                                                           "path_ST_WT", "path_TBM_KO", "path_TBM_WT",
                                                                           "folder_with_plots")))]
rm(list = elements_to_remove)
gc()

# visualize expression of Usp22 in LT-HSC, ST-HSC and MPP ####
LTHSC_WT <- load_process_data(path=path_HSC_WT, nFeature_lower_threshold=2500, nFeature_upper_threshold=10000, percent.mt_threshold=6, including_normalization=TRUE)
LTHSC_KO <- load_process_data(path=path_HSC_KO, nFeature_lower_threshold=2500, nFeature_upper_threshold=10000, percent.mt_threshold=6, including_normalization=TRUE)
STHSC_WT <- load_process_data(path=path_ST_WT, nFeature_lower_threshold=2500, nFeature_upper_threshold=10000, percent.mt_threshold=6, including_normalization=TRUE)
STHSC_KO <- load_process_data(path=path_ST_KO, nFeature_lower_threshold=2500, nFeature_upper_threshold=10000, percent.mt_threshold=6, including_normalization=TRUE)
MPP_WT <- load_process_data(path=path_MPP_WT, nFeature_lower_threshold=2500, nFeature_upper_threshold=10000, percent.mt_threshold=6, including_normalization=TRUE)
MPP_KO <- load_process_data(path=path_MPP_KO, nFeature_lower_threshold=2500, nFeature_upper_threshold=10000, percent.mt_threshold=6, including_normalization=TRUE)

LTHSC_WT <- Seurat::RunUMAP(LTHSC_WT, dims = 1:15, return.model = TRUE, n.neighbors = 25, min.dist = 0.7, metric="cosine")#25 0.7
STHSC_WT <- Seurat::RunUMAP(STHSC_WT, dims = 1:15, return.model = TRUE, n.neighbors = 30, min.dist = 0.5, metric="cosine")#25 0.7
MPP_WT <- Seurat::RunUMAP(MPP_WT, dims = 1:15, return.model = TRUE, n.neighbors = 30, min.dist = 0.7, metric="cosine")#25 0.7

# for each WT/KO pair: align KO-sample to WT-sample in PCA space and use the WT-umap-model to compute umaps of both conditions
LTHSC_KO <- map_sample(ref_object=LTHSC_WT, query_object=LTHSC_KO, ndims=15, reference_reduction="pca", normalization_method="SCT", label_to_transfer=NULL, reduction_model="umap")
STHSC_KO <- map_sample(ref_object=STHSC_WT, query_object=STHSC_KO, ndims=15, reference_reduction="pca", normalization_method="SCT", label_to_transfer=NULL, reduction_model="umap")
MPP_KO <- map_sample(ref_object=MPP_WT, query_object=MPP_KO, ndims=15, reference_reduction="pca", normalization_method="SCT", label_to_transfer=NULL, reduction_model="umap")

# visualize Usp22 expression
create_usp22_panel_for_all_stem_cell_populations(HSC_WT=LTHSC_WT, ST_WT=STHSC_WT, MPP_WT=MPP_WT, HSC_KO=LTHSC_KO, ST_KO=STHSC_KO, MPP_KO=MPP_KO, path_for_plot=paste(folder_with_plots, "Usp22_expression_in_stem_cells.pdf", sep="/"))

# clean workspace
elements_to_remove <- setdiff(ls(), lsf.str())
elements_to_remove <- c(elements_to_remove, "elements_to_remove")
elements_to_remove <- elements_to_remove[which(!(elements_to_remove %in% c("path_HSC_KO", "path_HSC_WT", "path_Linminus_KO", 
                                                                           "path_Linminus_WT", "path_LSK_KO", "path_LSK_WT",
                                                                           "path_MPP_KO", "path_MPP_WT", "path_ST_KO", 
                                                                           "path_ST_WT", "path_TBM_KO", "path_TBM_WT",
                                                                           "folder_with_plots")))]
rm(list = elements_to_remove)
gc()



