source("script2.R", local=FALSE)
path_to_working_directory <- getwd()
folder_with_plots <- "branch_specific_DEG_analysis_all_mice_20220811"
dir.create(paste(path_to_working_directory, folder_with_plots, sep="/"), showWarnings = FALSE)

MPP <- readRDS(file="MPP_all_experiments")
DefaultAssay(MPP) <- "RNA"
MPP <- Seurat::SCTransform(MPP, verbose = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3, conserve.memory = TRUE, do.correct.umi = TRUE)

# run DEG analysis for each branch in each experiment separately ####
experiment_wise_results <- list(experiment1=list(), experiment2=list())
for (index in (1:length(unique(MPP$experiment)))){
  current_exp <- unique(MPP$experiment)[index]
  for (index2 in (1:length(unique(MPP$branches)))){
    current_branch <- unique(MPP$branches)[index2]
    MPP_current_exp_branch <- MPP[, which((MPP$experiment==current_exp) & (MPP$branches==current_branch))]
    DEG_results <- Seurat::FindMarkers(MPP_current_exp_branch, group.by="condition", ident.1="KO", ident.2="WT", 
                                       test.use="wilcox", assay="SCT", slot="data", logfc.threshold = 0.25, pseudocount.use=0.1)
    significant_positive_genes <- rownames(DEG_results)[which((DEG_results$p_val_adj<0.05) & (DEG_results$avg_log2FC>0.25))]
    significant_negative_genes <- rownames(DEG_results)[which((DEG_results$p_val_adj<0.05) & (DEG_results$avg_log2FC<(-0.25)))]
    experiment_wise_results[[index]] <- list.append(experiment_wise_results[[index]], list())
    experiment_wise_results[[index]][[index2]] <- list.append(experiment_wise_results[[index]][[index2]], significant_positive_genes)
    experiment_wise_results[[index]][[index2]] <- list.append(experiment_wise_results[[index]][[index2]], significant_negative_genes)
    names(experiment_wise_results[[index]][[index2]]) <- c("positive", "negative")
  }
  names(experiment_wise_results[[index]]) <- unique(MPP$branches)
}


# for each branch, identify positively and negatively regulated genes that are differentially expressed in both experiments ####
#experiments <- names(experiment_wise_results)
branches <- unique(MPP$branches)
branch_wise_DEGs_shared_between_experiments <- list()
for (index in (1:length(branches))){
  current_branch <- branches[index]
  shared_positive_DEGs <- intersect(experiment_wise_results[[1]][[current_branch]][["positive"]], 
                                    experiment_wise_results[[2]][[current_branch]][["positive"]])
  shared_negative_DEGs <- intersect(experiment_wise_results[[1]][[current_branch]][["negative"]], 
                                    experiment_wise_results[[2]][[current_branch]][["negative"]])
  branch_wise_DEGs_shared_between_experiments <- list.append(branch_wise_DEGs_shared_between_experiments, list())
  branch_wise_DEGs_shared_between_experiments[[index]] <- list.append(branch_wise_DEGs_shared_between_experiments[[index]],
                                                                      shared_positive_DEGs)
  branch_wise_DEGs_shared_between_experiments[[index]] <- list.append(branch_wise_DEGs_shared_between_experiments[[index]],
                                                                      shared_negative_DEGs)
  names(branch_wise_DEGs_shared_between_experiments[[index]]) <- c("positive", "negative")
}
names(branch_wise_DEGs_shared_between_experiments) <- branches


# create heatmap of DEGs for each branch; show values of each experiment ####
for (index in (1:length(branches))){
  current_branch <- branches[index]
  data_exp1 <- MPP[, which((MPP$experiment=="1") & (MPP$branches==current_branch))]
  data_exp2 <- MPP[, which((MPP$experiment=="2") & (MPP$branches==current_branch))]
  
  genes_to_plot <- c(branch_wise_DEGs_shared_between_experiments[[current_branch]]$positive, branch_wise_DEGs_shared_between_experiments[[current_branch]]$negative)
  
  DEG_results1 <- Seurat::FindMarkers(data_exp1, features=genes_to_plot, group.by="condition", ident.1="KO", ident.2="WT", 
                                     test.use="wilcox", assay="SCT", slot="data", logfc.threshold = 0.0, pseudocount.use=1, min.pct=0.0)
  
  DEG_results2 <- Seurat::FindMarkers(data_exp2, features=genes_to_plot, group.by="condition", ident.1="KO", ident.2="WT", 
                                      test.use="wilcox", assay="SCT", slot="data", logfc.threshold = 0.0, pseudocount.use=1, min.pct=0.0)
  
  exp1_log2fcs <- DEG_results1[genes_to_plot, "avg_log2FC"]
  exp2_log2fcs <- DEG_results2[genes_to_plot, "avg_log2FC"]
  
  combined_fold_changes <- data.frame(experiment_1=exp1_log2fcs, experiment_2=exp2_log2fcs)
  rownames(combined_fold_changes) <- genes_to_plot
  
  mean_values <- apply(combined_fold_changes, 1, mean)
  combined_fold_changes_sorted <- combined_fold_changes[order(mean_values, decreasing=TRUE),]
  
  figure_height <- round(length(genes_to_plot))/10
  
  color <- colorRampPalette(c("blue", "white", "red"))(100)
  myBreaks <- c(seq(min(combined_fold_changes_sorted, na.rm = TRUE), 0, length.out=ceiling(100/2) + 1), 
                seq(max(combined_fold_changes_sorted, na.rm = TRUE)/100, max(combined_fold_changes_sorted, na.rm = TRUE), length.out=floor(100/2)))
  pheatmap::pheatmap(combined_fold_changes_sorted, na_col = "grey", color = color, breaks = myBreaks, cluster_cols = FALSE, cluster_rows = FALSE, 
                     angle_col=270, fontsize_row=6, fontsize_col=8, fontsize=8, 
                     main = "", 
                     filename = paste(folder_with_plots, paste0(paste("DEG_heatmap_all_replicates", current_branch, sep = "_"), ".pdf"), sep="/"), width = 3, height = figure_height)
  
  
}

# identify genes which are differentially expressed in all branches in both experiments ####
positive_intersection <- intersect(intersect(intersect(branch_wise_DEGs_shared_between_experiments$undefined$positive, branch_wise_DEGs_shared_between_experiments$CLP$positive), branch_wise_DEGs_shared_between_experiments$GMP$positive), branch_wise_DEGs_shared_between_experiments$MEP$positive)
negative_intersection <- intersect(intersect(intersect(branch_wise_DEGs_shared_between_experiments$undefined$negative, branch_wise_DEGs_shared_between_experiments$CLP$negative), branch_wise_DEGs_shared_between_experiments$GMP$negative), branch_wise_DEGs_shared_between_experiments$MEP$negative)

# compute experiment-wise fold changes for DEGs that are shared between all branches for pooled MPPs
experiment_wise_results_pooled <- list(experiment1=list(), experiment2=list())
for (index in (1:length(unique(MPP$experiment)))){
  current_exp <- unique(MPP$experiment)[index]
  MPP_current_exp <- MPP[, which(MPP$experiment==current_exp)]
  DEG_results <- Seurat::FindMarkers(MPP_current_exp, features=c(positive_intersection, negative_intersection) , group.by="condition", ident.1="KO", ident.2="WT", 
                                     test.use="wilcox", assay="SCT", slot="data", logfc.threshold = 0.0, pseudocount.use=1, min.pct=0.0)
  experiment_wise_results_pooled[[index]] <- DEG_results[c(positive_intersection, negative_intersection), "avg_log2FC"]
}

# visualize DEGs that are shared between branches in heatmap 
shared_DEGs_table <- data.frame(experiment_1=experiment_wise_results_pooled[[1]],
                                experiment_2=experiment_wise_results_pooled[[2]])
rownames(shared_DEGs_table) <- c(positive_intersection, negative_intersection)
row_means <- apply(shared_DEGs_table, 1, mean)
shared_DEGs_table_sorted <- shared_DEGs_table[order(row_means, decreasing = TRUE),]
figure_height <- (round(nrow(shared_DEGs_table_sorted))/10)+1
color <- colorRampPalette(c("blue", "white", "red"))(100)
myBreaks <- c(seq(min(shared_DEGs_table_sorted, na.rm = TRUE), 0, length.out=ceiling(100/2) + 1), 
              seq(max(shared_DEGs_table_sorted, na.rm = TRUE)/100, max(shared_DEGs_table_sorted, na.rm = TRUE), length.out=floor(100/2)))
pheatmap::pheatmap(shared_DEGs_table_sorted, na_col = "grey", color = color, breaks = myBreaks, cluster_cols = FALSE, cluster_rows = FALSE, 
                   angle_col=270, fontsize_row=6, fontsize_col=8, fontsize=8, 
                   main = "", 
                   filename = paste(folder_with_plots, paste0("shared_DEGs_pooled_MPPs", ".pdf"), sep="/"), width = 1.5, height = figure_height)
