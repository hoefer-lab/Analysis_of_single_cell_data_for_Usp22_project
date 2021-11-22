# helper functions for annotation
.read_ref_raw_counts <- function(){
  ref_data_immgen <- read.csv("GSE109125_Gene_count_table.csv", sep = ",", row.names = 1)
  ref_data_haemoshere <- read.csv("Haemopedia_Mouse_RNASeq_raw.txt", sep = "\t", row.names = 1)
  haemoshere_sample_table <- read.csv("Haemopedia_Mouse_RNASeq_samples.txt", sep = "\t", row.names = 1)
  ref_data_immgen2 <- read.csv("GSE122108_Gene_count_table.csv", sep = ",", row.names = 1)
  return(list(ref_data_immgen=ref_data_immgen, 
              ref_data_haemoshere=ref_data_haemoshere, 
              haemoshere_sample_table=haemoshere_sample_table, 
              ref_data_immgen2=ref_data_immgen2))
}
.clean_haemosphere_data <- function(haemoshere_sample_table, ref_data_haemoshere){
  haemosphere_colnames_bm <- rownames(haemoshere_sample_table)[which(grepl("BM", haemoshere_sample_table$tissue, fixed = TRUE))]
  haemosphere_colnames_untreated <- rownames(haemoshere_sample_table)[which(grepl("none", haemoshere_sample_table$treatment, fixed = TRUE))]
  haemosphere_colnames_use <- intersect(haemosphere_colnames_bm, haemosphere_colnames_untreated)
  ref_data_haemoshere <-  ref_data_haemoshere[,haemosphere_colnames_use]
  return(ref_data_haemoshere)
}
.translate_gene_ID <- function(gene_ID_vector){
  ens.mm.v101 <- AnnotationHub()[["AH83247"]]
  gene_symbols <- mapIds(ens.mm.v101, keys=gene_ID_vector,
                         keytype="GENEID", column="SYMBOL")
  gene_symbols_UNIQUEFEATURENAMES <- uniquifyFeatureNames(gene_ID_vector, 
                                                          gene_symbols)
  return(gene_symbols_UNIQUEFEATURENAMES)
}
.extract_T_NKT_cells_from_immgen_which_are_not_from_bone_marrow <- function(ref_data_immgen){
  ref_data_immgen_multiple_tissues <- ref_data_immgen
  toMatch <- c("NKT.", "preT.", "T.", "T8.", "Tgd.", "Treg.")
  columns_to_retain <- c()
  for (index in (1:length(toMatch))){
    expression_to_retain <- toMatch[index]
    col_add <- which(grepl(expression_to_retain, colnames(ref_data_immgen_multiple_tissues), fixed = TRUE))
    columns_to_retain <- append(columns_to_retain, col_add)
  }
  columns_to_retain <- unique(columns_to_retain)
  ref_data_immgen_multiple_tissues <- ref_data_immgen_multiple_tissues[, columns_to_retain]
  ref_data_immgen_multiple_tissues <- ref_data_immgen_multiple_tissues[, c(1:17, 18:24, 25:32, 33:36, 37:38, 39:40, 47:59, 60:63)]
  return(ref_data_immgen_multiple_tissues)
}
.create_same_feature_space_for_all_samples <- function(object_to_annotate, ref_data_immgen, ref_data_haemoshere, ref_data_immgen2, ref_data_immgen_multiple_tissues, assay_to_use){
  common_genes <- rownames(object_to_annotate@assays[[assay_to_use]]@counts)[which((rownames(object_to_annotate@assays[[assay_to_use]]@counts) %in% rownames(ref_data_immgen)) 
                                                                                   & (rownames(object_to_annotate@assays[[assay_to_use]]@counts) %in% rownames(ref_data_haemoshere))
                                                                                   & (rownames(object_to_annotate@assays[[assay_to_use]]@counts) %in% rownames(ref_data_immgen2))
                                                                                   & (rownames(object_to_annotate@assays[[assay_to_use]]@counts) %in% rownames(ref_data_immgen_multiple_tissues)))]
  sc_data <- object_to_annotate@assays[[assay_to_use]]@counts[common_genes,]
  ref_data_immgen <- ref_data_immgen[common_genes,]
  ref_data_immgen2 <- ref_data_immgen2[common_genes,]
  ref_data_haemoshere <- ref_data_haemoshere[common_genes,]
  ref_data_immgen_multiple_tissues <- ref_data_immgen_multiple_tissues[common_genes,]
  return(list(sc_data=sc_data,
              ref_data_immgen=ref_data_immgen,
              ref_data_immgen2=ref_data_immgen2,
              ref_data_haemoshere=ref_data_haemoshere,
              ref_data_immgen_multiple_tissues=ref_data_immgen_multiple_tissues))
}
.TPM_normalization <- function(genes_for_normalization, 
                               ref_data_immgen, 
                               ref_data_haemoshere, 
                               ref_data_immgen2, 
                               ref_data_immgen_multiple_tissues){
  ens.mm.v101 <- AnnotationHub()[["AH83247"]]
  common_genes_ensemble <- mapIds(ens.mm.v101, keys=genes_for_normalization,
                                  keytype="SYMBOL", column="GENEID")
  library(EDASeq)
  gene_length <- EDASeq::getGeneLengthAndGCContent(
    id = common_genes_ensemble,
    org = 'mmusculus',
    mode = c('biomart', 'org.db'))
  gene_length_vector <- gene_length[common_genes_ensemble, "length"]
  for (column in (1:ncol(ref_data_immgen))){
    ref_data_immgen[, column] <- ref_data_immgen[, column]/gene_length_vector
  }
  for (column in (1:ncol(ref_data_immgen2))){
    ref_data_immgen2[, column] <- ref_data_immgen2[, column]/gene_length_vector
  }
  for (column in (1:ncol(ref_data_haemoshere))){
    ref_data_haemoshere[, column] <- ref_data_haemoshere[, column]/gene_length_vector
  }
  for (column in (1:ncol(ref_data_immgen_multiple_tissues))){
    ref_data_immgen_multiple_tissues[, column] <- ref_data_immgen_multiple_tissues[, column]/gene_length_vector
  }
  # as some gene length were not available, NA have to be removed
  library(tidyr)
  ref_data_immgen <- ref_data_immgen %>% drop_na()
  ref_data_immgen2 <- ref_data_immgen2 %>% drop_na()
  ref_data_haemoshere <- ref_data_haemoshere %>% drop_na()
  ref_data_immgen_multiple_tissues <- ref_data_immgen_multiple_tissues %>% drop_na()
  # scaling is done for each dataset below
  for (column in (1:ncol(ref_data_immgen))){
    scaling_factor <- 1000000/sum(ref_data_immgen[, column])
    for (element in (1:nrow(ref_data_immgen))){
      ref_data_immgen[element, column] <- ref_data_immgen[element, column]*scaling_factor
    }
  }
  
  for (column in (1:ncol(ref_data_immgen2))){
    scaling_factor <- 1000000/sum(ref_data_immgen2[, column])
    for (element in (1:nrow(ref_data_immgen2))){
      ref_data_immgen2[element, column] <- ref_data_immgen2[element, column]*scaling_factor
    }
  }
  
  for (column in (1:ncol(ref_data_haemoshere))){
    scaling_factor <- 1000000/sum(ref_data_haemoshere[, column])
    for (element in (1:nrow(ref_data_haemoshere))){
      ref_data_haemoshere[element, column] <- ref_data_haemoshere[element, column]*scaling_factor
    }
  }
  
  for (column in (1:ncol(ref_data_immgen_multiple_tissues))){
    scaling_factor <- 1000000/sum(ref_data_immgen_multiple_tissues[, column])
    for (element in (1:nrow(ref_data_immgen_multiple_tissues))){
      ref_data_immgen_multiple_tissues[element, column] <- ref_data_immgen_multiple_tissues[element, column]*scaling_factor
    }
  }
  return(list(ref_data_immgen=ref_data_immgen,
              ref_data_immgen2=ref_data_immgen2,
              ref_data_haemoshere=ref_data_haemoshere,
              ref_data_immgen_multiple_tissues=ref_data_immgen_multiple_tissues))
}
.rename_samples <- function(ref_data_immgen_log2,
                            ref_data_haemoshere_log2,
                            ref_data_immgen_multiple_tissues_log2){
  colnames(ref_data_immgen_log2) <- c("B","B","B","GN","GN","LTHSC","LTHSC",
                                      "LTHSC","LTHSC","MMP2","MMP3","MMP3","MMP4","MMP4",
                                      "NK","NK","NK","NK","NK","NK","proB",
                                      "proB","proB","proB","proB","proB","STHSC","STHSC", 
                                      "DC","DC","DC","DC","DC","DC",
                                      "DC","DC","DC")
  colnames(ref_data_haemoshere_log2) <- c("Baso","Baso","Eo","Eo","Eo","Eo","Eo","Eo","Eo","Eo","Eo","Eo",
                                          "Eo","Eo","Eo","Eo","Eo","Eo","PreCFUE","PreCFUE","MEP","MEP","CFUE","CFUE",
                                          "EryBl","EryBl","EryBl","EryBl","EryBl","EryBl","Retic","Retic","Retic","Mono","Mono","InfMono",
                                          "InfMono","STHSC","STHSC","LSK","LSK","MPP","MPP","Neut","Neut","CMP","CMP","GMP",
                                          "GMP","GMP","GMP","GMP","GMP","GMP","CLP","CLP")
  colnames(ref_data_immgen_multiple_tissues_log2) <- c(rep("NKT", 17), rep("preT", 7), rep("T4", 8), rep("T8", 4), rep("T.DN", 2), rep("T.DP", 2), rep("T8", 11), rep("Tgd", 6))
  
  return(list(ref_data_immgen_log2=ref_data_immgen_log2,
              ref_data_haemoshere_log2=ref_data_haemoshere_log2,
              ref_data_immgen_multiple_tissues_log2=ref_data_immgen_multiple_tissues_log2))
}