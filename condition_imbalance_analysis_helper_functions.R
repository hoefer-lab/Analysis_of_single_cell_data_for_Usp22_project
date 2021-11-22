.create_matrix_with_number_of_observations_per_cell___bootstrap_iterations_correspond_to_rows <- function(seurat_sample, Indices_WT, Indices_KO, sample_size=min(c(length(Indices_WT), length(Indices_KO))), num_bootstrap_iter){
  matrix_samples <- matrix(0, nrow = num_bootstrap_iter, ncol = ncol(seurat_sample))
  colnames(matrix_samples) <- as.character(c(1:ncol(matrix_samples)))
  for (index in (1:num_bootstrap_iter)){
    current_vector_WT <- sample(Indices_WT, size = sample_size, replace = TRUE)
    current_vector_KO <- sample(Indices_KO, size = sample_size, replace = TRUE)
    matrix_samples[index, names(table(current_vector_WT))] <- as.vector(table(current_vector_WT))
    matrix_samples[index, names(table(current_vector_KO))] <- as.vector(table(current_vector_KO))
  }
  return(matrix_samples)
}

.compute_matrix_with_proportion_of_WT_cells_in_neighborhood___bootstrap_iterations_correspond_to_rows <- function(seurat_sample, matrix_samples, num_bootstrap_iter, nn_lists_in_rows, size_of_local_neighborhoods){
  matrix_with_scores <- matrix(0, ncol=ncol(matrix_samples), nrow=nrow(matrix_samples))
  for (index in (1:num_bootstrap_iter)){
    indices_of_available_cells <- which(matrix_samples[index,]!=0)
    for (inner_index in (1:ncol(matrix_samples))){
      cell <- inner_index
      counts <- matrix_samples[index, cell]
      if (counts!=0){
        potential_neighbors <- nn_lists_in_rows[cell,]
        # remove all neighbors that are not present in the current sample
        true_neighbors <- potential_neighbors[which(potential_neighbors %in% indices_of_available_cells)]
        # assign counts to all available neighbors
        counts_of_true_neighbors <- matrix_samples[index, true_neighbors]
        # count occurences of WT and KO until 50 neighbors were counted
        num_of_neighbors <- size_of_local_neighborhoods+1
        true_neighbors_weighted <- c()
        for (index2 in (1:length(counts_of_true_neighbors))){
          current_count <- counts_of_true_neighbors[index2]
          true_neighbors_weighted <- append(true_neighbors_weighted, rep(true_neighbors[index2], current_count))
        }
        identities_of_neighbors <- seurat_sample$identity[as.integer(true_neighbors_weighted)]
        num_KO <- as.integer(table(identities_of_neighbors[1:num_of_neighbors])["KO"])
        num_WT <- num_of_neighbors-num_KO
        # compute score and write into matrix_with_scores
        if (seurat_sample$identity[cell]=="WT"){
          condition_ratio <- ((num_WT-1)/(num_of_neighbors-1))
        } else {
          condition_ratio <- (num_WT/(num_of_neighbors-1))
        }
        matrix_with_scores[index, cell] <- condition_ratio
      }
    }
  }
  return(matrix_with_scores)
}


.improved___compute_matrix_with_proportion_of_WT_cells_in_neighborhood___bootstrap_iterations_correspond_to_rows <- function(seurat_sample, matrix_samples, num_bootstrap_iter, nn_lists_in_rows, size_of_local_neighborhoods){
  matrix_with_scores <- matrix(0, ncol=ncol(matrix_samples), nrow=nrow(matrix_samples))
  for (index in (1:num_bootstrap_iter)){
    indices_of_available_cells <- which(matrix_samples[index,]!=0)
    for (inner_index in (1:ncol(matrix_samples))){
      cell <- inner_index
      counts <- matrix_samples[index, cell]
      if (counts!=0){
        potential_neighbors <- nn_lists_in_rows[cell,]
        # remove all neighbors that are not present in the current sample
        true_neighbors <- potential_neighbors[which(potential_neighbors %in% indices_of_available_cells)]
        # assign counts to all available neighbors
        counts_of_true_neighbors <- matrix_samples[index, true_neighbors]
        # For the current cell, a vector with indices of its neighbors in the current iteration is computed; neighbors are ordered by increasing distance
        # Indices of neighbors that were sampled multiple times are included as often as they occurred -> "weighted";
        true_neighbors_weighted <- c()
        for (index2 in (1:length(counts_of_true_neighbors))){
          current_count <- counts_of_true_neighbors[index2]
          true_neighbors_weighted <- append(true_neighbors_weighted, rep(true_neighbors[index2], current_count))
        }
        # weighted indices of neighbors are translated into the conditions the cells belong to
        identities_of_neighbors <- seurat_sample$identity[as.integer(true_neighbors_weighted)]
        # In section below, the proportion of KO-cells among all cells in the local neighborhood is calculated
        num_KO <- as.integer(table(identities_of_neighbors[1:size_of_local_neighborhoods])["KO"])
        condition_ratio <- num_KO/size_of_local_neighborhoods
        matrix_with_scores[index, cell] <- condition_ratio
      }
    }
  }
  return(matrix_with_scores)
}



.create_list_with_dataframes___pseudotime_and_score_columns___each_dataframe_corresponds_to_a_bootstrap_iteration <- function(seurat_sample, matrix_with_scores, matrix_samples, pseudo_time_trajectory_to_use){
  list_sampling_results <- list()
  for (index in (1:nrow(matrix_with_scores))){
    # remove all entries from score- and sample-vector where weight = 0
    weight_vector_cleaned <- matrix_samples[index,][which(matrix_samples[index,]!=0)]
    score_vector_cleaned <- matrix_with_scores[index,][which(matrix_samples[index,]!=0)]
    pt_vector_cleaned <- seurat_sample@meta.data[[pseudo_time_trajectory_to_use]][which(matrix_samples[index,]!=0)]
    weighted_scores <- c()
    weighted_pt <- c()
    for (index2 in (1:length(score_vector_cleaned))){
      weighted_scores <- append(weighted_scores, rep(score_vector_cleaned[index2], weight_vector_cleaned[index2]))
      weighted_pt <- append(weighted_pt, rep(pt_vector_cleaned[index2], weight_vector_cleaned[index2]))
    }
    list_sampling_results <- list.append(list_sampling_results, data.frame(score=weighted_scores, pt=weighted_pt))
  }
  return(list_sampling_results)
}

.interpolate_values <- function(list_sampling_results, seurat_sample, pseudo_time_trajectory_to_use, size_moving_window=0.05){
  list_with_interpolated_values <- list()
  for (index in (1:length(list_sampling_results))){
    df_smoothing <- list_sampling_results[[index]][order(list_sampling_results[[index]][,"pt"], decreasing = FALSE),]
    # remove scores without pt (occurs because not every cell is assigned to trajectory of interest)
    df_smoothing <- df_smoothing[which(is.na(df_smoothing$pt)==FALSE),]
    df_smoothing$pt <- (df_smoothing$pt)/max(seurat_sample@meta.data[[pseudo_time_trajectory_to_use]], na.rm = TRUE)
    df_smoothing[,"moving_average"] <- slider::slide_index_dbl(.x=df_smoothing$score,
                                                              .i=df_smoothing$pt,
                                                              .f=mean, .before = size_moving_window, .after = size_moving_window, .complete = FALSE)
    # below duplicates can be removed because they were only important for computing the moving average to ensure that cells which were sampled multiple times obtain more weight in the statistic
    df_smoothing <- df_smoothing[which(!duplicated(df_smoothing$pt)),]
    # approximate points by function
    mean_function <- approxfun(x=df_smoothing$pt, y=df_smoothing$moving_average, method = "linear",
                               rule = 2, f = 0)
    list_with_interpolated_values[[index]] <- mean_function(seq(0,1,0.01))
  }
  return(list_with_interpolated_values)
}

.compute_median_and_ci_of_interpolated_values <- function(list_with_interpolated_values=list_with_interpolated_values){
  df_interpolated_values <- as.data.frame(list_with_interpolated_values)
  colnames(df_interpolated_values) <- c(1:ncol(df_interpolated_values))
  median_vector <- c()
  confidence_vector_upper <- c()
  confidence_vector_lower <- c()
  for (index in (1:nrow(df_interpolated_values))){
    median_vector <- append(median_vector, median(as.numeric(df_interpolated_values[index,])))
    confidence_vector_lower <- append(confidence_vector_lower, quantile(as.numeric(df_interpolated_values[index,]), c(.975, .025))[2])
    confidence_vector_upper <- append(confidence_vector_upper, quantile(as.numeric(df_interpolated_values[index,]), c(.975, .025))[1])
  }
  
  df_confidence <- data.frame(median=median_vector, lower=confidence_vector_lower, upper=confidence_vector_upper, pseudo_time=seq(0,1,0.01))
  return(df_confidence)
}
