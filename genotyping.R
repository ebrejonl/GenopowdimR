
library(purrr)
library(tidyverse)
library(fastcluster)
library(parallel)
library(parallelDist)
library(Rcpp)

sourceCpp("CPP2.cpp")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN FUNCTION : genotyping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
compute_genotype_labels <- function(df, loci_column_start, loci_column_end, NA_weight, n_thresholds, min_common_loci) {
  dist_matrix <- calculate_distances(as.matrix(df[, loci_column_start:loci_column_end]), NA_weight, min_common_loci)
  hc <- hclust(as.dist(dist_matrix))
  threshold_sequence <- seq(0, n_thresholds, by = 1)
  genotype_labels_df <- data.frame(matrix(ncol = length(threshold_sequence), nrow = nrow(df)))
  for (i in seq_along(threshold_sequence)) {
    genotype_labels <- cutree(hc, h = threshold_sequence[i] + 0.01)
    genotype_labels_df[, i] <- as.factor(genotype_labels)
  }
  colnames(genotype_labels_df) <- paste0("Genotype_",threshold_sequence)
  df <- cbind(df, genotype_labels_df)
  return(df)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~ Discovery curve ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
select_random_microsatellites <- function(df, num_loci) {
  loci_names <- unique(gsub("a|b", "", colnames(df)))
  selected_loci <- sample(loci_names, num_loci)
  selected_columns <- unlist(lapply(selected_loci, function(locus) {
    grep(locus, colnames(df), value = TRUE)
  }))
  return(df[, selected_columns])
}

bootstrap_unique_individuals <- function(df, num_loci, NA_weight, n_thresholds, n_bootstrap, min_common_loci) {
  threshold_sequence <- seq(0, n_thresholds, by = 1)
  n_thresholds_length <- length(threshold_sequence)
  unique_individuals <- matrix(nrow = n_bootstrap, ncol = n_thresholds_length)
  for (i in 1:n_bootstrap) {
    print(paste0("Bootstrap ", i))
    sampled_df <- select_random_microsatellites(df, num_loci)
    loci_column_start <- 1
    loci_column_end <- ncol(sampled_df)
    labeled_df <- compute_genotype_labels(sampled_df, loci_column_start, loci_column_end, NA_weight, n_thresholds, min_common_loci)
    for (threshold_index in seq_along(threshold_sequence)) {
      threshold_value <- threshold_sequence[threshold_index]
      col_index <- ncol(labeled_df) - n_thresholds_length + threshold_index
      unique_individuals[i, threshold_index] <- length(unique(labeled_df[, col_index]))
    }
  }
  ci_values <- apply(unique_individuals, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
  return(ci_values)
}




genotype_discovery_curve <- function(df, max_loci, NA_weight, n_thresholds, n_bootstrap, min_common_loci) {
  threshold_sequence <- seq(0, n_thresholds, by = 1)
  n_thresholds_length <- length(threshold_sequence)
  results <- data.frame(
    num_loci = rep(1:max_loci, each = n_thresholds_length),
    threshold = rep(threshold_sequence, times = max_loci),
    lower_CI = numeric(max_loci * n_thresholds_length),
    median = numeric(max_loci * n_thresholds_length),
    upper_CI = numeric(max_loci * n_thresholds_length)
  )
  for (num_loci in 1:max_loci) {
    print(paste0("locus ", num_loci))
    ci_values <- bootstrap_unique_individuals(df, num_loci, NA_weight, n_thresholds, n_bootstrap, min_common_loci)
    for (threshold_index in seq_along(threshold_sequence)) {
      threshold <- threshold_sequence[threshold_index]
      index <- (num_loci - 1) * n_thresholds_length + threshold_index
      results[index, 3:5] <- ci_values[, threshold_index]
    }
  }
  return(results)
}



#####~~~~~~~~~~~~~~~~~~~~~~~~~~~ Output the genotypes only seperated by spaces (for manual double checking ~~~~~~~~~~~~~~~~~~~~~~#####

# Function to add empty rows between groups of different genotypes
make_genet_file <- function(mydata, group_col) {
  # Step 1: Sort the dataframe by the 'group_col'
  mydata_sorted <- mydata %>% arrange(!!sym(group_col))
  # Initialize an empty dataframe to store the result
  result <- data.frame()
  # Get unique genotypes
  unique_genotypes <- unique(mydata_sorted[[group_col]])
  for (genotype in unique_genotypes) {
    # Subset the data frame by the current genotype
    subset_df <- mydata_sorted[mydata_sorted[[group_col]] == genotype, ]
    # Bind the subset to the result data frame
    result <- bind_rows(result, subset_df)
    # Add empty rows with NA values (preserve structure)
    empty_rows <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(mydata_sorted)))
    colnames(empty_rows) <- colnames(mydata_sorted)
    result <- bind_rows(result, empty_rows)
  }
  # Replace NA values with empty strings
  result[is.na(result)] <- ""
  return(result)
}
## insert empty rows in sorted genotype file
insert_empty_rows <- function(df) {
  empty_rows <- data.frame(genotype = "", value = "", stringsAsFactors = FALSE)
  rbind(df, empty_rows, empty_rows)
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Function to calculate probability of identity (Paetkau and Strobeck 1994) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## integration of PopGenUtils prob_identity package
library(adegenet)
library(devtools)
library("PopGenUtils")

# here, we use: 
#PID​(locus)=∑i=1n​pi4​+∑i=1n​∑j=i+1n​(2pi2​pj2​)
#
#The overall PIDPID​ for multiple loci is:
#
#PID(overall)=∏l=1LPID(locusl)PID​(overall)=∏l=1L​PID​(locusl​)

prob_identity <- function(mydata, loci_column_start,loci_column_end){ ## Loci only

microsat_data <-mydata[, loci_column_start:loci_column_end]
#rename column for adegenet
modify_colnames <- function(colnames) {
  colnames <- gsub("\\.", "", colnames) # Remove periods
  colnames <- gsub("\\_", "", colnames) # Remove all underscores
  colnames <- sub("(.*)(.)$", "\\1.\\2", colnames) 
  # Add period before the last character
  return(colnames)
}
colnames(microsat_data) <- modify_colnames(colnames(microsat_data))
  
# merge loci 
unique_loci <- unique(sub("\\.a|\\.b", "", colnames(microsat_data)))
merged_data <- lapply(unique_loci, function(locus) {
  microsat_data %>%
    select(matches(paste0(locus, "\\."))) %>%
    rowwise() %>%
    mutate(merged_value = paste(c_across(), collapse = ",")) %>%
    select(merged_value)
})
# Combine the merged data into a single dataframe
merged_data <- as.data.frame(do.call(cbind, merged_data))
colnames(merged_data) <- unique_loci
# get the data in genind format
genind_obj <- df2genind(merged_data, ploidy = 2,  sep=",", NA.char = "x")
genind_obj@all.names
pid_perm <- pid_permute(obj = genind_obj, 1000)
  return(pid_perm)
}