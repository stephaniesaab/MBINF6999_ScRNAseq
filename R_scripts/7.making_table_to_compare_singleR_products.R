#7.Making_table_to_compare_singleR_products
library(SingleCellExperiment)
library(tidyverse)
library(ggridges)
library(plotly)
library(htmlwidgets)

# ==============================================================================
# 1. CORE FUNCTION TO PROCESS A SINGLE SINGLECELLEXPERIMENT OBJECT
# ==============================================================================
analyze_singler_sce <- function(sce_obj, singleR_rds, sample_name, out_dir) {
  
  # Ensure output directory ends with a slash
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (!endsWith(out_dir, "/")) out_dir <- paste0(out_dir, "/")
  
  message(paste("--- Processing Sample:", sample_name, "---"))
  
  # Extract colData as a standard data frame for tidyverse operations
  cell_metadata <- as.data.frame(colData(sce_obj))
  
  # ----------------------------------------------------------------------------
  # Checking missingness & tracking stats
  # ----------------------------------------------------------------------------
  num_clusters  <- length(unique(cell_metadata$singleR_labels))
  na_labels     <- any(is.na(cell_metadata$singleR_labels))
  na_scores     <- any(is.na(singleR_rds$scores))
  min_score     <- min(singleR_rds$scores, na.rm = TRUE)
  max_score     <- max(singleR_rds$scores, na.rm = TRUE)
  median_score  <- median(singleR_rds$scores, na.rm = TRUE)
  na_delta      <- any(is.na(singleR_rds$delta.next))
  
  # Consolidating summary into a single-row tibble
  summary_stats <- tibble(
    Sample = sample_name,
    Num_Clusters = num_clusters,
    NA_Labels = na_labels,
    NA_Scores = na_scores,
    Min_Score = min_score,
    Max_Score = max_score,
    Median_Score = median_score,
    NA_Delta = na_delta
  )
  
  # ----------------------------------------------------------------------------
  # Plotting Distributions (Histograms)
  # ----------------------------------------------------------------------------
  png(paste0(out_dir, sample_name, "_hist_scores.png"), width = 1800, height = 1500, res = 300)
  hist(singleR_rds$scores, col = "skyblue", main = paste(sample_name, "- Histogram of Scores"), xlab = "Scores")
  dev.off()
  
  png(paste0(out_dir, sample_name, "_hist_delta_next.png"), width = 1800, height = 1500, res = 300)
  hist(singleR_rds$delta.next, col = "skyblue", main = paste(sample_name, "- Histogram of Delta scores"), xlab = "Delta scores")
  dev.off()
  
  # Get distribution of low scores (< 0.18)
  low_scores <- singleR_rds$scores[singleR_rds$scores < 0.18]
  
  png(paste0(out_dir, sample_name, "_hist_scores_below_0.18_lean.png"), width = 1800, height = 1500, res = 300)
  hist(low_scores, col = "skyblue", main = paste(sample_name, "- Scores < 0.18"), xlab = "Scores")
  dev.off()
  
  # ----------------------------------------------------------------------------
  # Saving violin and ridgeplots
  # ----------------------------------------------------------------------------
  singleR_df <- as.data.frame(singleR_rds)
  
  cluster_scores <- singleR_df %>%
    group_by(labels) %>%
    summarize(median_score = median(delta.next, na.rm = TRUE)) %>%
    ungroup()
  
  vio_box_plot_cluster_scores <- ggplot(cluster_scores, aes(x = labels, y = median_score, fill = labels)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.6) +
    geom_boxplot(width = 0.15, fill = "white", alpha = 0.8) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
    theme_minimal() +
    labs(
      title = paste(sample_name, "- SingleR Delta Scores across Labels"),
      x = "SingleR assigned label",
      y = "Median delta score of cluster"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  ggsave(paste0(out_dir, sample_name, "_violin_cluster_scores_lean.png"), plot = vio_box_plot_cluster_scores, width = 10, height = 7, dpi = 300)
  
  clean_cluster_scores <- singleR_df %>%
    group_by(labels) %>%
    filter(n() >= 3) %>%
    ungroup() %>%
    mutate(labels = as.factor(labels))
  
  ridgeplot_cluster_scores <- ggplot(clean_cluster_scores, aes(x = delta.next, y = labels, fill = labels)) +
    geom_density_ridges(stat = "density",
                        aes(height = after_stat(density)),
                        scale = 2, alpha = 0.7,
                        position = position_points_jitter(width = 0.005, height = 0),
                        point_shape = 21, point_size = 0.8, point_alpha = 0.2) +
    geom_vline(xintercept = 0.03, linetype = "dashed", color = "red") +
    coord_cartesian(xlim = c(0, 0.05))+
    theme_ridges() +
    labs(title = paste(sample_name, "- Ridgeplot Delta"), x = "Delta score", y = "SingleR assigned labels") +
    theme(legend.position = "none")
  
  ggsave(paste0(out_dir, sample_name, "_ridgeplot_cluster_scores_lean.png"), plot = ridgeplot_cluster_scores, width = 8, height = 11, dpi = 300)
  
  # ----------------------------------------------------------------------------
  # Contingency Tables & Interactive Matrices
  # ----------------------------------------------------------------------------
  contingency_table <- table(cell_metadata$singleR_labels, cell_metadata$cell_type)
  contigency_prop <- prop.table(contingency_table, margin = 2)
  
  df_contingency_table <- as.data.frame(contigency_prop)
  colnames(df_contingency_table) <- c("SingleR", "Manual", "Proportion")
  
  df_cell_counts <- as.data.frame(contingency_table)
  colnames(df_cell_counts) <- c("SingleR", "Manual", "Cell_count")
  
  df_contingency_table <- df_contingency_table %>%
    filter(Proportion > 0.01) %>%
    left_join(df_cell_counts, by = c("SingleR", "Manual")) %>%
    mutate(hover_text = paste0("Manual Label: ", Manual, "<br>",
                               "SingleR Label: ", SingleR, "<br>",
                               "Percentage: ", round(Proportion, 4), "<br>",
                               "Cell count: ", Cell_count))
  
  contingency_matrix <- ggplot(df_contingency_table, aes(x = Manual, y = SingleR, fill = Proportion, text = hover_text)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(title = paste(sample_name, "- SingleR vs. Manual Cluster Comparison"), x = "Manual Labels", y = "SingleR Labels")
  
  interactive_contingency_table <- ggplotly(contingency_matrix, tooltip = "text")
  
  htmlwidgets::saveWidget(
    interactive_contingency_table, 
    file = paste0(normalizePath(out_dir), "/", sample_name, "_interactive_contingency_table_lean.html"), 
    selfcontained = TRUE
  )
  
  return(list(
    stats = summary_stats,
    contingency = contingency_table
  ))
}

# ==============================================================================
# 2. RUNNING THE WRAPPER PIPELINE ACROSS MULTIPLE SAMPLES
# ==============================================================================

# Replace these placeholders with your actual loaded SingleCellExperiment data objects
# e.g., sce_list <- list("Sample_A" = sce_a, "Sample_B" = sce_b, "Sample_C" = sce_c)
sce_list     <- list("Sample1" = sce_obj1, "Sample2" = sce_obj2, "Sample3" = sce_obj3)
singler_list <- list("Sample1" = singleR_rds1, "Sample2" = singleR_rds2, "Sample3" = singleR_rds3)
output_directory <- "./singler_sce_analysis/"

all_compiled_stats <- list()
all_contingency_tables <- list()

for (sample in names(sce_list)) {
  result <- analyze_singler_sce(
    sce_obj = sce_list[[sample]], 
    singleR_rds = singler_list[[sample]], 
    sample_name = sample, 
    out_dir = output_directory
  )
  
  all_compiled_stats[[sample]] <- result$stats
  all_contingency_tables[[sample]] <- result$contingency
}

# Final table displaying your metrics side-by-side
master_stats_table <- bind_rows(all_compiled_stats)
print("--- Master Statistics Summary Table ---")
print(master_stats_table)

# Save Master Summary Dataframe
write.csv(master_stats_table, paste0(output_directory, "master_singleR_summary.csv"), row.names = FALSE)

# ==============================================================================
# 3. STATISTICAL EVALUATION: CHI-SQUARED TEST EXAMPLE
# ==============================================================================
message("\n--- Executing Chi-Squared Independence Test Example for Sample1 ---")

test_table <- all_contingency_tables[["Sample1"]]
chi_test_results <- chisq.test(test_table)
print(chi_test_results)