#3.graphing_new_labels_hpc.R

#Set output paths---------------------------------------------------------------

.libPaths(c("/home/ssaab/links/projects/def-itobias/R_libs/4.4.0", .libPaths()))
.libPaths()

out_dir <- "~/links/scratch/plots/"
#Load libraries ----------------------------------------------------------------
library(SingleR)
library(BiocParallel)
library(Seurat)
library(SingleCellExperiment)
library(org.Mm.eg.db)
library(SummarizedExperiment)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(ggridges)
library(plotly)
library(cowplot)
library(htmlwidgets) 

#Load objects-------------------------------------------------------------------

#Data Seurat object
srt_rds <- readRDS("~/links/scratch/data/merged_labels_seurat_object_full.rds")

#SingleR object
singleR_rds <- readRDS("~/links/scratch/data/singleR_fulldata_results.rds")

#Look at new labels on UMAP ----------------------------------------------------
dim_plot_singleR_clusters <- Seurat::DimPlot(srt_rds,
                reduction = "umap",
                group.by = "singleR_labels"
)
cat("singleR labeled ", length(unique(srt_rds$singleR_labels)), "number of clusters \n")



dim_plot_manual_clusters <- Seurat::DimPlot(srt_rds,
                reduction = "umap",
                group.by = "cell_type"
)

#Save SingleR plot
save_plot("~/links/scratch/plots/dim_plot_singleR.png", dim_plot_singleR_clusters, base_height = 7, base_width = 8)

# Save Manual plot
save_plot("~/links/scratch/plots/dim_plot_manual.png", dim_plot_manual_clusters, base_height = 7, base_width = 8)

#Checking missingness-----------------------------------------------------------
cat("are any singleR_labels NA ", any(is.na(srt_rds$singleR_labels)), "\n")
cat("are any singleR_scores NA ", any(is.na(singleR_rds$scores)), "\n")
cat("the minimum singleR score is ", min(singleR_rds$scores), "\n")
cat("the max singleR score is ", max(singleR_rds$scores), "\n")
cat("the median singleR score is ", median(singleR_rds$scores), "\n") 
cat("are any delta scores NA ", any(is.na(singleR_rds$delta.next)), "\n")

#Get distribution of the delta scores-------------------------------------------

#Histogram of the scores
png(paste0(out_dir, "hist_scores.png"), width = 1800, height = 1500, res = 300)
hist(singleR_rds$scores, 
     col = "skyblue", 
     main = "Histogram of Delta scores", 
     xlab = "Delta scores")
dev.off()

png(paste0(out_dir, "hist_delta_next.png"), width = 1800, height = 1500, res = 300)
hist(singleR_rds$delta.next,
     col = "skyblue",
     main = "Histogram of Delta scores",
     xlab = "Delta scores")
dev.off()

#Get distribution of low delta scores
print(table(singleR_rds$scores > 0.18, useNA = "always"))

delta_matrix <- singleR_rds$scores
low_scores <- delta_matrix[delta_matrix < 0.18]
cat("number of scores < 0.18 is ", length(low_scores), "\n")

png(paste0(out_dir, "hist_scores_below_0.18.png"), width = 1800, height = 1500, res = 300)
hist(low_scores,
     col = "skyblue",
     main = "Histogram of low Delta scores (below median)",
     xlab = "Delta scores")
dev.off()

#Saving violin and ridgeplots---------------------------------------------------

#Getting delta scores per cluster
#Get median delta score per cluster
cluster_scores <- singleR_rds %>%
  as.data.frame() %>%
  group_by(labels) %>%
  summarize(median_score = median(delta.next, na.rm = FALSE)) %>%
  as.data.frame()

#Distribution of scores and clusters
vio_box_plot_cluster_scores <- ggplot(cluster_scores, aes(x = labels, y = median_score, fill = labels)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme_minimal() +
  labs(
    title = "Distribution of SingleR Delta Scores across Labels",
    x = "SingleR assigned label",
    y = "Median delta score of cluster"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

#Save violin plot of scores and clusters
ggsave(paste0(out_dir, "violin_cluster_scores.png"), plot = vio_box_plot_cluster_scores, width = 10, height = 7, dpi = 300)

#Filter clusters with 0 or 1 cells, get ones with enough cells for variance calculation
clean_cluster_scores <- singleR_rds %>%
  as.data.frame() %>%
  group_by(labels) %>%
  filter(n() >= 3) %>% #Need at least three for variance
  ungroup() %>%
  mutate(labels = as.factor(labels))

cat("max delta.next in clean clusters: ", max(clean_cluster_scores$delta.next), "\n")
cat("median delta.next in clean clusters: ", median(clean_cluster_scores$delta.next), "\n")

ridgeplot_cluster_scores <- ggplot(clean_cluster_scores, aes(x = delta.next, y = labels, fill = labels)) +
  geom_density_ridges(stat = "density",
                      aes(height = after_stat(density)),
                      scale = 2,
                      alpha = 0.7,
                      position = position_points_jitter(width = 0.005, height = 0),
                      point_shape = 21,
                      point_size = 0.8,
                      point_alpha = 0.2) +
  geom_vline(xintercept = 0.03, linetype = "dashed", color = "red") +
  coord_cartesian(xlim = c(0, 0.05))+
  theme_ridges() +
  labs(x = "Delta score", y = "SingleR assigned labels") +
  theme(legend.position = "none")

#Save ridgeline plot of cleaned clusters (> 3 cells)
ggsave(paste0(out_dir, "ridgeplot_cluster_scores.png"), plot = ridgeplot_cluster_scores, width = 8, height = 11, dpi = 300)

#Saving interactive plot and contingency table----------------------------------
#Cross-tabulation table / contingency table -> statistical tool used to display frequency distribution of categorical variables
#Get the metadata into a dataframe
meta_data <- srt_rds@meta.data

#Get cross-tabulation frequency
contingency_table <- table(meta_data$singleR_labels, meta_data$cell_type)

#Convert to percentage to account for different sizes of cell populations
contigency_prop <- prop.table(contingency_table, margin = 2) #Account for differences in cell populations of clusters

#Make plot so label shows when you hover
df_contingency_table <- as.data.frame(contigency_prop)
colnames(df_contingency_table) <- c("SingleR", "Manual", "Proportion")

#Filtering low-proportion clusters
df_contingency_table <- df_contingency_table %>%
  filter(Proportion > 0.01)

df_cell_counts <- as.data.frame(contingency_table)
colnames(df_cell_counts) <- c("SingleR", "Manual", "Cell_count")
df_contingency_table <- left_join(df_contingency_table, df_cell_counts, by = c("SingleR", "Manual"))

print(head(df_cell_counts))

#Create string for hover box
df_contingency_table <- df_contingency_table %>%
  mutate(hover_text = paste0("Manual Label:", Manual, "<br>", #Adds line break
                             "SingleR Label:", SingleR, "<br>",
                             "Percentage:", Proportion, "<br>",
                             "Cell count:", Cell_count))

contingency_matrix <- ggplot(df_contingency_table, aes(x = Manual, y = SingleR, fill = Proportion, text = hover_text)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(title = "SingleR vs. Manual Cluster Comparison", x = "Manual Labels", y = "SingleR Labels")

#Making it interactive plot
interactive_contingency_table <- ggplotly(contingency_matrix, tooltip = "text") #Make it use manually created string when hovering

htmlwidgets::saveWidget(interactive_contingency_table, file=paste0(out_dir, "interactive_contingency_table.html"), selfcontained = TRUE)
