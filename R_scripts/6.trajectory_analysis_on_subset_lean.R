# 6.trajectory_analysis_on_subset_lean.R


#Set output paths---------------------------------------------------------------

out_dir <- "../figures/"
#Load libraries ----------------------------------------------------------------
library(SingleR)
library(BiocParallel)
library(Seurat)
library(SingleCellExperiment)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(monocle3)
library(SeuratWrappers)
#Load objects-------------------------------------------------------------------
#Register 4 cores for Windows (uses socket clusters instead of mclapply)
register(SnowParam(workers = 4)) #To let it use multiple cores

# Load subsetted Seurat object (clusters 2, 24, 28)
# srt_rds <- readRDS("../data/traj_srt_subset_lean.rds")
# 
# #Convert Seurt object to cell_data_set (cds)
# cds <- as.cell_data_set(srt_rds, assay = "RNA")
# cds <- estimate_size_factors(cds)
# 
# #Pre-processing ------------------------------------------
# #Redo UMAP and clustering with just the selected clusters
# 
# #PCA on gene expression matrix
# cds <- preprocess_cds(cds, num_dim = 50) #Default dim for PCA
# 
# #Run UMAP on PCs
# cds <- reduce_dimension(cds)
# 
# #Cluster nearby cells on UMAP -> clusters and partitions calculations
# cds <- cluster_cells(cds)
# 
# #Plot by Seurat clusters -> QC visualization of UMAP
# plot_cells(cds, color_cells_by = "seurat_clusters", show_trajectory_graph = FALSE)
# 
# #Plot partitions -> colour cells based on monocle's partitions
# plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
# 
# #Learn Trajectory graph
# #SimplePPT algorithm to fit a reversed graph embedding -> tree skeleton branching
# #Set use_partition = False so Monocle learns a single connected trajectory tree across
# #all cells -> Preventing a cluster from being completely disconnected
# 
# #Learn the trajectory graph
# cds <- learn_graph(cds, use_partition = FALSE)
# 
# plot_cells(cds) #Plot UMAP with generated black trajectory tree skeleton overlaid
# 
# #Order cells in pseudotime
# #Interactive test -> Can manually click on a node to set the root
# cds <- order_cells(cds)
# head(pseudotime(cds)) #Check that they're not NA
# 
# #Plot pseudotime trajectories
# plot_cells(cds,
#            color_cells_by = "pseudotime",
#            show_trajectory_graph = TRUE,
#            label_groups_by_cluster = FALSE,
#            label_leaves = FALSE,
#            label_branch_points = FALSE)
# 
# # save_monocle_objects(cds, "../data/cds_traj_v3", comment = "This is the 3rd version CDS with hand-chosen roots. stored 2026-07-23", verbose = TRUE)


#Loading objects -----------------------------------------------------
cds <- load_monocle_objects(directory_path = "../data/cds_traj_v3")
# Load pre-computed graph_test results from RDS
cds_pr_test_res <- readRDS("../data/cds_pr_test_res_v2.rds")


# saveRDS(cds, "../data/cds_traj_v3.rds")

#Graph test for differentially expressed genes along trajectory ----
#Spatial autocorrelation test to identify genes whos expression changes significantly
#as cells move along the trajectory graph
# cds_pr_test_res <- graph_test(cds,
#                               neighbor_graph = "principal_graph", #test gene expression changes across the learned trajectory skeleton
#                               cores = 4) #Use 4 CPU cores to speed up computation

# save_monocle_objects(cds_pr_test_res, "../data/cds_pr_test_res")
# saveRDS(cds_pr_test_res, "../data/cds_pr_test_res_v2.rds")

pr_deg_ids <- row.names(subset(cds_pr_test_res, q_value < 0.05)) #Filter test results to get vector of gene names that show statistically significant expression changes
print(length(pr_deg_ids)) #15240 #15316 after re-estimating size factor

#Cluster 2 = Mesenchymal stromal cell (11626 cells)
#Cluster 24 = Skeletal myocyte progenitor cell (2255 cells) --> Pick as root
#Cluster 28 = Immature Skeletal myocyte (1122 cells)
#Pathway is Skeletal myocyte progenitor cell -> Immature skeletal mycyte -> mature skeletal myocyte
#Root (t = 0) -> early stage progenitor cell
#Mesenchymal stromal cells (MSC) -> A distinct stromal lineage, may form disconnected partition or separation branch

#Setting trajectory with SMGC as root ------------------------------
# progenitor_cells <- colnames(cds)[colData(cds)$seurat_clusters == 24]
# cds <- order_cells(cds, root_cells = progenitor_cells)
# 
# plot_cells(
#   cds,
#   color_cells_by = "pseudotime",
#   show_trajectory_graph = TRUE,
#   label_cell_groups = TRUE
# )

#Regression analysis in Monocle: ====
#evaluate whether each gene depends on variables such as times, treatments, etc
# cds <- readRDS("../data/cds_traj_v2.rds")

#Extract genes with significant q-values and spatial autocorrelation
top_deg_res <- subset(cds_pr_test_res, q_value < 0.05 & morans_I > 0.25)
pr_deg_ids <- row.names(top_deg_res)
length(pr_deg_ids) #2052 #2575 after re-estimate size factor

#Filter down to trajectory-variable genes (pr_deg_ids from graph_test)
cds_subset <- cds[pr_deg_ids, ] #Subset CDS by genes

#Map gene_short_name to geneIDS
if (!"gene_short_name" %in% colnames(rowData(cds_subset))) {
  rowData(cds_subset)$gene_short_name <- rownames(cds_subset)
}

rowData(cds)$gene_short_name <- rownames(cds)

#Remove NA or infinite pseudotime
valid_cells <- colnames(cds_subset)[is.finite(pseudotime(cds_subset))]
cds_subset <- cds_subset[, valid_cells]

#Make dark/large font theme for monocle3 plots ----------
custom_monocle_theme <- theme(
  text = element_text(color = "black"),
  plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5, margin = margin(b = 10)),
  axis.title = element_text(size = 15, face = "bold", color = "black"),
  axis.text = element_text(size = 12, color = "black"),
  legend.title = element_text(size = 13, face = "bold", color = "black"),
  legend.text = element_text(size = 11, color = "black"),
  strip.text = element_text(size = 12, face = "bold", color = "black")
)


#Plot UMAP with traj path -----
p_umap_pseudo <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE, #Show black trajectory path
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 1.5
) + 
  custom_monocle_theme

ggsave(paste0(out_dir, "pseudo_umap_v4.png"), plot = p_umap_pseudo, width = 11, height = 8, dpi = 600)

# #Fitting GLMs on high-variance genes after =====
# gene_fits <- fit_models(
#   cds_subset,
#   model_formula_str = "~pseudotime",
#   cores = 4
# )

# saveRDS(gene_fits, "../data/gene_fits_traj_v1.rds")====


gene_fits <- readRDS("../data/gene_fits_traj_v2.rds")


#See which genes have time-dependent expression and extract coefficients
fit_coefs <- coefficient_table(gene_fits)


#Get Top 10 Genes that change effect (dynamic genes) (by Moran's I autocorrelation) -------------
top10_dynamic_ids <- cds_pr_test_res %>%
  filter(q_value < 0.05) %>%
  arrange(desc(morans_I)) %>%
  head(10) %>%
  rownames()

top5_dynamic_ids <- cds_pr_test_res %>%
  filter(q_value < 0.05) %>%
  arrange(desc(morans_I)) %>%
  head(5) %>%
  rownames()
top5_summary <- cds_pr_test_res |> 
  filter(q_value < 0.05) |> 
  arrange(desc(morans_I)) |> 
  head(5) |> 
  select(morans_I, q_value)

print(top5_summary)

p_umap_genes <- plot_cells(
  cds,
  genes = top5_dynamic_ids,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE, #Show black trajectory path
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 1.5
)


p_dynamic <- plot_genes_in_pseudotime(
  cds[top5_dynamic_ids, valid_cells],
  color_cells_by = "cell_type",
  min_expr = 0.1,
  label_by_short_name = TRUE,
  trend_formula = "~splines::ns(pseudotime, df = 3)"
) +
  custom_monocle_theme

ggsave(paste0(out_dir, "genes5_dynamic_pseudotime.png"), plot = p_dynamic, width = 11, height = 8, dpi = 600)
#Increase font size of everything and make legend dots bigger
p_dynamic_large <- p_dynamic + 
  theme_light(base_size = 18) +
  guides(color = guide_legend(override.aes = list(size = 5)))
# 
# #Save dynamic genes plot
# ggsave(
#   filename = file.path(out_dir, "top10_dynamic_genes_pseudotime_v2.png"),
#   plot = p_dynamic,
#   width = 18,
#   height = 6,
#   dpi = 300
# )

#Get Top 10 Statistically Significant Genes (by q-value) ----------
top10_sig_ids <- fit_coefs %>%
  filter(term == "pseudotime", !is.na(estimate)) %>%
  arrange(q_value) %>%
  head(10) %>%
  pull(gene_id)

p_sig <- plot_genes_in_pseudotime(
  cds[top10_sig_ids, valid_cells],
  color_cells_by = "cell_type",
  min_expr = 0.1,
  label_by_short_name = TRUE,
  trend_formula = "~splines::ns(pseudotime, df = 3)"
)

# Save significance plot
ggsave(
  filename = file.path(out_dir, "top10_sig_genes_pseudotime.png"),
  plot = p_sig,
  width = 12,
  height = 8,
  dpi = 300
)



#Filter for pseudotime terms with q_value < 0.05 ----
sig_gene_ids <- fit_coefs |> 
  filter(term == "pseudotime", q_value < 0.05, !is.na(estimate)) |> 
  pull(gene_id)

#Subset cds_subset to top genes first
top_genes<- head(sig_gene_ids, 20)

#Plot specific genes across pseudotime
plot_genes_in_pseudotime(cds_subset[top_genes, ],
                         color_cells_by = "cell_type",
                         min_expr = NULL,
                         label_by_short_name = FALSE,
                         trend_formula = "~splines::ns(pseudotime, df = 3)"
)

#Plot expression on UMAP embeddings --------------------
#Genes that score as highly significant according to graph_test():
p_umap_pseudo <- plot_cells(cds, 
     genes = pr_deg_ids,
     show_trajectory_graph = FALSE,
     label_cell_groups = FALSE,
     label_leaves = FALSE)

#Increase font size of everything and make legend bigger
p_umap_pseudo_large <- p_umap_pseudo + 
  theme_light(base_size = 18) +
  guides(color = guide_legend(override.aes = list(size = 5)))


#Group trajectory-variable genes into modules and plot heatmap --------------
#Can collect the trajectory-variable genes into modules:
gene_module_df <- find_gene_modules(
  cds[pr_deg_ids,], 
  resolution=c(10^seq(-6,-1))
  )

#Plot gene modules in heatmap -------------
plot_cells(
  cds, 
  genes = gene_module_df, 
  show_trajectory_graph = FALSE,
  label_cell_groups = FALSE
)
#Pull out the genes that have a significant time component
cell_type_terms |> filter(q_value < 0.05) |> select(gene_short_name, term, q_value, estimate)


#Can select a path with choose_cells() or by subsetting the cell dataset by cluster, cell type, or other annotation
plot_genes_in_pseudotime(cds_subset,
                         color_cells_by = "cell_type",
                         min_expr = 0)

#Analyzing branches in single-cell trajectories:
cds_subset <- choose_cells(cds) 
#Then call graph_test() on the subset -> Identifies genes with interesting patterns of expression that fall only within the region of the trajectory you selected
#Gives more refined and relevant set of genes

#Pseudotime-ordered-heatmap ====
#Aggregate gene expression by module and cell group or pseudotime bin
cell_group_df <- tibble::tibble(
  cell = colnames(cds_subset),
  cell_group = colData(cds_subset)$cell_type # or seurat_clusters
)

rowData(cds_subset)$gene_short_name <- row.names(rowData(cds_subset))
#Plot module expression heatmap across cell types or binned pseudotime
plot_genes_by_group(
  cds[pr_deg_ids, ],
  gene_module_df,
  group_cells_by = "cell_type",
  ordering_type = "maximal_on_diag"
)







save.image()