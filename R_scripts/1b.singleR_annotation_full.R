# 1b.singleR_annotation_full



.libPaths(c("/home/ssaab/links/projects/def-itobias/R_libs/4.4.0", .libPaths()))
.libPaths()

# Load Libraries ====
library(SingleR)
library(BiocParallel)
library(Seurat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(GenomicRanges)

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
cat(paste("Using", n_cores, "cores...\n"))

register(MulticoreParam(workers = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")), progressbar = TRUE))

#Set up objects ================================================================

#If objects end in "July 1" -> Removed promiscuous labels
#New atlas with gene symbols as rownames
atlas <- readRDS("/home/ssaab/links/scratch/data/mouse_atlas_sub_top50_markers_geneSymbols_July1.rds")

#Get the reference
embryo_atlas_20pct <- readRDS("/home/ssaab/links/scratch/data/embryo_atlas_25pct_singleR_reference_cells_July1.rds") #Adds 7GB

#Load Seurat object -> All of data -> THIS IS PROJECT DATA
srt_rds <- readRDS("/home/ssaab/links/scratch/data/merged_tumour_only.rds")

#Marker character list 
#Named list in the format SingleR wants it --> named vector list
marker_character_list <- readRDS("/home/ssaab/links/scratch/data/marker_character_list_July1.rds")
#Get log-normalized matrix =====================================================

#Get the log-normalized matrix from the SCE atlas
rownames(atlas) <- rowData(atlas)$gene_symbol #To set rownames to gene symbols
ref_matrix2 <- logcounts(atlas)

#Make SCE conversion (required for singleR)
test_sce <- as.SingleCellExperiment(srt_rds)

#Trying to get just the log-normalized counts from the SCE object
query_counts <- logcounts(test_sce)

#Get the cell type labels
head(colData(atlas)) #Column is label.main <factor>

#Extract labels column and convert to character vector for SingleR
ref_labels <- as.character(colData(atlas)$label.main)

#Check that it looks like a list of cell types
head(ref_labels) #Gives Erythroid, Erythroid, ExE ectoderm, Mesenchyme, etc

#Check that list names match reference labels
all(names(marker_character_list) %in% unique(ref_labels)) #TRUE -> All clear

#REMOVE UNNECESSARY OBJECTS ====================================================
#Keep specific objects for SingleR
keep_list <- c("query_counts", "ref_matrix2", "ref_labels", "marker_character_list", "atlas", "srt_rds")

#Remove everything else
rm(list = setdiff(ls(), keep_list))

#Run garbage collector
gc()

#Make all genes character strings
marker_character_list <- lapply(marker_character_list, as.character)
# RUNNING SINGLER cell label level =============================================

annotations <- SingleR(
  test = query_counts, #Target experiment counts matrix (Actual project data)
  ref = ref_matrix2, #Reference matrix for log expression counts
  labels = ref_labels, #list of labels
  genes = marker_character_list#Named list of genes
  # BPPARAM = bpparam #To make it parallel processing if on Compute Canada
)

head(annotations)
unique(annotations) 

saveRDS(annotations, file = "singleR_fulldata_results2.rds")
cat("SingleR job done, removed promiscuous labels")
