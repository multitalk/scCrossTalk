library(SeuratWrappers)
library(monocle3)
# obj_neu <- subset(obj_seurat, subset = celltype == "Neu")
obj_seurat <- readRDS("/path/to/obj_neu.rds")
obj_cds <- as.cell_data_set(obj_seurat)
obj_cds <- cluster_cells(cds = obj_cds, reduction_method = "UMAP")
cellcluster <- obj_seurat$subtype
obj_cds@clusters@listData[["UMAP"]][["clusters"]] <- factor(obj_mp$subtype)
obj_cds <- learn_graph(obj_cds, use_partition = FALSE)
obj_cds <- order_cells(obj_cds)
plot_cells(
  cds = obj_cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

# obj_tnk <- subset(obj_seurat, subset = celltype == "T_NK")
obj_seurat <- readRDS("/path/to/obj_tnk.rds")
obj_cds <- as.cell_data_set(obj_seurat)
obj_cds <- cluster_cells(cds = obj_cds, reduction_method = "UMAP")
cellcluster <- obj_seurat$subtype
obj_cds@clusters@listData[["UMAP"]][["clusters"]] <- factor(obj_mp$subtype)
obj_cds <- learn_graph(obj_cds, use_partition = FALSE)
obj_cds <- order_cells(obj_cds)
plot_cells(
  cds = obj_cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)