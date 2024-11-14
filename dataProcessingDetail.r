#UPLOAD LIBRARY
suppressWarnings(library(Seurat))
suppressWarnings(library(tidyverse))
suppressWarnings(library(ggplot2))
library(dplyr)

merged.pbmc <- readRDS("D:\\Tugas Akhir\\Seurat\\merged_pbmc.rds")

#Subset for T cell, NK cell and NKT cell aja
tcell.subset <- readRDS("D:\\Tugas Akhir\\Seurat\\t_nk_nkt_pbmc.rds")
tcell.subset <- subset(merged.pbmc, idents = c("T cell", "NK cell", "NKT cell"))
View(tcell.subset@meta.data)

#re-cluster for T cell subset
DefaultAssay(tcell.subset) <- "RNA"
tcell.subset <- NormalizeData(tcell.subset, verbose = FALSE)
tcell.subset <- FindVariableFeatures(tcell.subset, selection.method="vst", nfeatures = 2000, verbose = FALSE)
tcell.subset <- ScaleData(tcell.subset, verbose = FALSE)
tcell.subset <- RunPCA(tcell.subset, npcs = 30, verbose = FALSE)
tcell.subset <- RunUMAP(tcell.subset, reduction = "pca", dims = 1:30)
tcell.subset <- FindNeighbors(tcell.subset, dims = 1:30)
tcell.subset <- FindClusters(tcell.subset, resolution = 0.5)

subsubpopulation_T <- DimPlot(tcell.subset, reduction = "umap", label = TRUE)
ggsave("umap_subsubpopulation_T_NK_NKT.png", plot = subsubpopulation_T)


tcell.subset <- JoinLayers(tcell.subset)
markers_T <- FindAllMarkers(tcell.subset, test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25)
head(markers_T)
top_markers_T_clusters <- list()

# Iterasi untuk setiap cluster yang ada
for (c in unique(markers_T$cluster)) {
  # Filter marker untuk cluster tertentu
  markers_cluster <- markers_T %>%
    filter(cluster == c) %>%
    arrange(desc(avg_log2FC))  # Mengurutkan berdasarkan nilai p_value yang disesuaikan
  
  # Tampilkan 10 marker teratas untuk cluster tersebut
  top_markers_cluster <- markers_cluster %>%
    slice_head(n = 20)
  
  # Simpan hasil dalam list
  top_markers_T_clusters[[as.character(c)]] <- top_markers_cluster
}

top_markers_T_clusters
write.csv(top_markers_T_clusters, file = "top20_markers_per_cluster_T_NK_NKT.csv", row.names = FALSE)