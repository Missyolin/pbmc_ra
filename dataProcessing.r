#UPLOAD LIBRARY
suppressWarnings(library(Seurat))
suppressWarnings(library(tidyverse))
suppressWarnings(library(ggplot2))

#CHAPTER 1 - Upload Data
ctrl_fname <- Read10X_h5("D:\\Tugas Akhir\\Seurat\\ra\\sc5p_v2_hs_PBMC_10k_filtered_feature_bc_matrix.h5")
gene_expression_ctrl <- ctrl_fname[["Gene Expression"]]
ra_fname <- Read10X_h5("D:\\Tugas Akhir\\Seurat\\ra\\GSM4819747_RA_filtered_feature_bc_matrix.h5")

#Create Seurat Object
ra <- CreateSeuratObject(counts = ra_fname, project = "RA", min.cells = 3, min.features =200)
ctrl <- CreateSeuratObject(counts = gene_expression_ctrl, project = "CTRL", min.cells = 3, min.features = 200)

ra[["percent.mt"]] <- PercentageFeatureSet(ra, pattern = "^MT-")
ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^MT-")

#merged Seurat object with our data
pbmc.list <- list(ctrl, ra)
merged.pbmc <- merge(x= ctrl, y=ra, add.cell.ids = c("CTRL", "RA"))

View(merged.pbmc@meta.data)

#CHAPTER 2 - Adding MetaData
library(stringr)
sample <- names(merged.pbmc@active.ident)
sample_detect <- ifelse(str_detect(sample,"CTRL"), "CTRL", "RA")

merged.pbmc@meta.data$sample <- sample_detect

Idents(object = merged.pbmc) <- "sample"

#CHAPTER 3- Normalized Data and Integrate the two dataset
pbmc.list <- SplitObject(merged.pbmc, split.by = "sample")

#standard preprocessing on each object
for (i in 1: length(pbmc.list)) {
    pbmc.list[[i]] <- NormalizeData(pbmc.list[[i]], verbose = FALSE)
    pbmc.list[[i]] <- subset(pbmc.list[[i]], subset = nFeature_RNA > 300 & nFeature_RNA <2000 & percent.mt < 5 )
    pbmc.list[[i]] <- FindVariableFeatures(pbmc.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

# selected features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = pbmc.list)
pbmc.list <- lapply(X = pbmc.list, FUN = function (x) {
    x <- ScaleData(x, features = features, verbose= FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

#CHAPTER 4 - Integrate Data
#find anchors
anchors <- FindIntegrationAnchors(object.list = pbmc.list)

#integrate data
merged.pbmc <- IntegrateData(anchorset = anchors)

#CHAPTER 5 - Dimensional Reduction for Visualization
# run the standard pipeline for visualization and clustering
merged.pbmc <- ScaleData(merged.pbmc, verbose = FALSE)
merged.pbmc <- FindVariableFeatures(merged.pbmc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

merged.pbmc <- RunPCA(merged.pbmc, npcs = 30, verbose = FALSE)
merged.pbmc <- RunUMAP(merged.pbmc, reduction = "pca", dims = 1:30)

merged.pbmc <- FindNeighbors(merged.pbmc, reduction = "pca", dims = 1:30)

merged.pbmc <- FindClusters(merged.pbmc, resolution=0.5)

#visualization
plot1 <- DimPlot(merged.pbmc, reduction = "umap", label = TRUE)
plot2 <- DimPlot(merged.pbmc, reduction = "umap", split.by = "orig.ident")
plot1|plot2
plot2
ggsave("UMAP_cluster_unannotated_SplitbyCondition.png", plot = plot2)

#CHAPTER 6 - Finding Gene Marker for each cluster
library(dplyr)
markers <- FindAllMarkers(merged.pbmc, test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25)

top_markers_all_clusters <- list()

# Iterasi untuk setiap cluster yang ada
for (c in unique(markers$cluster)) {
  # Filter marker untuk cluster tertentu
  markers_cluster <- markers %>%
    filter(cluster == c) %>%
    arrange(desc(avg_log2FC))  # Mengurutkan berdasarkan nilai p_value yang disesuaikan
  
  # Tampilkan 10 marker teratas untuk cluster tersebut
  top_markers_cluster <- markers_cluster %>%
    slice_head(n = 20)
  
  # Simpan hasil dalam list
  top_markers_all_clusters[[as.character(c)]] <- top_markers_cluster
}

top_markers_all_clusters
write.csv(top_markers_all_clusters, file = "top20_markers_per_cluster.csv", row.names = FALSE)

featurePlot5 <- FeaturePlot(merged.pbmc, features = c("CD3D", "CD3E", "NKG7","KLRD1","CD19", "CD79A", "CD14", "CDKN1C", "PPBP", "PF4"), min.cutoff = "q10", pt.size =2)
featurePlot5
ggsave("Feature_Plot_5_major_cell_types.png", plot = featurePlot5)

# Save merged.pbmc as an RDS file
saveRDS(merged.pbmc, file = "merged_pbmc.rds")

#Read RDS
merged.pbmc <- readRDS("D:\\Tugas Akhir\\Seurat\\merged_pbmc.rds")
View(merged.pbmc@meta.data)

#Rename cluster to 6 big subpopulation (include NKT cell and Dendritic cell)
levels(merged.pbmc)
merged.pbmc <- RenameIdents(merged.pbmc, "0" = "T cell", "1" = "T cell", "2" = "Myeloid cell", "3" = "NKT cell", "4" = "T cell", "5" = "B cell", "6" = "NK cell", "7" = "T cell", "8" = "B cell", "9" = "T cell", "10" = "Dendritic cell", "11" = "Platelet")

table(Idents(object = merged.pbmc))

annotated_cluster <- DimPlot(merged.pbmc, label = TRUE)
ggsave("annotated_cluster_bigSubPopulation.png", plot = annotated_cluster)
