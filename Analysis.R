library(Seurat)
library(tidyverse)
library(patchwork)

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# 1. Download the 10x PBMC 3k dataset (about 7 MB)
url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
if (!dir.exists("data/filtered_gene_bc_matrices")) {
  download.file(url, "data/pbmc3k.tar.gz", mode = "wb")
  untar("data/pbmc3k.tar.gz", exdir = "data")
}

# 2. Build the Seurat object
counts <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(counts = counts, project = "pbmc3k",
                           min.cells = 3, min.features = 200)

# 3. Quality-control metrics
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

cat("Cells:", ncol(pbmc), " Genes:", nrow(pbmc), "\n")
print(summary(pbmc$nFeature_RNA))
print(summary(pbmc$percent.mt))

qc_plot <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3, pt.size = 0.1)
print(qc_plot)
ggsave("results/qc_violin.png", qc_plot, width = 10, height = 4, dpi = 300)

set.seed(42)

# 4. Filter low-quality cells and likely doublets
n_before <- ncol(pbmc)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                 percent.mt < 5)
cat("Cells before:", n_before, " after:", ncol(pbmc), "\n")

# 5. Normalise, find variable genes, scale, PCA
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), verbose = FALSE)

elbow <- ElbowPlot(pbmc, ndims = 30)
print(elbow)
ggsave("results/pca_elbow.png", elbow, width = 6, height = 4, dpi = 300)

# 6. Cluster and run UMAP (first 10 PCs, a standard choice for this dataset)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

umap <- DimPlot(pbmc, reduction = "umap", label = TRUE) +
  ggtitle("PBMC 3k: unsupervised clusters")
print(umap)
ggsave("results/umap_clusters.png", umap, width = 7, height = 5, dpi = 300)

cat("Cells per cluster:\n")
print(table(Idents(pbmc)))

# 7. Marker genes per cluster
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25)
write.csv(markers, "results/cluster_markers.csv", row.names = FALSE)

top3 <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(avg_log2FC, n = 3) %>%
  dplyr::ungroup()
print(as.data.frame(top3[, c("cluster", "gene", "avg_log2FC", "p_val_adj")]))

# 8. Canonical immune markers
canonical <- c("IL7R", "CCR7", "S100A4", "CD8A", "NKG7", "GNLY", "MS4A1",
               "CD14", "LYZ", "FCGR3A", "MS4A7", "FCER1A", "CST3", "PPBP")
dots <- DotPlot(pbmc, features = canonical) + RotatedAxis()
print(dots)
ggsave("results/canonical_markers_dotplot.png", dots, width = 9, height = 5, dpi = 300)

# 9. Assign cell-type labels (based on marker genes, see README)
new_ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T",
             "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new_ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new_ids)
pbmc$cell_type <- Idents(pbmc)

umap_labeled <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) +
  NoLegend() + ggtitle("PBMC 3k: manually annotated cell types")
print(umap_labeled)
ggsave("results/umap_annotated.png", umap_labeled, width = 7, height = 5, dpi = 300)

print(table(pbmc$cell_type))
saveRDS(pbmc, "data/pbmc_annotated.rds")

BiocManager::install(c("SingleR", "celldex"))

library(SingleR)
library(celldex)

BiocManager::install("alabaster.base", type = "binary")
BiocManager::install("celldex", type = "binary")
library(celldex)

# 10. Marker-signature check of the manual labels
sets <- list(
  Naive_CD4_T  = c("CCR7", "LEF1", "SELL", "TCF7", "IL7R"),
  Memory_CD4_T = c("S100A4", "IL7R", "CD40LG", "AQP3", "CD2"),
  CD8_T        = c("CD8A", "CD8B", "GZMK", "GZMH", "CCL5"),
  NK           = c("GNLY", "NKG7", "GZMB", "PRF1", "KLRD1"),
  B            = c("MS4A1", "CD79A", "CD79B", "CD74", "BANK1"),
  CD14_Mono    = c("CD14", "LYZ", "S100A8", "S100A9", "CST3"),
  FCGR3A_Mono  = c("FCGR3A", "MS4A7", "LST1", "AIF1", "IFITM3"),
  DC           = c("FCER1A", "CLEC10A", "CST3", "HLA-DPB1", "HLA-DQA1"),
  Platelet     = c("PPBP", "PF4", "GNG11", "SPARC", "CLU")
)
sets <- lapply(sets, function(g) intersect(g, rownames(pbmc)))
print(sapply(sets, length))  # genes found per signature

pbmc <- AddModuleScore(pbmc, features = sets, name = "sig_")
meta <- pbmc@meta.data
sig_cols <- paste0("sig_", seq_along(sets))
colnames(meta)[match(sig_cols, colnames(meta))] <- names(sets)

heat_df <- meta %>%
  group_by(cell_type) %>%
  summarise(across(all_of(names(sets)), mean)) %>%
  pivot_longer(-cell_type, names_to = "signature", values_to = "mean_score") %>%
  group_by(signature) %>%
  mutate(z = as.numeric(scale(mean_score))) %>%
  ungroup()

heat <- ggplot(heat_df, aes(signature, cell_type, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       name = "Z-scored\nmean score") +
  labs(title = "Marker signatures by annotated cell type",
       x = "Marker signature", y = "Annotated cell type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(heat)
ggsave("results/signature_check.png", heat, width = 8, height = 5.5, dpi = 300)

# Which signature scores highest in each cell?
meta$top_signature <- names(sets)[max.col(as.matrix(meta[, names(sets)]))]
print(table(Annotated = meta$cell_type, TopSignature = meta$top_signature))

# 11. Save outputs
write.csv(as.data.frame(table(Annotated = meta$cell_type,
                              TopSignature = meta$top_signature)),
          "results/signature_check_table.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")

# 12. Re-save plots with a white background
white <- theme(panel.background = element_rect(fill = "white", colour = NA),
               plot.background  = element_rect(fill = "white", colour = NA))

ggsave("results/qc_violin.png", qc_plot & white, width = 10, height = 4, dpi = 300, bg = "white")
ggsave("results/pca_elbow.png", elbow + white, width = 6, height = 4, dpi = 300, bg = "white")
ggsave("results/umap_clusters.png", umap + white, width = 7, height = 5, dpi = 300, bg = "white")
ggsave("results/canonical_markers_dotplot.png", dots + white, width = 9, height = 5, dpi = 300, bg = "white")
ggsave("results/umap_annotated.png", umap_labeled + white, width = 7, height = 5, dpi = 300, bg = "white")
ggsave("results/signature_check.png", heat + white, width = 8, height = 5.5, dpi = 300, bg = "white")

# 13. Regenerate two figures for readability
Idents(pbmc) <- "cell_type"

dots <- DotPlot(pbmc, features = canonical) + RotatedAxis() + white
print(dots)
ggsave("results/canonical_markers_dotplot.png", dots, width = 9, height = 5, dpi = 300, bg = "white")

umap_labeled <- DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5) +
  NoLegend() + ggtitle("PBMC 3k: manually annotated cell types") + white
print(umap_labeled)
ggsave("results/umap_annotated.png", umap_labeled, width = 7, height = 5, dpi = 300, bg = "white")