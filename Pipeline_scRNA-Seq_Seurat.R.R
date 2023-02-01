library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

##Création d'un objet Seurat à partir de la matrice d'expression
expression_matrix <- ReadMtx(
  mtx = "matrix.mtx.gz", features = "features.tsv.gz",
  cells = "barcodes.tsv.gz"
)
head(expression_matrix)
seurat_object <- CreateSeuratObject(counts = expression_matrix)

##Quality control

seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#QC Graphes

metadata <- seurat_object@meta.data

metadata %>%
  ggplot(aes(x=nCount_RNA)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

metadata %>%
  ggplot(aes(x=nFeature_RNA)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 300)

metadata %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
  geom_point() +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 300)


metadata %>%
  ggplot(aes(x=percent.mt)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 12)


#Sélection des cellules vivantes et avec un compte suffisant
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Log normalisation
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object

#Highly variable features
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat_object), 10)
plot1 <- VariableFeaturePlot(seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Scaling
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)

#Get expression matrix
exp_mat <- GetAssayData(object = seurat_object, slot = "counts")

#Réduction des dimensions par PCA
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
print(seurat_object[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(seurat_object, reduction = "pca")
DimHeatmap(seurat_object, dims = 1, cells = 500, balanced = TRUE)

#Choisir axes PCA
ElbowPlot(seurat_object)

## Clustering par UMAP
seurat_object <- FindNeighbors(seurat_object, dims = 1:5)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
head(Idents(seurat_object), 5)

seurat_object <- RunUMAP(seurat_object, dims = 1:5)
DimPlot(seurat_object, reduction = "umap")

#Find significant markers
markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_mark <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


#FeaturePlot(seurat_object, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))
markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(seurat_object, features = top10$gene) + NoLegend()

#assigner manuellement les noms des clusters en fonction des marqueurs identifiés
#remplacer dans new.cluster.ids les noms des types cellulaires identifiés
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")

names(new.cluster.ids) <- levels(seurat_object)
levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, new.cluster.ids)
DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

