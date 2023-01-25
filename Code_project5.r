library(Seurat)
library(ggplot2)
library(Matrix)

# Lecture de la matrice de comptage brute
setwd(dir = "D:/Users/VP/Documents/0 - Enfants/00 - Noémie/Projet Génomique médicale")
data_single_cell = ReadMtx('matrix.mtx', cells = 'barcodes.tsv', features = 'features.tsv')

# Création de l'objet Seurat
seu<- CreateSeuratObject(counts = data_single_cell)   #Est-ce qu'il faut mettre des min.cells, min.features ?

# Visualisation du dataset
table(seu@active.ident) 
head(seu@meta.data)
class(seu@meta.data)
colnames(seu@meta.data)

Seurat::FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Seurat::VlnPlot(seu, features = c("nCount_RNA","nFeature_RNA"))

# Calcul du pourcentage de gènes mitochondiraux et ribosomaux
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
FeatureScatter(seu, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(seu, feature1 = "nFeature_RNA", feature2 = "percent.ribo")

# Decide on thresholds to filter out cells :
seu_filtered_1 <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 9500 & nCount_RNA < 150000)
seu_filtered_1 <- subset(seu_filtered_1, subset = percent.mt < 25)

Seurat::FeatureScatter(seu_filtered_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


#### Normalization and scaling:
?Seurat::NormalizeData # log(1+ (count/colSums * scale.factor)) # natural logarithm
seu_norm_log <- Seurat::NormalizeData(seu_filtered_1,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)


# Find variable genes (for PCA):
seu_norm_log_vst <- Seurat::FindVariableFeatures(seu_norm_log,
                                    selection.method = "vst",
                                    nfeatures = 2000)

# what are the top 10 most variable genes?
top10_log_vst <- head(Seurat::VariableFeatures(seu_norm_log_vst), 10)
top10_log_vst

# plot as scatterplot of average expression level (x-axis) versus variance (y-axis)
vf_plot <- Seurat::VariableFeaturePlot(seu_norm_log_vst) # store unlabeled plot in object
vf_plot
Seurat::LabelPoints(plot = vf_plot, points = top10_log_vst, repel = TRUE) # add gene labels of top variable genes


# scaling (for PCA), only the variable genes are scaled by default
seu_norm_log_scaled <- Seurat::ScaleData(seu_norm_log_vst,
                         features = rownames(seu_norm_log))

# Run PCA and plot:
seu_log <- Seurat::RunPCA(seu_norm_log_scaled)
Seurat::DimPlot(seu_log, reduction = "pca")


# Heatmap of gene expression of top genes contributing to each of the 12 first PCs:
Seurat::DimHeatmap(seu_log, dims = 1:12, cells = 500, balanced = TRUE)
# Elbowplot:
Seurat::ElbowPlot(seu_log, ndims = 20)


# Clustering at several resolutions:
seu_int <- Seurat::FindNeighbors(seu_log, dims = 1:10)
seu_int <- Seurat::FindClusters(seu_int, resolution = seq(0.1, 0.8, by=0.1))

#UMAP for visualization
seu_int <- Seurat::RunUMAP(seu_int, dims = 1:10)

# Subdivision of clusters:
library(clustree)
clustree::clustree(seu_int@meta.data[,grep("RNA_snn_res", colnames(seu_int@meta.data))],
                   prefix = "RNA_snn_res.")

Idents(seu_int)
# Visualisation du clustering (la résolution peut être choisie parmi celles qui ont été claculées)
Seurat::DimPlot(seu_int, group.by = "RNA_snn_res.0.2")

#Visualisation des cellules où le gène UPK3B est exprimé
FeaturePlot(seu_int, features = c("UPK3B"))

# Recherche de gènes marqueurs pour chaque cluster
seu_int_1 <- Seurat::SetIdent(seu_int, value = seu_int$RNA_snn_res.0.1) # On choisit la résolution 0.1
markers_all_norm = FindAllMarkers(seu_int_1)

#On regarde les 5 premiers gènes marqueurs surexprimés de chaque cluster, pour ensuite essayer de leur attribuer un type cellulaire
markersnorm = markers_all_norm_1 %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)







#####AUTOMATIC ANNOTATION using reference gene expression profiles (with SingleR)
library(celldex)
library(SingleR)
library(dittoSeq) 

# go back to RNA assay as default:
DefaultAssay(seu_int) <- "RNA"

# Download reference data from the Human Primary Cell Atlas
ref2 <- celldex::HumanPrimaryCellAtlasData()
class(ref2)
?SummarizedExperiment
table(ref2$label.main)

?SingleR

#Annotation automatique => assigne à chaque cellule un type cellulaire suivant la similarité des profils d'expression avec les données de référence ref2
seu_int_SingleR2 <- SingleR::SingleR(test = Seurat::GetAssayData(seu_int, slot = "data"),
                                    ref = ref2,
                                    labels = ref2$label.main)
# the output contains annotation scores and assigned labels per cell:
head(seu_int_SingleR2)
dim(colData(ref2))
dim(seu_int_SingleR2)

# assess cell type annotation
SingleR::plotScoreHeatmap(seu_int_SingleR2)
SingleR::plotDeltaDistribution(seu_int_SingleR2)

# remove low frequency cells from annotation
singleR_labels <- seu_int_SingleR2$labels
t <- table(singleR_labels) # count the total number of cells per cell type
t
#On n'affiche pas sur les graphes suivants les types cellulaires qui présentent moins de 100 cellules
other <- names(t)[t < 100]
other
# assign to NA:
singleR_labels[singleR_labels %in% other] <- "Other" # or NA as in course website


# add the singleR annotation to the seurat's object meta.data:
seu_int$SingleR_annot <- singleR_labels

# UMAP with Seurat or dittoSeq:
Seurat::DimPlot(seu_int, group.by = "SingleR_annot")
dittoSeq::dittoDimPlot(seu_int, "SingleR_annot", size = 0.7) # can take Seurat, SingleCellExperiment, or SummarizedExperiment object.

# barplot of cell type numbers per sample:
dittoSeq::dittoBarPlot(seu_int, var = "SingleR_annot", group.by = "orig.ident")




##### RECLUSTERING DU CLUSTER 0 
#On regarde les sous-clusters dans le cluster 0
subset_clus0 =  subset(seu_int_1, idents = 0)
data_neighbours_sub = FindNeighbors(subset_clus0, dims = 1:10)
data_clusters_sub = FindClusters(data_neighbours_sub, resolution = 0.1)
data_umap_sub = RunUMAP(data_clusters_sub, dim = 1:10) 

DimPlot(data_umap_sub)

#On regarde où le gène UPK3B (marqueur des cellules tumorales) est exprimé, pour distinguer les clusters de cellules saines des cellules tumorales
FeaturePlot(data_umap_sub, features = c("UPK3B"))



#On peut refaire une annotation automatique sur ce sous-clustering
seu_int_SingleR2_sub <- SingleR::SingleR(test = Seurat::GetAssayData(data_umap_sub, slot = "data"),
                                         ref = ref2,
                                         labels = ref2$label.main)
singleR_labels <- seu_int_SingleR2_sub$labels
t <- table(singleR_labels) # count the total number of cells per cell type
t
# add the singleR annotation to the seurat's object meta.data:
data_umap_sub$SingleR_annot <- singleR_labels
Seurat::DimPlot(data_umap_sub, group.by = "SingleR_annot")
dittoSeq::dittoDimPlot(data_umap_sub, "SingleR_annot", size = 0.7) 


#On peut recalculer les gènes différentiellement exprimés dans chacun des sous-clusters
markers_all_sub_06 = FindAllMarkers(data_umap_sub)
markersnorm_sub_06 = markers_all_sub_06 %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)



cluster0.markers <- FindMarkers(data_umap_sub, ident.1 = 0, ident.2 = c(2, 3), min.pct = 0.25)
cluster2.markers <- FindMarkers(data_umap_sub, ident.1 = 2, ident.2 = c(0, 3), min.pct = 0.25)
cluster3.markers <- FindMarkers(data_umap_sub, ident.1 = 3, ident.2 = c(2, 0), min.pct = 0.25)
