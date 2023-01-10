
## Required packages 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("flowCore")

install.packages("remotes")
remotes::install_github("LTLA/bluster")
install.packages("bluster")
install.packages("igraph")


## Load libraries
{


library(flowCore)
library(data.table)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(dplyr)
library(uwot)
  
}


#Load data
{
tumor_sample = read.FCS("data/VP4_Tumor_CD45+ cells.fcs")
control_sample = read.FCS("data/Vp6_Control_CD45+ cells.fcs")
mix_sample = read.FCS("data/MixC_tumeur_CD45+ cells.fcs")
immune_control = read.FCS("data/cd45pos2_control_CD45+.fcs")
}


# Get expression matrix from cyto data
expr_matrix <- as.data.frame(immune_control@exprs)

# #transpose data frame 
# expr_matrix_t <- transpose(expr_matrix)
# 
# 
# #redefine row and column names
# rownames(expr_matrix_t) <- colnames(expr_matrix)
# colnames(expr_matrix_t) <- rownames(expr_matrix)


#Run UMAP on expression matrix and cast it into a data frame

umap_m <-as.data.frame(umap(expr_matrix))
colnames(umap_m) <- c("UMAP1","UMAP2")
colnames(umap_m)
typeof(umap_m)

#Plot UMAP

ggplot(data = umap_m, aes(x = UMAP1, y = UMAP2)) + geom_point()
umap_m %>% ggplot(aes(x = UMAP1, y = UMAP2)+ geom_point())+ labs(x = "UMAP1",y = "UMAP2",subtitle = "UMAP plot")

# Compute clusters based on a SNN Algorithm + Louvain clustering
graph <- bluster::makeSNNGraph(expr_matrix)
clust <- igraph::cluster_louvain(graph)
# Add cluster membership column to UMAP 
umap_m$clust <- factor(clust$membership)
# Plot UMAP and color it by cluster 
ggplot(umap_m, aes(UMAP1, UMAP2, colour = clust)) +  geom_point()

