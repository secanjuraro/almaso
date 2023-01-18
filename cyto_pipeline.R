###############################
## INSTALL REQUIRED PACKAGES ##
###############################

install.packages("remotes")
install.packages("bluster")
install.packages("igraph")


BiocManager::install("DESeq2")
BiocManager::install("FlowSOM")
BiocManager::install("flowCore")


remotes::install_github("LTLA/bluster")


###############################
### LOAD REQUIRED PACKAGES ####
###############################

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
  library(DESeq2)
  library(FlowSOM)
  library(ggcorrplot)
}

############################
### SET WORKING DIR      ###
############################

setwd("~/INSA/BIM 2022-2023/Projet 5BIM/real data")


############################
######## LOAD DATA #########
############################

{
tumor_sample = read.FCS("VP4_Tumor_CD45+ cells.fcs", column.pattern = "Time", invert.pattern = TRUE, truncate_max_range = FALSE)
control_sample = read.FCS("Vp6_Control_CD45+ cells.fcs", column.pattern = "Time", invert.pattern = TRUE, truncate_max_range = FALSE)
mix_sample = read.FCS("MixC_tumeur_CD45+ cells.fcs", column.pattern = "Time", invert.pattern = TRUE, truncate_max_range = FALSE)
immune_control = read.FCS("cd45pos2_control_CD45+.fcs", column.pattern = "Time", invert.pattern = TRUE, truncate_max_range = FALSE)

}

# Get expression matrix as a data frame
expr_matrix <- as.data.frame(immune_control@exprs)

############################
###### PREPROCESSING #######
############################


myTrans <- transformList(names, arcsinhTransform())
expr_matrix<- transform(expr_matrix, myTrans)


##################################
###### DIMENSION REDUCTION #######
##################################

library(plot3D)
library(ggfortify)
library(pca3d)
library(ade4)
library(factoextra)

##PCA##

# 1. Create PCA
pca <- prcomp(expr_matrix[7:19])

# 2. Evaluation of the number of axis
summary(pca)
fviz_eig(pca)
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4])
head(d)
ggplot(data = d, aes_string(x = "PC1", y = "PC2")) + geom_point(size = 3)

############################
####### CLUSTERING #########
############################

BiocManager::install("CytoTree")
library("CytoTree")
?runDiff


### 2 clustering methods:

### A. KNN + Louvain clustering 
### B. HSNE-based Gaussian Mean Shift clustering (Cytosplore)
### C. flowSOM

######### A.KNN + Leuvain clustering #########

## 1. Run UMAP on expression matrix and cast it into a data frame 
# NB: No dimensionality reduction

umap_m <-as.data.frame(umap(d))
colnames(umap_m) <- c("UMAP1","UMAP2")
?cluster_louvain

## 2. Compute clusters based on a SNN Algorithm + Louvain clustering
graph <- bluster::makeSNNGraph(d)
clust <- igraph::cluster_louvain(graph, resolution = 0.5)

## 3. Add cluster membership column to UMAP
umap_m$cluster_id <- factor(clust$membership)

## 4. Plot UMAP and color it by cluster 
ggplot(umap_m, aes(UMAP1, UMAP2, colour = cluster_id)) +  geom_point() +  labs(x = "UMAP1",y = "UMAP2",subtitle = "UMAP plot")


## 5. Get cell's ID (number) for each cluster
communities(clust) # for each cluster get the cells ID 

## 6. Create dataframe with cluster id for each cell 

temp_df_KNN <- data.frame(cell = numeric(), cluster_id = numeric())

for(i in names(clusters_cells)){
  temp_df <- as.data.frame(clusters_cells[i]) 
  temp_df$cluster_id <- i
  colnames(temp_df) <- c("cell", "cluster_id")
  temp_df_KNN <- rbind(temp_df_KNN,temp_df)
}
cluster_id <- temp_df_KNN %>% arrange(cell) %>% select(cluster_id)

df_KNN <- cbind(expr_matrix,cluster_id)



######### B.HSNE-based Gaussian Mean Shift clustering (Cytosplore) #########


## 1. Provide the directory of the fcs files --> already clustered from Cytosplore
dirFCS = "~/INSA/BIM 2022-2023/Projet 5BIM/real data/Clusters - vp6/HSNE-clusters"

## 2. Defining a function to read multiple fcs files from 'dir' into a single data.frame: 
# NB: The column in the output named 'fileName' tracks the original file where each cell came from.
## Defining a function to read multiple fcs files from a directory 'dir' into a single data.frame:
# Optionally perform remapping of column 'CSPLR_ST' holding cytosplore sample numbers to actual names

read.flowdat <- function(dir,path_CSPLR_ST=""){
  # Read:
  filepaths <- list.files(path=dir,pattern = ".fcs", full.names=TRUE)
  flowset <- read.flowSet(files=filepaths, transformation=FALSE, truncate_max_range = FALSE)
  # Transform to data frame:
  x <- as.data.frame(exprs(as(flowset,'flowFrame')),stringsAsFactors=FALSE)
  # Map column 'Original' to filename (in this case holding clusters of HSNE):
  filenames <- gsub("[.fcs]","",list.files(path=dir,pattern = ".fcs", full.names=FALSE))
  names(filenames) <- sort(unique(x$Original))
  x$fileName <- filenames[as.character(x$Original)]
  # Remove column 'Original':
  x <- x[,-which(colnames(x)=="Original")]
  # Optionally remap Cytosplore sample tags to original filename:
  if(file.exists(path_CSPLR_ST)){
    # Read:
    sampID <- gsub(".fcs","",basename(sapply(strsplit(readLines(path_CSPLR_ST),": "),function(x) x[1])))
    names(sampID) <- sapply(strsplit(readLines(path_CSPLR_ST),": "),function(x) x[2])
    x$sampleID <- sampID[as.character(x$CSPLR_ST)]
  }
  return(x)
}

## 3. Read fcs files 
df_HSNE <- read.flowdat(dir=dirFCS)
# Optional: Set columname 'fileName' to clusters_HSNE:
colnames(df_HSNE)[which(colnames(df_HSNE)=="fileName")] <- "clusters_HSNE"

## 4. Visualization
label_HSNE_dm <- df_HSNE%>%group_by(clusters_HSNE)%>%select(CSPLR_HsneDataX, CSPLR_HsneDataY)%>%summarize_all(mean)
ggplot(df_HSNE, aes(x=CSPLR_HsneDataX, y=CSPLR_HsneDataY, color=clusters_HSNE))+geom_point(size=0.1)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+geom_label_repel(aes(label=clusters_HSNE), data=label_HSNE_dm)+guides(colour=FALSE)


######### C. flowSOM ######

## 1. Run function FlowSOM, nClus is set to 10 by default 
flowsom <- FlowSOM(input = immune_control, 
                   transform = FALSE,
                   scale = FALSE,
                   nClus = 10) #we choose 14, since we also generated 14 clusters by HSNE
          

## 2. Get metaclustering per cell
clusters_flowsom <- as.factor(flowsom$map$mapping[,1])
levels(clusters_flowsom) <- flowsom$metaclustering

## 3. Add flowsom clusters to dataframe
df_FlowSOM <- cbind(d, clusters_flowsom)
head(df_FlowSOM)
write.csv(df_FlowSOM, "flowsom_clust_ctrl.csv", sep = ',')

############################
###### VISUALIZATION #######
############################


# select the columns for the UMAP calculation
# check different n_neighbours (controls how UMAP balances local versus global structure in the data) for your UMAP plot
# check min_dist (controls how tightly UMAP is allowed to pack points together, low values=clumpier embeddings) for your UMAP plot
umap <- umap(d, n_neighbors = 30, min_dist=0.001, verbose=TRUE)
umap<- as.data.frame(umap)
colnames(umap) <- c('umap_1', 'umap_2')
d <- cbind(d,umap)
df_FlowSOM <- cbind(df_FlowSOM,umap)
df_HSNE <- cbind(df_HSNE,umap)
df_KNN <- cbind(df_KNN, umap)

#visualize and label clusters on umap for each method
label_HSNE_umap <- df_HSNE%>%group_by(clusters_HSNE)%>%select(umap_1, umap_2)%>%summarize_all(mean)
label_flowsom_umap <- df_FlowSOM%>%group_by(clusters_flowsom)%>%select(umap_1, umap_2)%>%summarize_all(mean)

ggplot(df_HSNE, aes(umap_1, umap_2, colour = as.factor(clusters_HSNE))) +  geom_point() +  labs(x = "UMAP1",y = "UMAP2",subtitle = "UMAP plot HSNE")
ggplot(df_FlowSOM, aes(umap_1, umap_2, colour = as.factor(clusters_flowsom))) +  geom_point() +  labs(x = "UMAP1",y = "UMAP2",subtitle = "UMAP plot FlowSOM")


##visualize a specific marquer expression level in the UMAP 

marker_value <- df_FlowSOM$`SSC-H`
mid<-mean(marker_value)
ggplot(df_FlowSOM, aes(umap_1, umap_2, colour = as.numeric(marker_value))) +  geom_point() +  labs(x = "UMAP1",y = "UMAP2",subtitle = "UMAP plot HSNE")+scale_color_gradient2(midpoint=mid, low="blue", mid="white", high="red", space ="Lab" )
