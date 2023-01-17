#### Analyse des données cytométrie du projet 5BiM



############################
### SET WORKING DIR      ###
############################

setwd("~/INSAetude/5A/Projet_5BIM/Vrais_donnees_cyto")

###############################
## INSTALL REQUIRED PACKAGES ##
###############################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("remotes")
install.packages("igraph")
install.packages("agricolae")
install.packages("ggpubr")

BiocManager::install("DESeq2")
BiocManager::install("FlowSOM")
BiocManager::install("flowCore")
BiocManager::install("cytoFast")

BiocManager::install("bluster")
BiocManager::install("limma")
BiocManager::install("diffcyt")

install.packages("compareGroups")


BiocManager::install("MetaCyto")

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
  library(bluster)
  
  library(flowWorkspace)
  library(ggcyto)
  library(flowAI)
  library(cytutils)
  library(patchwork)
  
  library(cytofast)
  library(MetaCyton)
  library(limma)
  
  library(agricolae)
  library(ggpubr)
  
  library(compareGroups)
  library(xml2)
  library(htmltools)
  library(knitr)
  library(kableExtra)
  library(data.table)
  
  library(diffcyt)
}


############################
######## LOAD DATA #########
############################

{
  tumor_sample = read.FCS("VP4_Tumor_CD45+ cells.fcs", column.pattern = "Time", invert.pattern = TRUE)
  control_sample = read.FCS("Vp6_Control_CD45+ cells.fcs", column.pattern = "Time", invert.pattern = TRUE)
  mix_sample = read.FCS("MixC_tumeur_CD45+ cells.fcs", column.pattern = "Time", invert.pattern = TRUE)
  immune_control = read.FCS("cd45pos2_control_CD45+.fcs", column.pattern = "Time", invert.pattern = TRUE)
  
}

control_sample = read.FCS("cd45pos2_control_CD45+.fcs", column.pattern = "Time", invert.pattern = TRUE, truncate_max_range = FALSE)


# few command to play with data
control_sample
names(control_sample)
colnames(control_sample)
description(control_sample)
summary(control_sample)
keyword(control_sample)
each_col(control_sample, median)
colnames(control_sample)


############################
###### PREPROCESSING #######
############################


#### ça vient de là: https://github.com/hally166/Cytometry-R-scripts/blob/master/compensate_transform_clean.R 
#Compensation
spillover(control_sample)
data_comp <-compensate(control_sample, spillover(control_sample)$SPILL)
data_comp

#Cleaning
data_comp_clean <- flow_auto_qc(control_sample)          #data_comp)
data_comp_clean
keyword(data_comp_clean) <- keyword(control_sample)
data_comp_clean

#Transformation
trans <- estimateLogicle(data_comp_clean, colnames(data_comp_clean[,3:10]))
data_comp_clean_trans <- transform(data_comp_clean, trans)

trans <- estimateLogicle(control_sample, colnames(control_sample[,3:10]))
control_sample <- transform(control_sample, trans)

#Visualise the results
??ggcyto
autoplot(data_comp_clean)
autoplot(data_comp_clean_trans)
autoplot(data_comp_clean_trans, x="PE-Cy7-A", y="PerCP-Cy5-5-A", bins = 256)
autoplot(data_comp_clean_trans, x="Time", y="FSC-A", bins = 128)
autoplot(transform(data_comp,trans), x="Time", y="FSC-A", bins = 128)


# Get expression matrix as a data frame
expr_matrix <- as.data.frame(control_sample@exprs)
expr_matrix <- expr_matrix[,7:19]; expr_matrix

############################
####### CLUSTERING #########
############################

### 2 clustering methods:

### A. KNN + Leuvain clustering 
### B. HSNE-based Gaussian Mean Shift clustering (Cytosplore)
### C. flowSOM



######### A.KNN + Leuvain clustering #########

## 1. Run UMAP on expression matrix and cast it into a data frame 
# NB: No dimensionality reduction

umap_m <-as.data.frame(umap(expr_matrix))
colnames(umap_m) <- c("UMAP1","UMAP2")

## 2. Compute clusters based on a SNN Algorithm + Louvain clustering
graph <- bluster::makeSNNGraph(expr_matrix)
clust <- igraph::cluster_louvain(graph)

## 3. Add cluster membership column to UMAP
umap_m$cluster_id <- factor(clust$membership)

## 4. Plot UMAP and color it by cluster 
ggplot(umap_m, aes(UMAP1, UMAP2, colour = cluster_id)) +  geom_point() +  labs(x = "UMAP1",y = "UMAP2",subtitle = "UMAP plot")


## 5. Get cell's ID (number) for each cluster
communities(clust) # for each cluster get the cells ID 





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

set.seed(47)
colnames(control_sample)

## 1. Create a flowSOM object, nClus is set to 10 by default 
fSOM <- FlowSOM(input = control_sample,
                  compensate = FALSE,
                  transform = FALSE,
                  scale = FALSE,
                  # SOM options:
                  colsToUse = c(7:19),
                  # Metaclustering options:
                  nClus = 10) #we choose 14, since we also generated 14 clusters by HSNE


## 2. Get metaclustering per cell
clusters_flowsom <- as.factor(fSOM$map$mapping[,1]); clusters_flowsom
levels(clusters_flowsom) <- fSOM$metaclustering

## 3. Add flowsom clusters to dataframe
df_FlowSOM <- cbind(expr_matrix, clusters_flowsom)


## 4. get a pdf file with a lot of informations
FlowSOMmary(fsom = fSOM,plotFile = "FlowSOMmary.pdf")

Clust_cv <- GetClusterCVs(fSOM);Clust_cv
clusters <- GetClusters(fSOM); clusters
PlotFlowSOM(fSOM, view= "grid")
Plot2DScatters(fSOM, channelpairs = )

############################
###### VISUALIZATION #######
############################


# select the columns for the UMAP calculation
# check different n_neighbours (controls how UMAP balances local versus global structure in the data) for your UMAP plot
# check min_dist (controls how tightly UMAP is allowed to pack points together, low values=clumpier embeddings) for your UMAP plot
umap <- umap(expr_matrix, n_neighbors = 30, min_dist=0.001, verbose=TRUE)
umap<- as.data.frame(umap)
colnames(umap) <- c('umap_1', 'umap_2')
expr_matrix <- cbind(expr_matrix,umap)
df_FlowSOM <- cbind(df_FlowSOM,umap)
df_HSNE <- cbind(df_HSNE,umap)


#visualize and label clusters on umap for each method
label_HSNE_umap <- df_HSNE%>%group_by(clusters_HSNE)%>%select(umap_1, umap_2)%>%summarize_all(mean)
label_flowsom_umap <- df_FlowSOM%>%group_by(clusters_flowsom)%>%select(umap_1, umap_2)%>%summarize_all(mean)

ggplot(df_HSNE, aes(umap_1, umap_2, colour = as.factor(clusters_HSNE))) +  geom_point() +  labs(x = "UMAP1",y = "UMAP2",subtitle = "UMAP plot HSNE")
ggplot(df_FlowSOM, aes(umap_1, umap_2, colour = as.factor(clusters_flowsom))) +  geom_point() +  labs(x = "UMAP1",y = "UMAP2",subtitle = "UMAP plot FlowSOM")


##visualize a specific marquer expression level in the UMAP 

marker_value <- df_FlowSOM$`SSC-H`
mid<-mean(marker_value)
ggplot(df_FlowSOM, aes(umap_1, umap_2, colour = as.numeric(marker_value))) +  geom_point() +  labs(x = "UMAP1",y = "UMAP2",subtitle = "UMAP plot HSNE")+scale_color_gradient2(midpoint=mid, low="blue", mid="white", high="red", space ="Lab" )




####################### wesh les résultats sont trop moches :
## on va chercher quelle est le marqueur le plus exprimé dans le cluster 3: 
# Get the count matrix
percentages <- GetFeatures(fsom = fSOM, 
                           files = "Vp6_Control_CD45+ cells.fcs", 
                           type = "percentages")

PlotMarker(fSOM, c(7:19))



# Perform the statistics
MC_stats <- GroupStats(percentages[["metacluster_percentages"]], control_sample)
C_stats <- GroupStats(percentages[["cluster_percentages"]], control_sample)

# comme on a : df_FlowSOM = (expr_matrix + clusters_flowsom)
cluster3 <- df_FlowSOM[df_FlowSOM$clusters_flowsom==3,]; head(cluster3)
markerFinder(control_sample)


##### USE OF STAT: PAIREWISE WILCOX:
df_FlowSOM2 <- df_FlowSOM %>% group_by(clusters_flowsom) %>% summarise(across(everything(), mean, na.rm=TRUE))  
df_FlowSOM2
res_wilcox_col1<-pairwise.wilcox.test(df_FlowSOM2$`FJComp-APC-A`, df_FlowSOM2$clusters_flowsom,p.adjust.method = "BH"); res_wilcox_col1
#res_wilcox_col1<- apply(res_wilcox_col1, function(x){res_wilcox_col1[x,]==0 res_wilcox_col1[x,]=0 })

###### USE OF STAT: TUCKEY.HSD() ######
colnames(df_FlowSOM)
df_FlowSOM
d_counts<-group_by(df_FlowSOM, clusters_flowsom) %>%summarise(count = n())

List <- names(df_FlowSOM); List 
List<-List[-(14)]# select just the variables
List
model1 <- lapply(List, function(x) { lm(substitute(i~clusters_flowsom, list(i = as.name(x))), data = df_FlowSOM)})
lapply(model1, summary)
letters = lapply(model1, function(m) HSD.test((m), "clusters_flowsom", group = TRUE, console = FALSE))



### Graph boxplot pour premier marqueur
levels(df_FlowSOM$clusters_flowsom)
ggboxplot(df_FlowSOM, x = "clusters_flowsom", y = "`FJComp-APC-A`", 
          color = "clusters_flowsom", ylab = "FJComp-APC-A", xlab = "Cluster")

### Graph bowplot avec 4 marqueurs
List2<-names(df_FlowSOM)[1:4]; List2
lapply(List2, function(z){ ggboxplot( df_FlowSOM, x="clusters_flowsom", y = substitute(i,list(i=as.name(z))) , color= "clusters_flowsom", xlab="Clusters"  )})





# Make an ANOVA!
data <- data.frame(group = rep(c("P1", "P2", "P3"), each = 40),
                   values = c(rnorm(40, 0, 3),rnorm (40, 0, 6),rnorm (40, 1, 5)))
head(data)
model <- aov(`FJComp-APC-A`~clusters_flowsom, data=df_FlowSOM)
summary(model)
#Make a tuckey test
TukeyHSD(model, conf.level=.95)
plot(TukeyHSD(model, conf.level=.95), las = 2)




###### USE COMPAREGROUPS
df_FlowSOM2 <- df_FlowSOM %>% group_by(clusters_flowsom) %>% summarise(across(everything(), mean, na.rm=TRUE))  
df_FlowSOM2

compareGroups(clusters_flowsom ~ . , data=df_FlowSOM2) #, max????)

#%>% remove_rownames %>% column_to_rownames(var="clusters_flowsom")


##########################################
##########  diffcyt   ####################
##########################################
suppressPackageStartupMessages(library(diffcyt))

# Meta-data: marker information
# source: Bruggner et al. (2014), Table 1
# column indices of all markers, lineage markers, and functional markers
cols_markers <- c(1:19)
cols_lineage <- c(7:19)
cols_func <- setdiff(cols_markers, cols_lineage)

# channel and marker names
channel_name <- colnames(control_sample); channel_name
marker_name <- gsub("\\(.*$", "", channel_name); marker_name

# marker classes
# note: using lineage markers for 'cell type', and functional markers for 
# 'cell state'
marker_class <- rep("none", ncol(control_sample)); marker_class
marker_class[cols_lineage] <- "type"
marker_class[cols_func] <- "state"
marker_class <- factor(marker_class, levels = c("type", "state", "none")); marker_class

marker_info <- data.frame(
  channel_name, marker_name, marker_class, stringsAsFactors = FALSE
)
marker_info


colnames(df_FlowSOM)
# Create design matrix or createFormula()
design <- createDesignMatrix(df_FlowSOM, cols_design=c("clusters_flowsom")) 
# Create contrast matrix
contrast <- createContrast(c(rep(1,10))); contrast
# check
nrow(contrast) == ncol(design)
data.frame(parameters= colnames(design), contrast)

# differential testing
out_DA <- diffcyt(
  d_input = control_sample, 
  experiment_info = df_FlowSOM, 
  marker_info = marker_info, 
  design = design, 
  contrast = contrast, 
  analysis_type = "DA", 
  seed_clustering = 123
)

d_counts<-group_by(df_FlowSOM, clusters_flowsom) %>%summarise(count = n())
res_DA <- testDA_edgeR(d_counts, design, contrast)

df_FlowSOM["clusters_flowsom"]<-as.numeric(unlist(df_FlowSOM["clusters_flowsom"]))
lmFit(df_FlowSOM, design)
