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
  library(factoextra)
  library(flowAI)
  library(PeacoQC)
  library(NOISeq)
  library(edgeR)
}

############################
### SET WORKING DIR      ###
############################

# setwd("~/INSA/BIM 2022-2023/Projet 5BIM/real data")


############################
######## LOAD DATA #########
############################

file_name <- "Vp6_Control_CD45+ cells.fcs"
path <- file.path("data", file_name); path
immune_control = read.FCS(path, column.pattern = "Time", invert.pattern = TRUE, truncate_max_range = FALSE)
fs_immune_control = read.flowSet(path, column.pattern = "Time", invert.pattern = TRUE, truncate_max_range = FALSE)

# file_name <- "cd45pos2_control_CD45+.fcs"
# path <- file.path("Vrais_donnees_cyto", file_name); path
# immune_control = read.FCS(path, column.pattern = "Time", invert.pattern = TRUE, truncate_max_range = FALSE)
# fs_immune_control = read.flowSet(path, column.pattern = "Time", invert.pattern = TRUE, truncate_max_range = FALSE)


{
  
  
  #################################
  ######   PRE PROCESSING   #######
  #################################
  
  #' @title Apply PeacoQC analysis to a flowSet
  #' @name peaco_QC
  #' @description Apply the function PeacoQC from the package PeacoQC to a flowSet to determine peaks and to remove anomalies caused by e.g. clogs, changes in speed etc.
  #' @param fsc a flowSet on which will be applied the function PeacoQC
  #' @param file_name a string that holds the name of the fsc file -- not the path
  #' @return A cleaned flowSet
  
  peaco_QC <- function(fsc, file_name) {
    
    channels <- colnames(fsc)  # channels on which the QC will be applied
    
    # Apply PeacoQC function from the PeacoQC package on a flowFrame from the flowSet fsc
    # For more informations, please refer to the PeacoQC docstring: https://bioconductor.org/packages/release/bioc/manuals/PeacoQC/man/PeacoQC.pdf
    peaco_res <- PeacoQC(fsc[[1]], channels, determine_good_cells="all",
                         save_fcs=TRUE, suffix_fcs = "", output_directory=".",
                         name_directory="PeacoQC_results", report=TRUE)
    
    path <- file.path("PeacoQC_results/fcs_files/", file_name)    # The result is stored under the name of the fsc => file_name
    fsc_peaco_qc = read.flowSet(path, truncate_max_range = FALSE ) # Read the result stored in path as a flowSet
    
    return(fsc_peaco_qc)
  }
  
  #' @title Normalize flow cytometry expression matrix
  #' @name norm_cyto
  #' @description Transform negative values into zeros, keep only useful columns, do logCPM normalization
  #' @param exp_matr an expression matrix after QC
  #' @return A normalized expression matrix
  
  norm_cyto <- function(exp_matr) {
    
    #Attribute zero to negative values
    exp_matr[exp_matr <0 ] <- 0
    #Remove column Original_ID that was added by PeacoQC
    exp_matr <- as.data.frame(exp_matr) %>% select(-('Original_ID'))
    #Perform logCPM transformation
    exp_matr_cpm <- cpm(exp_matr)
    
    return(exp_matr_cpm)
  }
  
  ##################################
  ###### DIMENSION REDUCTION #######
  ##################################
  
  
  #' @title Perform PCA on cells
  #' @name PCA
  #' @description Perform PCA on chosen markers to keep reduced number of dimensions
  #' @param expr_mat a matrix or data frame that corresponds to the expression matrix of an FCS file
  #' @param first_mark a column number of the expression matrix
  #' @param last_mark a column number of the expression matrix
  #' @return A prcomp PCA object
  
  PCA <- function(expr_mat, first_mark, last_mark) {
    if (first_mark > 0 && last_mark > 0 && first_mark <= dim(expr_mat)[2] && last_mark <= dim(expr_mat)[2]){
      pca <<- prcomp(expr_mat[,first_mark:last_mark], scale = TRUE) #Keep pca in the environment
      viz <- fviz_eig(pca) #Create visualization
    }else{
      stop("First and last markers numbers do not fit expression matrix dimensions")
    }
    return(viz)
  }
  
  #' @title Choose number of PC axis to use for clustering
  #' @name  choose_dims_PCA
  #' @description Build data frame with chosen number of PC axes (minimum 2 axes)
  #' @param pca a PCA object returned by the PCA function
  #' @param number_axes the number of chosen PC axes, between 2 and the maximum number of PC axes
  #' @return A data frame containing the PC axes
  
  choose_dims_PCA <- function(pca, number_axes) {
    #Create new data frame with 2 PC dimensions
    pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
    if (number_axes > 2 && number_axes < dim(pca$rotation)[1]){
      #Add columns to data frame corresponding to chosen PC axes
      for (i in (3:number_axes)){
        pca_df <- cbind(pca_df,PC = pca$x[,i]) #Add columns
        colnames(pca_df)[i] <- paste0('PC',as.character(i)) #Set column names
      }
    }else{
      print('Number of axes was incorrect, 2 PC axes were kept')
    }
    return (pca_df)
  }
  
  ############################
  ####### CLUSTERING #########
  ############################
  
  #' @title Cluster cells by markers
  #' @name clustering_cyto
  #' @description Computes a SNN graph and the Louvain method to find cells clusters
  #' @param expr_df a data frame that corresponds to the expression matrix of an FCS file
  #' @param resolution resolution parameter for finding the clusters
  #' @return A data frame that corresponds to the expression matrix with a cluster number associated to each cell
  
  clustering_cyto <- function(df_pca,expr_df, resolution) {
    SNN_graph <- bluster::makeSNNGraph(df_pca) # Build the SNN graph for community detection
    louvain_clusters <- igraph::cluster_louvain(SNN_graph, resolution = 0.5) # Implementation of the Louvain method to find clusters
    clusters_id <<- communities(louvain_clusters) # Get the cluster for each cell
    
    df_SNN <- data.frame(cell = numeric(),cluster = numeric()) # dataframe with each cell associated to a cluster
    
    for(i in names(clusters_id)){
      temp_df <- as.data.frame(clusters_id[i]) # get cells in each cluster as a dataframe
      temp_df$cluster_id <- i # associate each cell to the corresponding cluster
      colnames(temp_df) <- c("cell", "cluster_id")
      df_SNN <- rbind(df_SNN,temp_df)  # bind each temporary dataframe to the previous cluster one
    }
    clusters <- df_SNN %>% arrange(cell) %>% select(cluster_id) # order dataframe by cells and get clusters
    df_SNN <- cbind(expr_df,clusters)   # bind the cluster column for all cells to expression dataframe
    
    return(df_SNN)
  }
  
  #' @title Run UMAP
  #' @name RunUMAP_cyto
  #' @description Plots the UMAP with the clusters and saves a data frame with the UMAP coordinates and the cluster associated to each cell
  #' @param df_SNN  A data frame that corresponds to the expression matrix with a cluster number associated to each cell
  #' @return A plot of the UMAP with the cells coloured by clusters
  
  
  RunUMAP_cyto <- function(df_SNN, df_pca){
    umap_m <<-as.data.frame(umap(df_pca))
    colnames(umap_m) <- c("UMAP1","UMAP2")
    umap_m$cluster <- as.numeric(df_SNN$cluster_id)
    plot_umap <- ggplot(umap_m, aes(UMAP1, UMAP2, colour = factor(cluster))) +  geom_point() +  labs(x = "UMAP1",y = "UMAP2",subtitle = "UMAP plot")
    
    umap_m <<- umap_m
    return(plot_umap)
  }
  
  #######################################
  ####### DIFFERENTIAL EXPRESSION #######
  #######################################
  
  
  #' @title Find the expression of all markers within a cluster
  #' @name findMarkers_cyto
  #' @description Computes a Wilcox test to compare the expression of all markers across all clusters. A p-value and a log2 FC are associated to each comparison 
  #' @param df_SNN  A data frame that corresponds to the expression matrix with a cluster number associated to each cell
  #' @return A data frame with all the comparisons made associated to a p value and a FC 
  
  findMarkers_cyto <- function(df_SNN) {
    df_W <- df_SNN %>% select(-(contains("FSC") | contains("SSC"))) %>% group_by(cluster_id) %>% arrange(as.numeric(cluster_id)) #Get data frame for the wilcox test. Delete columns corresponding to the FSC and SSC
    markers <- colnames(df_SNN %>% select(-(contains("FSC") | contains("SSC"))) %>% select(-cluster_id))  # Get markers from df_FC (does not contain column cluster_id)
    df_FC <- df_W %>% select(-(contains("FSC") | contains("SSC"))) %>% group_by(cluster_id) %>% summarise(across(everything(),mean)) #%>% select(-cluster_id) #Get data frame to calculate the fold change (average expression level for each luster)
    
    clusters_id <- order(levels(factor(df_W$cluster_id))) # Get the number of clusters and order them 
    
    df_markers <- data.frame(marker = character(),cluster = numeric(),pvalue = numeric(),FC = numeric() ) #Initialize dataframe with all information regarding the expression of markers in the clusters 
    
    for(i in markers){ # i browses each marker 
      for(j in 1:(length(clusters_id))){ # j browses each cluster 
        df_W_b <- df_W %>% filter(cluster_id != j) # filtered data frame without the cluster selected by j for the wilcoxon test
        df_FC_b <- df_FC %>% filter(cluster_id!=j) %>% select(-cluster_id) %>% summarise(across(everything(),mean)) # filtered data frame without the cluster selected by j for the fold change 
        cluster_a <- df_W %>% filter(cluster_id == j) %>% pull(i) # get expression values of i for j
        cluster_b <- numeric() # create vector to keep the expressions values of i for all clusters except i 
        for(k in 1:(length(clusters_id)-1)){ # k browses all clusters different from i 
          cluster_b <- append(cluster_b,df_W_b %>% filter(cluster_id == k ) %>% pull(i)) # get expression values of i for the rest of clusters
        } # close k loop
        pvalue <- wilcox.test(cluster_a,cluster_b)$p.value  # get pvalue from wilcox test 
        fold_change <- log2(abs(df_FC[[i]][j])/df_FC_b[[i]])  # calculate log2 fold change with the average expression of i in j and the rest of clusters
        info <- c(i,j,pvalue,fold_change) 
        df_markers <- rbind(df_markers,info) 
      } # close j loop
    } # close i loop
    colnames(df_markers) <- c("marker", "cluster","p_value","FC") 
    df_markers$adj_pvalue <- p.adjust(df_markers$p_value,"BH") #adjust p_value 
    df_markers$cluster <- as.numeric(df_markers$cluster) 
    df_markers <- df_markers %>% group_by(cluster,marker) %>% arrange(adj_pvalue,.by_group = TRUE)  #arrange dataframe by pvalue 
    
    return(df_markers)
  }
  
  #' @title Add the name of proteins associated to the markers in alphabetical order
  #' @name add_ProtMarkerNames
  #' @description Get protein names, order by markers in alphabetical order and bind protein names to dataframe obtained after wilcox test
  #' @param exp_matrix A matrix corresponding to the expression matrix
  #' @param df_wilcox  A data frame returned by the findMarkers_cyto function
  #' @param fsc A flowSet
  #' @return A data frame with a column of protein names associated to markers
  
  add_ProtMarkerNames <- function(exp_matr, df_wilcox, fsc){
    
    mark <- as.data.frame(colnames(exp_matr))
    colnames(mark)<- "Marker"
    mark <- mark %>% filter(!(grepl('FSC|SSC|Original', Marker))); mark
    
    antigene_list <- as.vector(fsc@frames[[file_name]]@parameters@data[["desc"]])
    antigene_list <-antigene_list[7:length(antigene_list)]
    mark <- cbind(mark, antigene_list)
    colnames(mark)[2] <- "Antigene"
    
    mark <- mark %>% arrange(Marker)
    antigene_sort <- unname(unlist(as.data.frame(mark['Antigene'])))
    
    antigene <- rep(antigene_sort, length(unique(as.numeric(df_SNN$cluster_id))))
    df_wilcox_antigene <- cbind(df_wilcox, antigene)
    colnames(df_wilcox_antigene)[6] <- "Antigene"
    
    
    return (df_wilcox_antigene)
  }
  
  ############################
  ###### VISUALIZATION #######
  ############################
  
  #Heatmap
  
  {
    #' @title Make heatmap of marker expression in clusters
    #' @name  heatmap_cyto
    #' @description Build heatmap of mean expression of each marker in each cluster
    #' @param df_clust a data frame containing the expression matrix and the clusters associated to each cell in the column clusters_id
    #' @param fsc a flowSet
    #' @param file_name name of the flowSet
    #' @return A heatmap containing the mean expression of each marker in each cluster
    
    
    heatmap_cyto <- function(df_clust, fsc, file_name) {
      if ('cluster_id' %in% colnames(df_clust)){
        #Create new data frame with mean expression of markers by cluster
        df_clust_mean <- df_clust %>% select(-(contains("FSC") | contains("SSC"))) %>% group_by(cluster_id) %>% summarise(across(everything(), mean, na.rm=TRUE))  #%>% remove_rownames %>% column_to_rownames(var="cluster_id")
  
        #Scale values 
        for( i in colnames(df_clust_mean %>% select(-c(cluster_id)))){
          df_clust_mean[[i]] <- scale(df_clust_mean[[i]])
        }
        
        #Change marker names to antigene names
        antigene_list <- as.vector(fsc@frames[[file_name]]@parameters@data[["desc"]]) # we take the gene names column stored in the "desc" variable
        antigene_list <-antigene_list[7:length(antigene_list)] # Only the columns 7 to 19 have gene associated
        colnames(df_clust_mean)[2:14] <- antigene_list
        
        #Melt dataframe to use ggplot
        melt_df_clust <<- reshape2::melt(df_clust_mean) 
        colnames(melt_df_clust) <- c("cluster","antigene","value")
        melt_df_clust$cluster <- as.numeric(melt_df_clust$cluster)
        
        #Heatmap with ggplot
        heat <- ggplot(melt_df_clust,aes(factor(cluster),antigene)) + geom_tile(aes(fill=value)) + scale_fill_gradient(low = "gold",high= "red") + ggtitle("Antigene expression within each cluster")
      }else {
        stop('Clusters  associated to each cells must be contained in column "clusters_id')
      }
      return (heat)
    }
    
    
    #' @title Visualize the expression of a marker
    #' @name marker_expression
    #' @description Plots the UMAP with a color gradient that correspond to the expression of a marker
    #' @param df_SNN  A data frame that corresponds to the expression matrix with a cluster number associated to each cell
    #' @param protein  A string with the name of the protein
    #' @param marker  A string with the name of the fluorescent marker that stained for the protein
    #' @param df_umap A dataframe with the coordinates of the umap for all cells
    #' @return A plot of the UMAP with the cells coloured by clusters
    
    marker_expression <- function(df_SNN,protein,marker,umap_m){
      marker <- paste0(marker)
      protein <- paste0(protein)
      marker_value <- df_SNN[[marker]]
      mid<-mean(marker_value)
      title <- paste(" Expression of ", protein, sep ="")
      m <- as.numeric(marker_value)
      plot_marker <- ggplot(umap_m, aes(UMAP1,UMAP2, colour = m )) +  geom_point() +  labs(x = "UMAP1",y = "UMAP2",subtitle = title ) + scale_color_gradient2(midpoint=mid, low="blue", mid="white", high="red", space ="Lab" )
      return(plot_marker)
      
    }
    
    #' @title Visualize the expression of all markers on different plots
    #' @name all_markers_expression
    #' @description Use the marker_expression function on each marker and plot all expression maps
    #' @param df_SNN  A data frame that corresponds to the expression matrix with a cluster number associated to each cell
    #' @param df_umap A dataframe with the coordinates of the umap for all cells
    #' @param fsc a flowSet
    #' @param file_name name of the flowSet
    #' @return A list of ggplots
    
    all_markers_expression <- function(df_SNN,umap_m, fsc, file_name){
      myplots <- list()  # new empty list
      
      mark <- as.data.frame(colnames(df_SNN)) # dataframe of all markers names
      colnames(mark)<- "Marker"
      mark <- mark %>% filter(!(grepl('FSC|SSC|Original|cluster', Marker))); mark  # we only take the fluorescent markers
      
      antigene_list <- as.vector(fsc@frames[[file_name]]@parameters@data[["desc"]])# we take the gene names column stored in the "desc" variable
      antigene_list <-antigene_list[7:length(antigene_list)] # Only the columns 7 to 19 have gene associated
      
      mark <- cbind(mark, antigene_list) # we bind both columns together
      colnames(mark)[2] <- "Antigene"
      
      for (i in antigene_list) {
        m  <- mark %>% filter(Antigene == i ) %>% select(Marker) # we select the marker name associated to antigene i 
        m  <- paste0(m)
        p1 <- marker_expression(df_SNN, i, m, umap_m)  # plot the expression marker map for the marker i
        myplots[[i]] <- p1  # add each plot into plot list
      }
      grid.arrange(grobs = myplots, ncol = 5) # plot all the plots all together
      return(myplots)   # return the list of plots
    }
    
  }
  
}


#######################
###### TESTS ##########
#######################


# Compensation with flowCore, if data not alraedy compensated
if(exists("spillover(immune_control)")){
  fs_immune_control <- compensate(fs_immune_control, spillover(immune_control))
}

# Quality control
fsc_peaco_QC <- peaco_QC(fs_immune_control, file_name)
exp_matr <- fsc_peaco_QC@frames[[file_name]]@exprs

# Normalisation
exp_matr_cpm <- norm_cyto(exp_matr)
head(exp_matr_cpm)

# Perform PCA
PCA(exp_matr_cpm, 7, 19)
df_pca <- choose_dims_PCA(pca, 5)

# PCA Visualisation
summary(pca)
ggplot(data = df_pca, aes_string(x = "PC1", y = "PC2")) + geom_point(size = 3, color = 'orange')
fviz_contrib(pca, "var")

# Perform Clustering
df_SNN <- clustering_cyto(df_pca,exp_matr_cpm)

# Clustering Visualization
plt_umap <- RunUMAP_cyto(df_SNN, df_pca)
plt_umap

# Visualize markers in clusters
## Heatmap:
heatmap_cyto(df_SNN,fs_immune_control,file_name)
## Visualization of 1 marker on the UMAP: FJComp-APC-A <=> SiglecH
plt_marker <- marker_expression(df_SNN,"SiglecH","FJComp-APC-A",umap_m); plt_marker
## Visualization of all markers on separated UMAP:
myplots <- all_markers_expression(df_SNN,umap_m,fs_immune_control,file_name)

# Perform differential analysis
df_wilcox <- findMarkers_cyto(df_SNN)
df_wilcox <- add_ProtMarkerNames(exp_matr, df_wilcox, fs_immune_control)
head(df_wilcox)

test1 <- df_wilcox %>% filter(adj_pvalue < 0.05 & FC > 0)
write.csv(test1, file = 'df_wilcox_pca.csv')