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

setwd("~/INSA/BIM 2022-2023/Projet 5BIM/real data")


############################
######## LOAD DATA #########
############################

file_name <- "cd45pos2_control_CD45+.fcs"
immune_control = read.FCS(file_name, column.pattern = "Time", invert.pattern = TRUE, truncate_max_range = FALSE)
fs_immune_control = read.flowSet(file_name, column.pattern = "Time", invert.pattern = TRUE, truncate_max_range = FALSE)


############################
###### PREPROCESSING #######
############################


  # 1 . Compensation with flowCore, if data not alraedy compensated
  if(exists("spillover(immune_control)")){
    fs_immune_control <- compensate(fs_immune_control, spillover(immune_control))
  }

  # 2. Quality control
{
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
      pca <<- prcomp(expr_mat[,first_mark:last_mark]) #Keep pca in the environment
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
    #' @description Computes a KNN graph and the Louvain method to find cells clusters 
    #' @param expr_df a data frame that corresponds to the expression matrix of an FCS file
    #' @param resolution resolution parameter for finding the clusters
    #' @return A data frame that corresponds to the expression matrix with a cluster number associated to each cell
    
    clustering_cyto <- function(df_pca,expr_df, resolution) {
      KNN_graph <- bluster::makeSNNGraph(df_pca) # Build the KNN graph for community detection
      louvain_clusters <- igraph::cluster_louvain(KNN_graph, resolution = 0.5) # Implementation of the Louvain method to find clusters 
      clusters_id <<- communities(louvain_clusters) # Get the cluster for each cell
      
      df_KNN <- data.frame(cell = numeric(),cluster = numeric()) # dataframe with each cell associated to a cluster
      
      for(i in names(clusters_id)){
        temp_df <- as.data.frame(clusters_id[i]) # get cells in each cluster as a dataframe 
        temp_df$cluster_id <- i # associate each cell to the corresponding cluster
        colnames(temp_df) <- c("cell", "cluster_id")
        df_KNN <- rbind(df_KNN,temp_df)  # bind each temporary dataframe to the previous cluster one 
      }
      clusters <- df_KNN %>% arrange(cell) %>% select(cluster_id) # order dataframe by cells and get clusters
      df_KNN <- cbind(expr_df,clusters)   # bind the cluster column for all cells to expression dataframe 
      
      return(df_KNN)
    } 
    
    #' @title Run UMAP 
    #' @name RunUMAP_cyto
    #' @description Plots the UMAP with the clusters and saves a data frame with the UMAP coordinates and the cluster associated to each cell
    #' @param df_KNN  A data frame that corresponds to the expression matrix with a cluster number associated to each cell
    #' @return A plot of the UMAP with the cells coloured by clusters 
    
    
    RunUMAP_cyto <- function(df_KNN, df_pca){
      df_umap <- df_KNN %>% select(-(contains("FSC") | contains("SSC"))) %>% select(-c(cluster_id))
      umap_m <<-as.data.frame(umap(df_umap)) 
      colnames(umap_m) <- c("UMAP1","UMAP2") 
      umap_m$cluster <- df_KNN$cluster_id
      plot_umap <- ggplot(umap_m, aes(UMAP1, UMAP2, colour = cluster)) +  geom_point() +  labs(x = "UMAP1",y = "UMAP2",subtitle = "UMAP plot")
      
      
      #df_umap <<- df_umap
      return(plot_umap)
    }
    
    #' @title Find the expression of all markers within a cluster
    #' @name findMarkers_cyto
    #' @description Computes a Wilcox test to compare the expression of all markers across all clusters. A p-value and a log2 FC are associated to each comparison 
    #' @param df_KNN  A data frame that corresponds to the expression matrix with a cluster number associated to each cell
    #' @return A data frame with all the comparisons made associated to a p value and a FC 
    
    findMarkers_cyto <- function(df_KNN) {
      df_W <<- df_KNN %>% select(-(contains("FSC") | contains("SSC"))) %>% group_by(cluster_id) %>% arrange(as.numeric(cluster_id)) #Get data frame for the wilcox test. Delete columns corresponding to the FSC and SSC
      df_FC <<- df_W %>% select(-(contains("FSC") | contains("SSC"))) %>% group_by(cluster_id) %>% summarise(across(everything(),mean)) %>% select(-cluster_id) #Get data frame to calculate the fold change (average expression level for each luster)
      
      markers <- colnames(df_FC) # Get markers from df_FC (does not contain column cluster_id)
      clusters_id <- order(levels(factor(df_W$cluster_id))) # Get the number of clusters and order them 
      
      df_markers <- data.frame(marker = character(),cluster_a = numeric(), cluster_b = numeric(), pvalue = numeric(),FC = numeric() ) #Initialize dataframe with all information regarding the expression of markers in the clusters 
      for(i in markers){   # i browses each marker
        for(j in 1:(length(clusters_id)-1)){    # j browses each cluster
          for(k in (j+1):length(clusters_id)){  # once a cluster chosen with j, k browses the rest of the clusters. We do not compare a marker within the same cluster 
            cluster_a <- df_W %>% filter(cluster_id == j) %>% pull(i)   # get expression values of i for j
            cluster_b <- df_W %>% filter(cluster_id == k) %>% pull(i)   # get expression values of i for k
            pvalue <- wilcox.test(cluster_a,cluster_b)$p.value          # get pvalue from wilcox test 
            fold_change <- log2(abs((df_FC[[i]][j])/(df_FC[[i]][[k]]))) # calculate log2 fold change with the average expression of i in j and k 
            info <- c(i,j,k,pvalue,fold_change)  # create vector with all informations 
            df_markers <- rbind(df_markers,info)   # fill the df_Wilcox
          } # close k loop
        } #close j loop
      } #close i loop
      
      colnames(df_markers) <- c("marker", "cluster_a","cluster_b","p_value","FC") 
      df_markers$adj_pvalue <- p.adjust(df_markers$p_value,"BH") #adjust p_value 
      df_markers$cluster_a <- as.numeric(df_markers$cluster_a) 
      df_markers$cluster_b <- as.numeric(df_markers$cluster_b)
      df_markers <- df_markers %>% group_by(cluster_a,marker) %>% arrange(adj_pvalue,.by_group = TRUE)  #arrange dataframe by pvalue 
      
      
      return(df_markers)
      
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
      #' @return A heatmap containing the mean expression of each marker in each cluster
      
      
      heatmap_cyto <- function(df_clust) {
        if ('cluster_id' %in% colnames(df_clust)){
          #Create new data frame with mean expression of markers by cluster
          df_clust_mean <- df_clust %>% select(-(contains("FSC") | contains("SSC"))) %>% group_by(cluster_id) %>% summarise(across(everything(), mean, na.rm=TRUE))  %>% remove_rownames %>% column_to_rownames(var="cluster_id")
          #Create heatmap
          heat <- heatmap(as.matrix(df_clust_mean),Rowv = NA, Colv = NA, xlab = "Marqueur", ylab="Cluster",verbose = TRUE)
        }else {
          stop('Clusters  associated to each cells must be contained in column "clusters_id')
        }
        return (heat)
      }    
      
      
      #' @title Visualize the expression of a marker
      #' @name marker_expression
      #' @description Plots the UMAP with a color gradient that correspond to the expression of a marker 
      #' @param df_KNN  A data frame that corresponds to the expression matrix with a cluster number associated to each cell
      #' @param marker  A string with the name of the marker
      #' @param df_umap A dataframe with the coordinates of the umap for all cells 
      #' @return A plot of the UMAP with the cells coloured by clusters 
      
      marker_expression <- function(df_KNN,marker,umap_m){
        marker <- paste0(marker)
        marker_value <- df_KNN[[marker]]
        mid<-mean(marker_value)
        plot_marker <- ggplot(umap_m, aes(V1,V2, colour = as.numeric(marker_value))) +  geom_point() +  labs(x = "UMAP1",y = "UMAP2",subtitle = "UMAP plot ")+scale_color_gradient2(midpoint=mid, low="blue", mid="white", high="red", space ="Lab" )
        
        return(plot_marker)
        
      }
    }
    

}
  
    
#######################
###### TESTS ##########
#######################
    

fsc_peaco_QC <- peaco_QC(fs_immune_control, file_name)

exp_matr <- fsc_peaco_QC@frames[[file_name]]@exprs
exp_matr_tmm = tmm(exp_matr, long = 1000, lc = 0, k = 0, refColumn = 1, logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10)
exp_matr_cpm = cpm(exp_matr_tmm)

PCA(exp_matr_cpm, 7, 19)
summary(pca)

df_pca <- choose_dims_PCA(pca, 5)

df_KNN <- clustering_cyto(df_pca,exp_matr_cpm)

plt_umap <- RunUMAP_cyto(df_KNN, df_pca)
plt_umap

df_wilcox <- findMarkers_cyto(df_KNN)

test1 <- df_wilcox %>% filter(adj_pvalue < 0.05 & FC > abs(2.5))
write.csv(test1, file = 'df_wilcox_pca.csv')


heatmap_cyto(df_KNN)

plt_marker <- marker_expression(df_KNN,"FSC-H",umap_m)
plt_marker
head(umap_m)
