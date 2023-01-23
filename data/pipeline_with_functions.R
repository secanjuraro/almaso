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
  library(ggfortify)
  library(factoextra)
  library(flowAI)
  library(ggcyto)
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

# Get expression matrix as a data frame
expr_matrix <- as.data.frame(immune_control@exprs)

############################
###### PREPROCESSING #######
############################


# 1 . Compensation with flowCore, if data not alraedy compensated
fs_comp <- compensate(fs_immune_control, spillover(cyto))

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
    fsc_peaco_qc = read.flowSet(path, truncate_max_range = FALSE) # Read the result stored in path as a flowSet
    
    return(fsc_peaco_qc)
  }    
}

fsc_peaco_QC <- peaco_QC(fs_immune_control, file_name)

# 3. Scaling normalization

exp_matr <- fsc_peaco_QC@frames[[file_name]]@exprs

exp_matr_tmm = tmm(exp_matr, long = 1000, lc = 0, k = 0, refColumn = 1, logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10)
exp_matr_cpm = cpm(exp_matr_tmm)

##################################
###### DIMENSION REDUCTION #######
##################################

{
  #' @title Perform PCA on cells
  #' @name PCA
  #' @description Perform PCA on chosen markers to keep reduced number of dimensions
  #' @param expr_mat a matrix or data frame that corresponds to the expression matrix of an FCS file
  #' @param first_mark a column number of the expression matrix
  #' @param last_mark a column number of the expression matrix
  #' @return A prcomp PCA object
  
  PCA <- function(expr_mat, first_mark, last_mark) {
    if (first_mark > 0 && last_mark > 0 && first_mark < dim(expr_mat)[2] && last_mark < dim(expr_mat)[2]){
      pca <<- prcomp(expr_matrix[first_mark:last_mark]) #Keep pca in the environment
      viz <- fviz_eig(pca) #Create visualization
    }else{
      stop("First and last markers numbers do not fit expression matrix dimensions")
    }
    return(viz)
  }    
}

PCA(exp_matr_cpm, 7, 19)
summary(pca)

{
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
}

df_pca <- choose_dims_PCA(pca, 5)
head(df_pca)

############################
####### CLUSTERING #########
############################






############################
###### VISUALIZATION #######
############################

#Heatmap

{
  #' @title Make heatmap of marker expression in clusters
  #' @name  heatmap 
  #' @description Build heatmap of mean expression of each marker in each cluster
  #' @param df_clust a data frame containing the expression matrix and the clusters associated to each cell in the column clusters_id
  #' @return A heatmap containing the mean expression of each marker in each cluster
  
  
  heatmap <- function(df_clust) {
    if ('clusters_id' %in% colnames(df_clust)){
      #Create new data frame with mean expression of markers by cluster
      df_clust_mean <- df_clust %>% select(-(contains("FSC") | contains("SSC"))) %>% group_by(clusters_id) %>% summarise(across(everything(), mean, na.rm=TRUE))  %>% remove_rownames %>% column_to_rownames(var="clusters_id")
      #Create heatmap
      heat <- heatmap(as.matrix(df_clust_mean),Rowv = NA, Colv = NA, xlab = "Marqueur", ylab="Cluster",verbose = TRUE)
    }else {
      stop('Clusters  associated to each cells must be contained in column "clusters_id')
    }
    return (heat)
  }    
}
