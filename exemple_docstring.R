{
  #' @title Cluster cells by markers 
  #' @name clustering
  #' @description Computes a KNN graph and the Louvain method to find cells clusters 
  #' @param expr_df a data frame that corresponds to the expression matrix of an FCS file
  #' @param resolution resolution parameter for finding the clusters
  #' @return A data frame that corresponds to the expression matrix with a cluster number associated to each cell
  
  clustering <- function(expr_df,resolution) {
    KNN_graph <- bluster::makeSNNGraph(d) # Build the KNN graph for community detection
    louvain_clusters <- igraph::cluster_louvain(graph, resolution = 0.5) # Implementation of 
    clusters_id <- communities(louvain_clusters) 
    return(expr_df)
  }    
}