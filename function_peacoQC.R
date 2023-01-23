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
    
    return(fsc_peaco_QC)
  }    
}

# A améliorer ci-besoin:
# rajouter un if(flowFrame = flowSet)
# prend-t-on les colonnes SSC et FSc?? Là j'ai tout pris... 