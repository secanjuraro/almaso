# almaso

### Description
This R project proposes a pipeline to analyze flow cytometry data with a scRNA-Seq approach. 

In this repository you will find:

*	A text document *Projet_almaso_recap.pdf* that reports all the packages we tried, and all the choices we made during the development of the pipeline.
*	A *Pipeline_scRNA-Seq_Seurat.R* file with a basic Seurat pipeline to analyse scRNA-Seq data.
*	A R project directory named “*Flow_Cytometry_Project_INSA_almaso/*” composed of:
       *	A *renv/* directory : it is a reproducable environment. It stores all packages with the right version. A snapshot of all package version is available in the file *renv.lock*. 
       *	A *almaso_pipeline.R* file that stores all functions composing our pipeline. At the end of the file, a test code is written to run all functions once.
       *    A directory  *data/* that contains all the flow cytometry data files (.fcs files). 
       *	A *Vignette_almaso_Rmarkdown.Rmd* file and its HTML knit version: it is a tutorial on how to use the pipeline we created. It has cache and file repositories. 
       *	The *PeacoQC_result/* repository is the result of the function peaco_QC run in the vignette. 
       *	The *df_wilcox_pca.csv* file is the result of the last function (run in the vignette). It is table of all pvalue and log Fold Change for all clusters against all markers for the data file used in the vignette.


A tutorial can be found at: [vignette](almaso/R/Vignette_almaso_Rmarkdown.Rmd) the file Rmarkdown or HTML file at almaso/R/Vignette_almaso_Rmarkdown.html
