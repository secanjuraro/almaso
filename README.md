# almaso

### Description
This R project proposes a pipeline to analyze flow cytometry data with a scRNA-Seq approach. 

A tutorial can be found at: [vignette](almaso/R/Vignette_almaso_Rmarkdown.Rmd) the file Rmarkdown or HTML file at almaso/R/Vignette_almaso_Rmarkdown.html


The .R file *almaso_pipeline.R* gather all functions composing our pipeline. 

The .R file *Seurat.R* is a basic Seurat pipeline to analyse scRNA-Seq data.

The directory *data* contains all the flow cytometry data files (.fcs files). 

The .dox file *Project_almaso_recap.dox* reports all the packages we tried and all choices made during the development of the pipeline. 

The directory *almaso* is the directory of the package we created to keep track of all packages necessary. Named almaso, it contains in the almaso/R/ directory the vignette and the same exact *almaso_pipeline.R* file (necessary for the functions in the vignette to be sourced). 
