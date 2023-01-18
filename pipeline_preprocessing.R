if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("FlowSOM")

require(flowAI)
library(flowCore)
library(flowClean)
library(ggcyto)
library(flowStats)
library(PeacoQC)
library(FlowSOM)
library(devtools)
library(NOISeq)
library(edgeR)

#######SCRIPT PREPROCESSING#######

#Le pre processing contiendra une étape de compensation, une étape de transformation
#une étape de quality control et une étape de normalisation

cyto = read.FCS("VP4_Tumor_CD45+ cells.fcs", truncate_max_range = FALSE)
fsc <- read.flowSet("VP4_Tumor_CD45+ cells.fcs", truncate_max_range = FALSE)
head(fsc)

#### Compensation avec flowCore

fs_comp <- compensate(fsc, spillover(cyto))

#### Transformation with flowCore
names <- c("FJComp-APC-A", "FJComp-APC-Cy7-A", "FJComp-Alexa Fluor 700-A", "FJComp-DAPI-A", 
           "FJComp-FITC-A", "FJComp-Hoechst Red-A", "FJComp-PE-A", "FJComp-PE-Cy5-A",
           "FJComp-PE-Cy7-A", "FJComp-PE-Texas Red-A", "FJComp-PerCP-Cy5-5-A", "FJComp-V450-A", "FJComp-V500-A")

#Log
myTrans <- transformList(names, logTransform())
fs_log<- transform(fsc, myTrans)

autoplot(fs_log, "FJComp-V500-A")

#Arcsinh

myTrans <- transformList(names, arcsinhTransform())
fs_arc<- transform(fsc, myTrans)

autoplot(fs_arc, "FJComp-V500-A")

#Logicle
myTrans <- transformList(names, logicleTransform())
fs_logicle<- transform(fsc, myTrans)

autoplot(fs_logicle, "FJComp-V500-A")


####Quality control

#flowAI
resAI <- flow_auto_qc(fs_log)

#flowClean
head(fs_log)
clean_fs <- clean(fsc, vectMarkers = c(6:17),  filePrefixWithDir="sample_out", ext="fcs")
exprs(clean_fs)

clean_fs_fc <- clean_fs["GoodVsBad"< 10000]
head(clean_fs_fc)
dim(clean_fs_fc)

#peacoQC

write.flowSet(fs_arc, "cyto_comp_transf")
cyto_ct = read.FCS("cyto_comp_transf/VP4_Tumor_CD45+ cells.fcs", truncate_max_range = FALSE)

peaco_fs <- PeacoQC(cyto_ct, names, determine_good_cells="all",
        plot=20, save_fcs=TRUE, output_directory=".",
        name_directory="PeacoQC_results1", report=TRUE)

fs_peaco_qc = read.flowSet("PeacoQC_results1/fcs_files/VP4_Tumor_CD45+ cells_QC.fcs", truncate_max_range = FALSE)

#### Normalization

exp_matr <- fs_peaco_qc@frames[["VP4_Tumor_CD45+ cells_QC.fcs"]]@exprs

tmm = tmm(exp_matr, long = 1000, lc = 0, k = 0, refColumn = 1, logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10)
cpm = cpm(tmm)
