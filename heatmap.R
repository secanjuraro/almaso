
############################
######## LOAD DATA #########
############################


immune_control = read.FCS("cd45pos2_control_CD45+.fcs", column.pattern = "Time", invert.pattern = TRUE)
fsc_immune <- read.flowSet("cd45pos2_control_CD45+.fcs", truncate_max_range = FALSE)


############################
###### PREPROCESSING #######
############################

names <- colnames(expr_matrix[7:19])

#Log
myTrans <- transformList(names, logTransform())
fs_log<- transform(fsc_immune, myTrans)

autoplot(fs_log, "FJComp-V500-A")

#Arcsinh

myTrans <- transformList(names, arcsinhTransform())
fs_immune_arc<- transform(fsc_immune, myTrans)


# Get expression matrix as a data frame
expr_matrix <- as.data.frame(fs_immune_arc@frames[["cd45pos2_control_CD45+.fcs"]]@exprs)
expr_matrix <- expr_matrix[1:19]



############################
####### CLUSTERING #########
############################

######### flowSOM ######

## 1. Run function FlowSOM, nClus is set to 10 by default 
flowsom <- FlowSOM(input = immune_control, 
                   transform = FALSE,
                   scale = FALSE,
                   colsToUse = c(7:19), #provide the columns for the clustering
                   nClus = 10, #we choose 14, since we also generated 14 clusters by HSNE
                   seed = 100)

## 2. Get metaclustering per cell
clusters_flowsom <- as.factor(flowsom$map$mapping[,1])
levels(clusters_flowsom) <- flowsom$metaclustering

## 3. Add flowsom clusters to dataframe
df_FlowSOM <- cbind(expr_matrix, clusters_flowsom)

##### HEATMAP #####

df_FlowSOM2 <- df_FlowSOM %>% select(-(contains("FSC") | contains("SSC"))) %>% group_by(clusters_flowsom) %>% summarise(across(everything(), mean, na.rm=TRUE))  %>% remove_rownames %>% column_to_rownames(var="clusters_flowsom")
heatmap(as.matrix(df_FlowSOM2),Rowv = NA, Colv = NA, xlab = "Marqueur", ylab="Cluster",verbose = TRUE)
