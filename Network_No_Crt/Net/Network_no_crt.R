# Try Diff - 0108
##===============================================================================##
##                                0.Data pre                                     ## 
##===============================================================================##
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA")
source("Function_Source.R")
sample_index =read_excel("Samples_RNA-Seq.xlsx")
control_index = dplyr::filter(sample_index,TRT == "a") %>% dplyr::select('Tube ID') %>% unlist(use.names = F)
treatment_index = dplyr::filter(sample_index,TRT == "b") %>%  dplyr::select('Tube ID') %>% unlist(use.names = F)

setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/counts/")
#test = read.csv("M6228-1.counts.raw.txt",sep = "\t",header = F) %>% dplyr::slice(1:(n()-5))
#test2 = read.csv("M6228-2.counts.raw.txt",sep = "\t",header = F) %>% dplyr::slice(1:(n()-5))
#sum(test$V2) > sum(test2$ V2) # use 1
# pre index 
raw_data_index = list.files()
raw_data_index = raw_data_index[-c(1,which(raw_data_index == "M6228-2.counts.raw.txt"),length(raw_data_index))]

# pre raw data - ungrouped
data_expr_all_raw = data.frame()
for (i in seq_along(raw_data_index)){
  tmp_name = substr(raw_data_index[i],1,5)
  #tmp_rownames = read.csv(raw_data_index[i],sep = "\t",header = F) %>% dplyr::slice(1:(n()-5)) %>%  dplyr::select(V1)
  tmp_read = read.csv(raw_data_index[i],sep = "\t",header = F) %>% dplyr::slice(1:(n()-5)) %>% column_to_rownames('V1')
  if (i == 1){
    data_expr_all_raw = data.frame(tmp_name = tmp_read)
    colnames(data_expr_all_raw) = tmp_name
  } else {
    data_expr_all_raw = cbind(data_expr_all_raw,tmp_read)
    colnames(data_expr_all_raw)[i] = tmp_name
  }
}
data_expr_all_with0 = data_expr_all_raw[,c(which(substr(raw_data_index,2,5) %in% control_index),which(substr(raw_data_index,2,5) %in% treatment_index))]
dim(data_expr_all_with0)

##########################################################################
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net")
networkData_final  =  DataPre_C(data_expr_all_with0, cousin = 0.4, n1 = 9, n2 = 10, perct = 0.5,
                                thres_rmzero = 5,count_rmzero = 9,Correct=F)
network_final = data.frame(networkData_final[[1]])
#dim(datExpr_treatment)
datExpr_control = t(network_final[,which(substr(names(network_final),2,5) %in% control_index)])
datExpr_treatment = t(network_final[,which(substr(names(network_final),2,5) %in% treatment_index)])

# we have data pre ready!!!
# load("data_expr_allprepare_with_corrections_top50.RData")
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
##===============================================================================##
##                                1.Weighted                                     ## 
##===============================================================================##
## pick soft thresholds
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft_b_cl = pickSoftThreshold(datExpr_control, powerVector = powers, corFnc = "bicor",verbose = 0)
### 10 works good. 10 - 0.832 and corresponging mean connectivity
#softPower_b = min(sft_b_cl$fitIndices[,1][which(sft_b_cl$fitIndices[,2] > 0.9)])
# pre_checked
softPower_b = sft_b_cl$powerEstimate
MeanK_b = sft_b_cl$fitIndices[softPower_b,5]
# Plot the results of threshold picking:
pdf(file = "soft_b_threshold.pdf", width = 12, height = 9)
sizeGrWindow(9,5);cex1 = 0.9;par(mfrow = c(1,2))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft_b_cl$fitIndices[,1], -sign(sft_b_cl$fitIndices[,3])*sft_b_cl$fitIndices[,2],
     xlab="Soft Threshold bicor (power) control",ylab="Scale Free Topology Model Fit bicor, unsigned R^2 control",type="n",
     main = paste("Scale independence bicor control"));text(sft_b_cl$fitIndices[,1], -sign(sft_b_cl$fitIndices[,3])*sft_b_cl$fitIndices[,2],
                                                            labels=powers,cex=cex1,col="red");abline(h=0.80,col="red")
# this line corresponds to using an R^2 cut-off of h
# Mean connectivity as a function of the soft-thresholding power
plot(sft_b_cl$fitIndices[,1], sft_b_cl$fitIndices[,5],
     xlab="Soft Threshold bicor (power) control",ylab="Mean Connectivity bicor control",type="n",
     main = paste("Mean connectivity bicor control"));text(sft_b_cl$fitIndices[,1], sft_b_cl$fitIndices[,5], labels=powers, cex=cex1,col="red");abline(h=MeanK_b,col="red")
dev.off()
#save
save(sft_b_cl,softPower_b,MeanK_b,file = "SoftThres_bicor_control.RData")
load("SoftThres_bicor_control.RData")
print("Step2 - soft thre plotted and Rdata saved")
#======================================================================================
#                           2 . soft-threshold and dissimilarity                ######
#=======================================================================================
# ref - cl; test - ht
# dim(datExpr14_ht);dim(datExpr14_cl)
# Peaerson Cor
adjacency_control = adjacency(datExpr_control,power=softPower_b,type="unsigned",corFnc = "bicor");
diag(adjacency_control)=0
dissTOM_control = 1-TOMsimilarity(adjacency_control, TOMType="unsigned")
geneTree_control = hclust(as.dist(dissTOM_control), method ="average")
#
adjacency_treatment = adjacency(datExpr_treatment,power=softPower_b,type="unsigned",corFnc = "bicor");
diag(adjacency_treatment)=0
dissTOM_treatment = 1-TOMsimilarity(adjacency_treatment, TOMType="unsigned")
geneTree_treatment = hclust(as.dist(dissTOM_treatment), method="average")
# save the matrix
save(adjacency_control,dissTOM_control,geneTree_control,
     adjacency_treatment,dissTOM_treatment,geneTree_treatment,file = "AllMatrix.RData")
print("Step3 - adj matrix created and rdata saved")
#================================================================ =======================
#                                 3.  plot trees                                  ######
#========================================================================================
# cor
pdf("dendrogram_bicor.pdf",     height=6,width=16)
par(mfrow=c(1,2))
plot(geneTree_control,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity bicor (control))",labels=FALSE,hang=0.04);
plot(geneTree_treatment,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity bicor (treatment)",labels=FALSE,hang=0.04);
dev.off()
print("Step4 - dissmi plottd and rdata saved")
#=========================================================================================
#                                4.cutting and merging                              ######
#=========================================================================================
# set the minimum module size relatively high:
minModuleSize = 20;
# Module identification using dynamic tree cut:
dynamicMods_control = cutreeDynamic(dendro = geneTree_control, 
                                    distM = dissTOM_control,
                                    cutHeight=0.995, deepSplit = 1, 
                                    #pamRespectsDendro = FALSE,
                                    minClusterSize = minModuleSize);
table(dynamicMods_control)
# Convert numeric lables into colors
dynamicColors_control = labels2colors(dynamicMods_control)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,16)
pdf("dendrogram_bicor_control_nomg.pdf",height=8,width=16)
plotDendroAndColors(geneTree_control, 
                    dynamicColors_control, 
                    "Dynamic_Tree_Cut_bicor_control_nomg",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors bicor control nomg")
dev.off()
print("Step5 - cutting finished")
### Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList_control = moduleEigengenes(datExpr_control,colors = dynamicColors_control)
MEs_control = MEList_control$eigengenes
#greyMEName = paste(moduleColor.getMEprefix(), "grey", sep = "") 
#if (greyMEName %in% colnames(MEList14_b_cl$eigengenes))  { print("grey found")
#  MEs14_b_cl = removeGreyME(MEList14_b_cl$eigengenes)}
# Calculate dissimilarity of module eigengenes
MEDiss_control = 1 - cor(MEs_control)
# Cluster module eigengenes
METree_control = hclust(as.dist(MEDiss_control), method = "average");
# Plot the result of module eigengenes
sizeGrWindow(8,16)
pdf("Clustering_of_module_eigengenes_bicor_control.pdf",height=8,width=16)
plot(METree_control, main = "Clustering of module eigengenes bicor control",
     xlab = "", sub = "")
## We choose a height cut of 0.2, corresponding to correlation of 0.80, to merge
MEDissThres = 0.3
# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")
dev.off()
# Call an automatic merging function
merge_control = mergeCloseModules(datExpr_control, 
                                  dynamicColors_control, 
                                  cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors_control = merge_control$colors
# Eigengenes of the new merged modules:
mergedMEs_control = merge_control$newMEs;
# Rename to moduleColors
moduleColors_control = mergedColors_control
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors())
moduleLabels_control = match(moduleColors_control, colorOrder) -1;
MEs_control = mergedMEs_control
length(table(moduleLabels_control))
# Save module colors and labels for use in subsequent parts
save(MEs_control, moduleLabels_control, 
     moduleColors_control, geneTree_control, file = "module_colorsNlabels_control.RData")
load("module_colorsNlabels_control.RData")
print("Step5 - mergeing finished")
#=================================================================================================
#                              5. plotting heatmap                                            ###
#=================================================================================================
#pdf("Network_heatmap_plot_all_genes.pdf",height=8,width=16)
#i = 8
#for (i in c(6:12)){
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
#plotTOM_b = dissTOM14_b_cl^i
# Set diagonal to NA for a nicer plot
#diag(plotTOM_b) = NA
# Call the plot function
#TOMplot(plotTOM_b, geneTree14_b_cl, moduleColors14_b_cl, 
#        main = "Network heatmap plot bicor - all genes")
#}
#dev.off()
#print("Step6 - heapmap created")
#===============================================================================================
#                           7. plot cross-condition dendrogram                               ###
#===============================================================================================
pdf("Gene_dendrogram_cross_condition_bicor.pdf",height=8,width=16)
plotDendroAndColors(geneTree_control, moduleLabels_control, "Modules", dendroLabels=F, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors bicor (control)")
plotDendroAndColors(geneTree_treatment, moduleLabels_control, "Modules", dendroLabels=F,hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors bicor (treatment)")
dev.off()
print("Step7 - cross condition dendrogram created")
#=================================================================================================
#         8. Qualitatively and quantitatively measure network preservation at the module level  ##
#=================================================================================================
######### To quantify this result module preservation statistics ######################
# data pre and check data structure
setLabels = c("control", "treatment")
multiExpr = list(control=list(data = adjacency_control),
                 treatment=list(data = adjacency_treatment))
multiColor = list(control= moduleColors_control)
names(multiExpr) = setLabels
# permutation
mp = modulePreservation(multiExpr,
                        multiColor,
                        referenceNetworks=1,
                        verbose=3,
                        corFnc = "bicor",
                        networkType="unsigned", nPermutations = 2000,
                        maxGoldModuleSize = 500, maxModuleSize = 500,
                        calculateQvalue = T,
                        calculateCor.kIMall = T,
                        calculateClusterCoeff = T,
                        indent = 3)
stats = mp$preservation$Z$ref.control$inColumnsAlsoPresentIn.treatment
Results_mp = stats[order(-stats[,2]),c(1:2)]
save(mp, file = "modulePreservation_bicor_methionine.RData")
load("modulePreservation_bicor_methionine.RData")
print("Step8 - mp finished and data saved")
################ output - shortest - only p summmary  ######################
write.csv(Results_mp,"module_size_and_preservation_statistics_bicor.csv")
################ output - shortest - only p summmary  ######################
# specify the reference and the test networks
ref=1; test = 2
### print results - short version
statsObs = cbind(mp$quality$observed[[ref]][[test]][,-1], mp$preservation$observed[[ref]][[test]][,-1])
statsZ =   cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][,-1])
# Compare preservation to quality:
print(cbind(statsObs[,c("medianRank.pres", "medianRank.qual")],
            signif(statsZ[,c("Zsummary.pres", "Zsummary.qual")],2)))
### print results - full
Obs.PreservationStats= mp$preservation$observed[[ref]][[test]]
Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
# Look at the observed preservation statistics
# Obs.PreservationStats
# Z statistics from the permutation test analysis Z.PreservationStats
write.csv(Obs.PreservationStats,"Obs_PreservationStats_bicor.csv")
write.csv(Z.PreservationStats,"Z_PreservationStats_bicor.csv")
print("Step8 - preservation statistics calculated and saved")

#===========================================================================================
#                                9. mp visualization                                      ##
#===========================================================================================
####### ###### ######   present Z summary ###### ###### ###### ###### ###### 
# Let us now visualize the data.
modColors = rownames(Obs.PreservationStats)
moduleSize = Obs.PreservationStats$moduleSize
# we will omit the grey module (background genes) and the gold module (random sample of genes)
selectModules = !(modColors %in% c("grey", "gold"))
# Text labels for points
point.label = modColors[selectModules]
#Composite preservation statistics
medianRank=Obs.PreservationStats$medianRank.pres
Zsummary=Z.PreservationStats$Zsummary.pres
#
pdf("medianRank_Zsummary_versus_module_size_bicor.pdf",height = 8, width = 16)
par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))
# plot medianRank versus module size
plot(moduleSize[selectModules],medianRank[selectModules],col=1,
     bg=modColors[selectModules], pch = 21,main="medianRank Preservation bicor",
     cex = 2, ylab ="medianRank",xlab="Module size", log="x")
labelPoints(moduleSize[selectModules],medianRank[selectModules],point.label,cex=1,offs=0.03)
# plot Zsummary versus module size
plot(moduleSize[selectModules],Zsummary[selectModules], col = 1,
     bg=modColors[selectModules],pch = 21,main="Zsummary Preservation bicor",
     cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
labelPoints(moduleSize[selectModules],Zsummary[selectModules],point.label,cex=1,offs=0.03)
# Add threshold lines for Zsummary
abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)
dev.off()
print("Step9 - preservation statistics calculated vis and saved")
############################################################################################################################
####### ###### ######   present individual Z  ###### ###### ###### ###### ###### 
# Re-initialize module color labels and sizes
ref = 1;test = 2
# Module labels and module sizes are also contained in the results
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][,1];
# Exclude improper modules / leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs= match(modColors[plotMods], standardColors());

# Compare preservation to quality:
print(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
            signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)))

# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData_b = cbind(mp$preservation$observed[[ref]][[test]][,2],
                   mp$preservation$Z[[ref]][[test]][,2])
# Main titles for the plot
mains = c("Preservation Median rank bicor", "Preservation Zsummary bicor")
# Start the plot: open a suitably sized graphical window and set sectioning and margins.
# Plot each Z statistic in a separate plot.
pdf("all_module_preservation_statistics_nested_bicor.pdf",height = 8, width = 16)
par(mfrow = c(4,5))
for (s in 1:ncol(statsZ)){
  min = min(statsZ[plotMods, s], na.rm = TRUE)
  max = max(statsZ[plotMods, s], na.rm = TRUE)
  if (min > -max/5) min = -max/5
  plot( moduleSizes[plotMods], statsZ[plotMods,s], col = 1, bg = modColors[plotMods], pch = 21,
        main = colnames(statsZ)[s],
        cex = 1.7,
        ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
        ylim = c( min - 0.1 * (max-min), max + 0.1 * (max-min) ),
        xlim = c(20, 1000))
  labelPoints(moduleSizes[plotMods], statsZ[plotMods,s],labs,cex = 0.7, offs = 0.04)
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}
dev.off()
print("Step9 - all_module_preservation_statistics finished and data saved")
