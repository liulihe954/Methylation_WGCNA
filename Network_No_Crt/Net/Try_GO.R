# Founction preparation
setwd('/ufrc/penagaricano/lihe.liu/Methylation_WGCNA')
source("Function_Source.R")
setwd('/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net')
load("data_expr_all_with0prepare_no_corrections_top50.RData")
load("permutedStats-actualModules.RData")
load("modulePreservation_methionine.RData")
load("SoftThres_control.RData")
load("modulePreservation_methionine.RData")
load("module_colorsNlabels_control.RData")

load("Enrich_Ensentials.RData")
# Container pre
Total_list_out_entrez = Sig_list_out_entrez
Sig_list_out_entrez = Total_list_out_entrez
TestingSubsetNames = nonpres_modulenames_b

load("Ensembl2Entrez_Convert.RData")
# Run loops
#===========================================================================================
#                             13. Gene Ontology enrichment                                ##
#===========================================================================================
Enrich_Results_thres005_1102 = Go_Enrich_Plot(total_genes_all = Total_list_out_ens,
                                              sig_genes_all = Sig_list_out_ens,
                                              TestingSubsetNames = TestingSubsetNames,
                                              GOthres = 0.05,
                                              keyword = "GO_Enrichment_0113")

##################################################
### =======             GO           ========== ##
##################################################te
load("GO_Enrichment_0113.RData")

#GO_results_b = GO_results_b_raw
# get loop index
all_module = character()
for (i in seq_along(names(GO_results_b))){
  all_module[i] = unlist(strsplit(names(GO_results_b)[i]," "))[1]
}
# get names
biomart="ensembl";dataset="btaurus_gene_ensembl";attributes = c("go_id","namespace_1003")
database = useMart(biomart);genome = useDataset(dataset, mart = database);gene = getBM(attributes,mart = genome)
namespace_index = dplyr::filter(gene,go_id != "",namespace_1003 != "")

# loop
all_go_results = list()
for (i in seq_along(all_module)){
  tmp_name = all_module[i]
  tmp_go_results = Parse_Results(GO_results_b[i], keyword= "-")
  tmp_go_results = dplyr::select(tmp_go_results,-ExternalLoss_total,-InternalLoss_sig) %>% 
    dplyr::arrange(pvalue_r) %>% 
    dplyr::left_join(namespace_index, by =c("GOID" ="go_id")) %>% 
    dplyr::rename(go_id = GOID)
  all_go_results[[i]] = tmp_go_results
}
names(all_go_results) = all_module


require(openxlsx)
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/enrich_results")
write.xlsx(all_go_results,file = "GO_Results_all_0113.xlsx")


tmp = read.xlsx('Supplementary_File.xlsx',sheet = 1)
write.xlsx(tmp,file = 'Supplementary_File1_MappingStats.xlsx')
write.xlsx(tmp,file = 'Supplementary_File2_ModulePreservationStats.xlsx')
write.xlsx(tmp,file = 'Supplementary_File3_GeneMeasurements.xlsx')
write.xlsx(tmp,file = 'Supplementary_File4_GeneCpGCountByRegion.xlsx')
write.xlsx(tmp,file = 'Supplementary_File5_DiffMethCpGAssociations.xlsx')

# 
# library(openxlsx)
# library(tidyverse)
# GO_OUT = read.xlsx('GO_Results_all_0113.xlsx') %>% 
#   dplyr::filter(pvalue<= 0.001)
# 
# ReduceDim_GO_Plot(GO_OUT,
#                   GOthres = 0.001,
#                   label_sizeCC = 0.4,
#                   label_sizeBP = 0.4,
#                   label_sizeMF = 0.4,
#                   Database = "org.Bt.eg.db",
#                   measure="Jiang",combine=NULL,
#                   Dataset_Name = "Methylation_GO_Enrich")
# 
# Enrich_Out = GO_OUT
# GOthres = 0.001
# label_sizeCC = 0.4
# label_sizeBP = 0.4
# label_sizeMF = 0.4
# Database = "org.Bt.eg.db"
# measure="Jiang"
# combine=NULL
# 
#   # load libraries + download ref database
#   library(GOSemSim);library(corrplot);library(tidyverse)
#   do.call(library,list(Database))
#   semData_BP <- godata(paste(Database), ont="BP", computeIC=T)
#   semData_MF <- godata(paste(Database), ont="MF", computeIC=T)
#   semData_CC <- godata(paste(Database), ont="CC", computeIC=T)
#   # selection + formating: for each category we have one vector containing all the sig GO terms
#   BP_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "biological_process") %>% 
#     dplyr::select(go_id) %>% unlist();attributes(BP_List) = NULL # name is an attribute and we dont them, so set null
#   CC_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "cellular_component") %>% 
#     dplyr::select(go_id) %>% unlist();attributes(CC_List) = NULL
#   MF_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "molecular_function") %>% 
#     dplyr::select(go_id) %>% unlist();attributes(MF_List) = NULL
# 
# 
#   goSimMatrix_BP = GOSemSim::mgoSim(BP_List,
#                                     BP_List,
#                                     semData=semData_BP,measure=measure,combine = combine)
#   suspectID_BP = rownames(goSimMatrix_BP)[is.na(goSimMatrix_BP[,1])]
#   if (length(suspectID_BP) != 0){BP_List_new = setdiff(BP_List,suspectID_BP)
#   message(length(suspectID_BP)," invalid ID captured in BP: ",suspectID_BP,", thus been removed!")
#   } else {BP_List_new = BP_List;message("Nice! All IDs are valid in BP!")}
#   goSimMatrix_BP_new = GOSemSim::mgoSim(BP_List_new,
#                                         BP_List_new,
#                                         semData=semData_BP,measure=measure,combine = combine)
#   colnames(goSimMatrix_BP_new) = paste(BP_List_new,Enrich_Out$GO_Name[(Enrich_Out$go_id %in% BP_List_new)])
#   rownames(goSimMatrix_BP_new) = paste(Enrich_Out$GO_Name[(Enrich_Out$go_id %in% BP_List_new)],BP_List_new)
#   # For CC
#   goSimMatrix_CC = GOSemSim::mgoSim(CC_List,
#                                     CC_List,
#                                     semData=semData_CC,measure=measure,combine = combine)
#   suspectID_CC = rownames(goSimMatrix_CC)[is.na(goSimMatrix_CC[,1])]
#   if (length(suspectID_CC) != 0){CC_List_new = setdiff(CC_List,suspectID_CC)
#   message(length(suspectID_CC)," invalid ID captured in CC: ",suspectID_CC,", thus been removed!")
#   } else {CC_List_new = CC_List;message("Nice! All IDs are valid in CC!")}
#   goSimMatrix_CC_new = GOSemSim::mgoSim(CC_List_new,
#                                         CC_List_new,
#                                         semData=semData_CC,measure=measure,combine =combine)
#   colnames(goSimMatrix_CC_new) = paste(CC_List_new,Enrich_Out$GO_Name[(Enrich_Out$go_id %in% CC_List_new)])
#   rownames(goSimMatrix_CC_new) = paste(Enrich_Out$GO_Name[(Enrich_Out$go_id %in% CC_List_new)],CC_List_new)
#   # For MF
#   goSimMatrix_MF = GOSemSim::mgoSim(MF_List,
#                                     MF_List,
#                                     semData=semData_MF,measure=measure,combine = combine)
#   suspectID_MF = rownames(goSimMatrix_MF)[is.na(goSimMatrix_MF[,1])]
#   if (length(suspectID_MF) != 0){MF_List_new = setdiff(MF_List,suspectID_MF)
#   message(length(suspectID_MF)," invalid ID captured in MF: ",suspectID_MF,", thus been removed!")
#   } else {MF_List_new = MF_List;message("Nice! All IDs are valid in MF!")}
#   goSimMatrix_MF_new = GOSemSim::mgoSim(MF_List_new,
#                                         MF_List_new,
#                                         semData=semData_MF,measure=measure,combine = combine)
#   colnames(goSimMatrix_MF_new) = paste(MF_List_new,Enrich_Out$GO_Name[(Enrich_Out$go_id %in% MF_List_new)])
#   rownames(goSimMatrix_MF_new) = paste(Enrich_Out$GO_Name[(Enrich_Out$go_id %in% MF_List_new)],MF_List_new)
#   # Now we take the results and plot
#   pdf(paste("Semantic_Similarity_Measure_",Dataset_Name,"_",formatC(GOthres, format = "e", digits = 0),".pdf",sep = ""))
#   
#   
#   library(ggcorrplot)
#   library(reshape2)
#   
#   goSimMatrix_CC_new_melt <- melt(goSimMatrix_CC_new)
#   ggplot(data = goSimMatrix_CC_new_melt, aes(Var2, Var1, fill = value))+
#     geom_tile(color = "white")+
#     scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
#                          midpoint = 0, limit = c(-1,1), space = "Lab", 
#                          name="Pearson\nCorrelation") +
#     theme_minimal()+ 
#     theme(axis.text.x=element_blank(),
#           axis.title.x=element_blank(),
#           axis.title.y=element_blank())+
#     coord_fixed()
#   
#   goSimMatrix_BP_new_melt <- melt(goSimMatrix_BP_new)
#   ggplot(data = goSimMatrix_BP_new_melt, aes(Var2, Var1, fill = value))+
#     geom_tile(color = "white")+
#     scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
#                          midpoint = 0, limit = c(-1,1), space = "Lab", 
#                          name="Pearson\nCorrelation") +
#     theme_minimal()+ 
#     theme(axis.text.x=element_blank(),
#           axis.title.x=element_blank(),
#           axis.title.y=element_blank())+
#     coord_fixed()
#   
# 
#   
#   corrplot(goSimMatrix_CC_new,title = "Semantic_Similarity_Measure_CC",
#            tl.col = "black", tl.cex = label_sizeCC, 
#            method = "shade", order = "hclust", 
#            hclust.method = "centroid", is.corr = FALSE,mar=c(0,0,1,0))
#   
#   corrplot(goSimMatrix_BP_new,title = "Semantic_Similarity_Measure_BP",
#            tl.col = "black", tl.cex = label_sizeBP, 
#            method = "shade", order = "hclust", 
#            hclust.method = "centroid", is.corr = FALSE,mar=c(0,0,1,0))
#   
#   
#   
#   corrplot(goSimMatrix_MF_new,title = "Semantic_Similarity_Measure_MF",
#            tl.col = "black", tl.cex = label_sizeMF, 
#            method = "shade", order = "hclust", 
#            hclust.method = "centroid", is.corr = FALSE,mar=c(0,0,1,0))



