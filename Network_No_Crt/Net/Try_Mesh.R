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

# MeSH db pre
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Mesh_db")
keyword_outer = "MeshDB2"
DB = paste(keyword_outer,".RData",sep = "")
load(DB)
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net")


test_all = unlist(Total_list_out_entrez[1])
test_sig = unlist(Sig_list_out_entrez[6]);attributes(test_sig) = NULL
######
meshParams <- new("MeSHHyperGParams", geneIds = test_sig,
                  universeGeneIds = test_all,
                  annotation = "MeSH.Bta.eg.db", category = "G", database = "gene2pubmed",
                  pvalueCutoff = 0.05, pAdjust = "none")
meshR <- meshHyperGTest(meshParams)
out = data.frame(meshR@ORA$MESHID, meshR@ORA$MESHTERM,
                 meshR@ORA$Size, meshR@ORA$Count, signif(meshR@ORA$Pvalue,2))
colnames(out) = c("MeSH Term ID", "MeSH Term Name",
                  "NT.Genes", "DE Genes", "P-value")
print(unique(out), row.names = F)
all = summary(meshR)

# Run loops
#==============================================================================================
#                                      11. Mesh enrichment                                   ##
#==============================================================================================
#TestingSubsetNames
MESH_Enrichment_1102 = MESH_Enrich(total_genes_all = Total_list_out_entrez,
                                   sig_genes_all = Sig_list_out_entrez,
                                   TestingSubsetNames = TestingSubsetNames,
                                   Meshthres = 0.05,
                                   Sig_list_out = Sig_list_out,
                                   MeshCate = c("D","G"),
                                   dataset="MeSH.Bta.eg.db",
                                   keyword = "MESH_Enrichment_0115")

############################################################
### =======                Mesh                ========== ##
############################################################
load("MESH_Enrichment_0115.RData")
# get loop index
#Mesh_results_b = Mesh_results_b_raw
all_module = character()
for (i in seq_along(names(Mesh_results_b))){
  all_module[i] = unlist(strsplit(names(Mesh_results_b)[i]," "))[1]
}
# loop
all_mesh_results = list()
for (i in seq_along(all_module)){
  tmp_name = all_module[i]
  tmp_results = Parse_Results(Mesh_results_b[i], keyword= "-")
  if (!(dim(tmp_results)[1] == 0)){
    tmp_results = dplyr::select(tmp_results,-ExternalLoss_total,-InternalLoss_sig) %>% dplyr::arrange(pvalue_r) 
  }
  all_mesh_results[[i]] = tmp_results
  names(all_mesh_results)[i] = all_module[i]
}


require(openxlsx)
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/enrich_results")
write.xlsx(all_mesh_results,file = "Mesh_Results_all_0323.xlsx")
