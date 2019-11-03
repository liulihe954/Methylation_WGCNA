# Founction preparation
source("Functions_Source.R")
load("module_colorsNlabels_control.RData")
load("modulePreservation_bicor_methionine.RData")
load("data_expr_allprepare with corrections_top50.RData")
load("Ensembl2Entrez_Convert.RData")
# Container pre
Total_list_out_entrez = Sig_list_out_entrez
Sig_list_out_entrez = Total_list_out_entrez
TestingSubsetNames = nonpres_modulenames_b
Sig_list_out = Sig_list_out_ens
# MeSH db pre
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Mesh_db/")
keyword_outer = "MeshDB"
DB = paste(keyword_outer,".RData",sep = "")
load(DB)
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/")
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
                                   keyword = "MESH_Enrichment_1102")

