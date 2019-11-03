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
# Run loops
#===========================================================================================
#                             13. Gene Ontology enrichment                                ##
#===========================================================================================
Enrich_Results_thres005_1102 = Go_Enrich_Plot(total_genes_all = Total_list_out_ens,
                                              sig_genes_all = Sig_list_out_ens,
                                              TestingSubsetNames = TestingSubsetNames,
                                              GOthres = 0.05,
                                              keyword = "GO_Enrichment_1102")

