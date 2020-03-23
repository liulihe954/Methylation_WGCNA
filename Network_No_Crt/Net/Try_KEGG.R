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
Sig_list_out = Sig_list_out_ens
load("Ensembl2Entrez_Convert.RData")
# Run loops
Kegg_Enrichment_pval005_1102 = Kegg_Enrich_Plot(sig_genes_all = Sig_list_out_entrez,
                                                total_genes_all = Total_list_out_entrez,
                                                TestingSubsetNames = TestingSubsetNames,
                                                KEGGthres = 0.05, 
                                                species = "bta", 
                                                id.type = "kegg",
                                                Sig_list_out =Sig_list_out,
                                                keyword = "Kegg_Enrichment_0113")
############################################################
### =======                KEGG                ========== ##
############################################################
load("Kegg_Enrichment_0113.RData")

#KEGG_results_b = KEGG_results_b_raw

# get loop index
all_module = character()
for (i in seq_along(names(KEGG_results_b))){
  all_module[i] = unlist(strsplit(names(KEGG_results_b)[i]," "))[1]
}
# loop
all_kegg_results = list()
for (i in seq_along(all_module)){
  tmp_name = all_module[i]
  tmp_results = Parse_Results(KEGG_results_b[i], keyword= "-")
  if (!(dim(tmp_results)[1] == 0)){
    tmp_results = dplyr::select(tmp_results,-ExternalLoss_total,-InternalLoss_sig) %>% dplyr::arrange(pvalue_r) 
  }
  all_kegg_results[[i]] = tmp_results
  names(all_kegg_results)[i] = all_module[i]
}


require(openxlsx)
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/enrich_results")
write.xlsx(all_kegg_results,file = "KEGG_Results_all_0323.xlsx")



# test = load("Msig_Enrichment_0124.RData")
# Msig_Results_b_raw = Results_b_raw
# list = ls()[grep('raw', ls())];list = list[1:5]
# 
# 
# 
# i = 5
# for (i in seq_along(list)){
#   file_tmp = paste(list[i],'with_qvalue.xlsx',sep = '')
#   Add_Q(list[i],file = file_tmp)
# }
# 
# Add_Q = function(results_all,
#                  qthres = 0.05,
#                  file = "test.xlsx"){
#   library(tidyverse)
#   #
#   # get loop index
#   all_module = character()
#   for (i in seq_along(names(results_all))){
#     all_module[i] = unlist(strsplit(names(results_all)[i]," "))[1]
#   }
#   # loop
#   all_results_out = list()
#   for (i in seq_along(all_module)){
#     tmp_name = all_module[i]
#     tmp_results = Parse_Results(results_all[i], keyword= "-")
#     if (!(dim(tmp_results)[1] == 0)){
#       tmp_results = dplyr::select(tmp_results,-ExternalLoss_total,-InternalLoss_sig)
#     }
#     all_results_out[[i]] = tmp_results
#     names(all_results_out)[i] = all_module[i]
#   }
#   #
#   all_results_out_q = list()
#   for (i in seq_along(names(all_results_out))){
#     tmp = all_results_out[[i]] %>%
#       mutate(qvalue = p.adjust(pvalue_r)) %>%
#       dplyr::filter(qvalue <= qthres) %>%
#       arrange(qvalue)
#     all_results_out_q[[i]] = tmp
#   }
#   #return(all_results_out_q)
#   require(openxlsx)
#   setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/enrich_results")
#   write.xlsx(all_go_results_q,file = file)
# }
# 
# 
