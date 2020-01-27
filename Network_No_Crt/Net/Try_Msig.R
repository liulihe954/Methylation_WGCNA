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

###

# specify database location; Read in database
library(msigdbr)
Msig_db_destination = '/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Msig_db'
setwd(Msig_db_destination)
load('Msigdb_bta.RData')
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net")


test = Msig_Enrich(m_df_all = m_df_all,
                   total_genes_all=Total_list_out_entrez,
                   sig_genes_all=Sig_list_out_entrez,
                   TestingSubsetNames = TestingSubsetNames,
                   Sigthres = 0.05,
                   Sig_list_out,
                   DB_List = DB_List,
                   keyword = "Msig_Enrichment_0124")

############################################################
### =======               Msig                 ========== ##
############################################################
load("Msig_Enrichment_0124.RData")

# get loop index
all_module = character()
for (i in seq_along(names(Results_b))){
  all_module[i] = unlist(strsplit(names(Results_b)[i]," "))[1]
}

#names(Interpro_results_b[[1]])
# loop
all_Msig_results = list()
for (i in seq_along(all_module)){
  tmp_name = all_module[i]
  tmp_results = Parse_Results(Results_b[i], keyword= "-")
  if (!(dim(tmp_results)[1] == 0)){
    tmp_results = dplyr::select(tmp_results,-ExternalLoss_total,-InternalLoss_sig) %>% dplyr::arrange(pvalue_r) 
  }
  all_Msig_results[[i]] = tmp_results
  names(all_Msig_results)[i] = all_module[i]
}

require(openxlsx)
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/enrich_results")
write.xlsx(all_interpro_results,file = "Msig_Results_all_0124.xlsx")


