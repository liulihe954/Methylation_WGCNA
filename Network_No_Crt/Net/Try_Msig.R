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
                   keyword = "Msig_Enrichment")


