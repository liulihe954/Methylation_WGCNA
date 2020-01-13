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

# preservation stats pre
ref=1; test = 2
Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
Zsummary=Z.PreservationStats$Zsummary.pres
# label pre
nonpres_index_b = (which(Zsummary < 2))
nonpres_modulenames_b = rownames(Z.PreservationStats)[nonpres_index_b]
nonpres_modulenames_b = nonpres_modulenames_b[-grep("grey",nonpres_modulenames_b)]
# convert ensembl to entrez(ncbi)
TestingModAssign = moduleColors_control
bg_gene = rownames(network_final)
TestingSubsetNames = nonpres_modulenames_b 

Convert = ConvertNformat(bg_gene,
                         TestingSubsetNames,
                         TestingModAssign,
                         keyword = "Ensembl2Entrez_Convert")
load("Ensembl2Entrez_Convert.RData")

save(Sig_list_out_entrez,
     Total_list_out_entrez,
     nonpres_modulenames_b,
     Sig_list_out_ens,
     file = "Enrich_Ensentials.RData")
