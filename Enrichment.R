# Founction preparation
source("Function_Source.R")
load("module_colorsNlabels_control.RData")
load("modulePreservation_bicor_methionine.RData")
load("data_expr_allprepare with corrections_top50.RData")
# preservation stats pre
ref=1; test = 2
Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
Zsummary=Z.PreservationStats$Zsummary.pres
# label pre
nonpres_index_b = (which(Zsummary < 2))
nonpres_modulenames_b = rownames(Z.PreservationStats)[nonpres_index_b]
nonpres_modulenames_b = nonpres_modulenames_b[-grep("gold",nonpres_modulenames_b)]
# convert ensembl to entrez(ncbi)
TestingModAssign = moduleColors_control
table(TestingModAssign)
bg_gene = rownames(networkData_final)
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
