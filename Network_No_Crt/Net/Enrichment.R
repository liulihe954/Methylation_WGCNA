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


# local plotting
#setwd('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/enrich_results')
GO_results_antiquewhite2 = read_xlsx('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/enrich_results/GO_Results_all_0113.xlsx')
GO_results_brown2 = read_xlsx('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/enrich_results/GO_Results_all_0113.xlsx',2)
#setwd('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/Net')
test = ReduceDim_GO_Plot(GO_results_brown2,
                         GOthres = 0.01,
                         label_sizeCC = 0.4,
                         label_sizeBP = 0.4,
                         label_sizeMF = 0.4,
                         Database = "org.Bt.eg.db",
                         measure="Jiang",combine=NULL,
                         Dataset_Name = 'GO')
                 
