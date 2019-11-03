# Founction preparation
source("Function_Source.R")
load("module_colorsNlabels_control.RData")
load("modulePreservation_bicor_methionine.RData")
load("data_expr_allprepare with corrections_top50.RData")
load("Enrich_Ensentials.RData")
# Container pre
Total_list_out_entrez = Sig_list_out_entrez
Sig_list_out_entrez = Total_list_out_entrez
TestingSubsetNames = nonpres_modulenames_b
# Run loops
load("Ensembl2Entrez_Convert.RData")

# Run loops
#===========================================================================================
#                             14. Interpro enrichment                                    ##
#===========================================================================================
Interpro_Enrich_Results_thres005_1102 = 
  InterPro_Enrich(total_genes_all = Total_list_out_ens,
                  sig_genes_all = Sig_list_out_ens,
                  TestingSubsetNames = TestingSubsetNames,
                  IPthres = 0.05,
                  biomart="ensembl",
                  dataset="btaurus_gene_ensembl",
                  Identifier = "ensembl_gene_id",
                  attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
                  keyword = "Interpro_Enrichment_1102")

