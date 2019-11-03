# Founction preparation
source("Function_Source.R")
load("module_colorsNlabels_control.RData")
load("modulePreservation_bicor_methionine.RData")
load("data_expr_allprepare with corrections_top50.RData")
load("Ensembl2Entrez_Convert.RData")
load("Enrich_Ensentials.RData")
#### Read in database
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Reactome_db/")
# lowest_path
NCBI2Reactome_lowest_path = read.csv("NCBI2Reactome.txt",sep = "\t",header = F)
NCBI2Reactome_lowest_path_bt = dplyr::filter(NCBI2Reactome_lowest_path, V6 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V4,V5,V6) %>% 
  dplyr::rename(EntrezID = V1,ReactomeID = V2,Reactome_Description = V4, Source = V5,Species = V6)
#head(NCBI2Reactome_lowest_path_bt,10)
# all_path
NCBI2Reactome_all_path = read.csv("NCBI2Reactome_All_Levels.txt",sep = "\t",header = F)
NCBI2Reactome_all_path_bt = 
  dplyr::filter(NCBI2Reactome_all_path,V6 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V4,V5,V6) %>% 
  dplyr::rename(EntrezID = V1,
                ReactomeID = V2,
                Reactome_Description = V4, 
                Source = V5, 
                Species = V6)
#head(NCBI2Reactome_all_path_bt)
# all_react
NCBI2Reactome_all_react = read.csv("NCBI2Reactome_PE_Reactions.txt",sep = "\t",header = F)
NCBI2Reactome_all_react_bt = 
  dplyr::filter(NCBI2Reactome_all_react,V8 == "Bos taurus") %>% 
  dplyr::select(V1,V4,V6,V2,V3,V7,V8) %>% 
  dplyr::rename(EntrezID = V1,ReactomeID = V4, 
                Reactome_Description = V6,
                ProteinID = V2,
                Protein_Description = V3,
                Source = V7, Species = V8)
#head(NCBI2Reactome_all_react_bt,50)

# turn data input as charactor
NCBI2Reactome_all_react_bt[] <-   lapply(NCBI2Reactome_all_react_bt, function(x) if(is.factor(x)) as.character(x) else x)
NCBI2Reactome_lowest_path_bt[] <- lapply(NCBI2Reactome_lowest_path_bt, function(x) if(is.factor(x)) as.character(x) else x)
NCBI2Reactome_all_path_bt[] <-   lapply(NCBI2Reactome_all_path_bt, function(x) if(is.factor(x)) as.character(x) else x)
setwd("/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/")
### Container pre
Total_list_out_entrez = Sig_list_out_entrez
Sig_list_out_entrez = Total_list_out_entrez
TestingSubsetNames = nonpres_modulenames_b
Sig_list_out = Sig_list_out_ens
# Run loops
#===========================================================================================
#                             12. Reactome  enrichment                                    ##
#===========================================================================================
## all react
Reactome_Enrich_all_react_1102 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                 sig_genes_all=Sig_list_out_entrez,
                                                 TestingSubsetNames = TestingSubsetNames,
                                                 InputSource=  NCBI2Reactome_all_react_bt,
                                                 Sig_list_out = Sig_list_out,
                                                 Reacthres = 0.05,
                                                 keyword = "Reactome_Enrichment_all_react_1102")
## lowest path
Reactome_Enrich_lowest_path_1102 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                   sig_genes_all=Sig_list_out_entrez,
                                                   TestingSubsetNames = TestingSubsetNames,
                                                   InputSource=  NCBI2Reactome_lowest_path_bt,
                                                   Sig_list_out = Sig_list_out,
                                                   Reacthres = 0.05,
                                                   keyword = "Reactome_Enrich_lowest_path_1102")
## all path
Reactome_Enrich_all_path_1102 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                sig_genes_all=Sig_list_out_entrez,
                                                TestingSubsetNames = TestingSubsetNames,
                                                InputSource=  NCBI2Reactome_all_path_bt,
                                                Sig_list_out = Sig_list_out,
                                                Reacthres = 0.05,
                                                keyword = "Reactome_Enrich_all_path_1102")


