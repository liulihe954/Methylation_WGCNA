######=========================##########
##            Data Pre                 ##
######=========================##########
#setwd('/Users/liulihe95/Desktop/Methionine')
# Founction preparation
setwd('/ufrc/penagaricano/lihe.liu/Methylation_WGCNA')
source("Function_Source.R")
setwd('/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net')
#setwd('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/Net')
load("data_expr_all_with0prepare_no_corrections_top50.RData")
load("permutedStats-actualModules.RData")
load("modulePreservation_methionine.RData")
load("SoftThres_control.RData")
load("modulePreservation_methionine.RData")
load("module_colorsNlabels_control.RData")
load("Enrich_Ensentials.RData")
load("Ensembl2Entrez_Convert.RData")
load('network_final.RData')
# load('MethEval_all.RData')

# Gather Info: KME and Meth
# datKME_tmp = signedKME(datExpr_control, MEs_control)
# datKME = datKME_tmp %>% 
#   dplyr::mutate(ensembl_gene_id = rownames(datKME_tmp)) %>% 
#   dplyr::mutate(MdouleAssign = moduleColors_control) %>% 
#   dplyr::left_join(Meth_prmt, by= c("ensembl_gene_id" = "ensembl_gene_id"))
#%>% dplyr::filter(meth != 1)
# head(datKME)
# summary(datKME$meth)
# table(datKME$ensembl_gene_id%in% Meth_prmt$ensembl_gene_id)

######=========================##########
##            Index Pre                ##
######=========================##########
# ref=1; test = 2
# Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
# Zsummary=Z.PreservationStats$Zsummary.pres
# #
# nonpres_index_b = (which(Zsummary < 2))
# nonpres_modulenames_b = rownames(Z.PreservationStats)[nonpres_index_b]
# Mod_Index_NonPre  = nonpres_modulenames_b[-grep("grey",nonpres_modulenames_b)]
# Mod_Index_Pre = rownames(Z.PreservationStats)[-nonpres_index_b]

######=========================##########
##        diff C prop vs M.M           ##
######=========================##########
# 
# pdf('PDF_Results.pdf')
# for (i in seq_along(Mod_Index_NonPre)){
#   text = Mod_Index_NonPre[i]
#   sub = paste('kME',text,sep = '')
#   print(ggplot(datKME,aes(x=get(sub),y=meth)) + 
#           geom_point(colour="grey") +
#           geom_point(data = subset(datKME,MdouleAssign == text), 
#                      aes(x=get(sub),y=meth),colour="red", size=1)+
#           ggtitle(paste("Methylation VS ModuleMembership",text,sep='-')) +
#           xlab("InModule_Con") + ylab("MethC_Prop")+
#           theme(plot.title = element_text(hjust = 0.5)))
# }
# dev.off()

######=========================##########
##        Myth_extent porp            ##
######========================##########
# gene pre
setwd('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/Net/')
library(readxl)
library(tidyverse)
DiffC2Gene_raw = read_xlsx('DiffC_Gene.xlsx')
DiffC2Gene.extend = DiffC2Gene_raw %>%
  dplyr::filter(Gene != '-')


# genome pre
library(biomaRt)
# genome <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",  dataset = "btaurus_gene_ensembl", host = "grch37.ensembl.org")
# gene = getBM(c("ensembl_gene_id","external_gene_name", "start_position", "end_position", "chromosome_name"), mart = genome)
# Gene_genome = length(gene$ensembl_gene_id)
# gene_pos_info_bta = dplyr::select(gene,ensembl_gene_id,start_position,end_position,chromosome_name) %>%
#   arrange(ensembl_gene_id) %>%
#   dplyr::mutate_at(vars(chromosome_name),add)

# # function pre
# AGCTcount = function(ENS,
#                      genome = genome,
#                      type = "ensembl_gene_id", 
#                      seqType = "gene_exon_intron",
#                      upstream = 5000,
#                      downstream = 5000,
#                      find = 'CG'){
#   library(tidyverse)
#   library(biomaRt)
#   add = function(x, sep = ''){paste("chr",x,sep = sep)}
#   # generates 5' to 3' sequences of the requested type on the correct strand
#   seq1 = try(getSequence(id = ENS, 
#                          type = type, 
#                          seqType = seqType,
#                          upstream = upstream,
#                          mart = genome),TRUE)
#   seq2 = try(getSequence(id = ENS, 
#                          type = type, 
#                          seqType = seqType,
#                          downstream = downstream,
#                          mart = genome),TRUE)
#   
#   if(isTRUE(class(seq1)=="try-error"|class(seq2)=="try-error")) {
#     return('Find Error')
#     next 
#   } else { 
#     seq1 = getSequence(id = ENS, 
#                        type = type, 
#                        seqType = seqType,
#                        upstream = upstream,
#                        mart = genome)
#     seq_p1 = c(unlist(seq1[1]));attributes(seq_p1) = NULL
#     seq2 = getSequence(id = ENS, 
#                        type = type, 
#                        seqType = seqType,
#                        downstream = downstream,
#                        mart = genome)
#     seq_p2 = c(unlist(seq2[1]));attributes(seq_p2) = NULL
#     seq_all = paste(seq_p1,
#                     substring(seq_p2,nchar(seq_p2) - downstream + 1,nchar(seq_p2)),
#                     sep = '')
#     # nchar(seq_all)
#     Find_loc = data.frame(str_locate_all(seq_all,find))
#     total_c = nrow(Find_loc)
#     return(total_c)
#   }
# }
AGCTcount_api = function(ENS,
                         upstream = 5500,
                         downstream = 4000,
                         find = 'CG'){
  library(httr);library(jsonlite);library(xml2);library(stringr)
  server = "https://rest.ensembl.org"
  #ext <- "/sequence/id/ENSBTAG00000000070?expand_3prime=5000;expand_5prime=5000"
  ext = paste( "/sequence/id/",ENS,"?expand_3prime=",upstream,";expand_5prime=",downstream,sep = "")
  #ext <- "/sequence/id/ENSBTAG00000000070?"
  r = GET(paste(server, ext, sep = ""), content_type("text/plain"))
  stop_for_status(r)
  seq_all = content(r)
  seq_up = substr(seq_all,1,upstream)
  # nchar(seq_all)
  Find_loc1 = data.frame(str_locate_all(seq_all,find))
  Find_loc2 = data.frame(str_locate_all(seq_up,find))
  total_c1 = nrow(Find_loc1)
  total_c2 = nrow(Find_loc2)
  return(list(Total = total_c1,
              Upstream = total_c2))
  
}

# get all genes _ index
biomart="ensembl";dataset="btaurus_gene_ensembl";Identifier = "external_gene_name";attributes = c("ensembl_gene_id")
database = useMart(biomart)
genome = useDataset(dataset, mart = database)
gene = getBM(attributes,mart = genome)
gene_all_v2 = unique(gene$ensembl_gene_id)


# 
Genes_C_count_all_api = data.frame(Gene = c(),Total = c(),Upstream = c())
for (i in seq_along(gene_all_v2)){
  message('Working on ',gene_all_v2[i])
  Genes_C_count_all_api[i,1] = gene_all_v2[i]
  Count = AGCTcount_api(gene_all_v2[i])
  Genes_C_count_all_api[i,2] = as.character(Count[[1]])
  Genes_C_count_all_api[i,3] = as.character(Count[[2]])
  print(Count)
}

save(Genes_C_count_all_api,
     file = 'Genes_C_count_all_Final_api.RData')

# load("Genes_C_count_all_Final_api.RData")
library(readxl)
test = read_xlsx('DiffC_Gene.xlsx')
test2 = test %>% dplyr::filter(Gene != "-") #%>% dplyr::filter(Region != 'DOWNSTREAM')

Associ_out2 =  test2 %>% 
  #dplyr::select(-Distance,-Transcript,-`Exon/Intron`,-TSSDistance,-PercRegion,-PercArea) %>% 
  group_by(Gene) %>% 
  dplyr::count(Region) %>%
  tidyr::spread(key = Gene, value = n) %>% 
  transpose_df() %>% 
  replace(is.na(.), 0) %>% 
  mutate_at(vars(-Region), as.numeric) %>% 
  dplyr::rename(Gene = Region)

Associ_out_count2 = Associ_out2 %>% 
  mutate(Count1 = `1st_EXON`+ GENE_BODY + INTRON) %>% 
  mutate(Count2 = PROMOTER + TSS + UPSTREAM)

plot(Associ_out_count2$Count1)


load('Total_C_count_raw.rda')
Total_C_count = Total_C_count_raw %>% 
  mutate(Count1_all = `1st_EXON`+ GENE_BODY + INTRON) %>% 
  mutate(Count2_all = PROMOTER + TSS + UPSTREAM)
Total_C_count_2join = Total_C_count %>% 
  dplyr::select(Gene,Count1_all,Count2_all)

"/" <- function(x,y) ifelse(y==0,0,base:::"/"(x,y)) # special division
Associ_out_count_final2 = Associ_out_count2 %>% 
  dplyr::left_join(Total_C_count_2join,by = c('Gene'= 'Gene')) %>% 
  dplyr::select(Gene,Count1,Count2,Count1_all,Count2_all) %>% 
  mutate(Prop_Body = Count1/Count1_all) %>% 
  mutate(Prop_Prpt = Count2/Count2_all) %>% 
  mutate(Prop_All = (Count1+Count2)/(Count1_all+Count2_all))


# Count1 = 1st_EXON+ GENE_BODY + INTRON
# Count2 = PROMOTER + TSS + UPSTREAM
Genes_meth_prop2 = Associ_out_count_final2

Genes_meth_select2 = Genes_meth_prop2 %>% 
  dplyr::filter((Count1 >= 30 | Count2 >= 5)) %>% 
  #dplyr::filter(Prop_All>=quantile(Prop_All,0.5)) %>% 
  #dplyr::filter(Prop_All >= 0.1) %>% 
  arrange(Prop_All)

length(which )

summary(sort(Genes_meth_prop2$Prop_Body,decreasing = T))

summary(sort(Genes_meth_prop$Prop_Body,decreasing = T))

table(Genes_meth_prop$Prop_Prpt)
length(which(Genes_meth_prop$Prop_Prpt >= .1))
summary(Genes_meth_prop2$Prop_All)
Diff_Meth_Gene_index2 = unique(Genes_meth_select2$Gene)
length(Diff_Meth_Gene_index2)


