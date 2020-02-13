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
# setwd('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/Net/')
library(readxl)
library(tidyverse)
# DiffC2Gene_raw = read_xlsx('DiffC_Gene.xlsx')
# DiffC2Gene.extend = DiffC2Gene_raw %>% 
#   dplyr::filter(Gene != '-')

# genome pre
library(biomaRt)
genome <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",  dataset = "btaurus_gene_ensembl", host = "grch37.ensembl.org")
# gene = getBM(c("ensembl_gene_id","external_gene_name","description", "start_position", "end_position", "chromosome_name"), mart = genome)
# gene_pos_info_bta = dplyr::select(gene,ensembl_gene_id,start_position,end_position,chromosome_name) %>% 
#   arrange(ensembl_gene_id) %>% 
#   dplyr::mutate_at(vars(chromosome_name),add)

# function pre
AGCTcount = function(ENS,
                     genome = genome,
                     type = "ensembl_gene_id", 
                     seqType = "gene_exon_intron",
                     upstream = 5000,
                     downstream = 5000,
                     find = 'CG'){
  library(tidyverse)
  library(biomaRt)
  add = function(x, sep = ''){paste("chr",x,sep = sep)}
  # generates 5' to 3' sequences of the requested type on the correct strand
  seq1 = try(getSequence(id = ENS, 
                         type = type, 
                         seqType = seqType,
                         upstream = upstream,
                         mart = genome),TRUE)
  seq2 = try(getSequence(id = ENS, 
                         type = type, 
                         seqType = seqType,
                         downstream = downstream,
                         mart = genome),TRUE)
  
  if(isTRUE(class(seq1)=="try-error"|class(seq2)=="try-error")) {
    return('Find Error')
    next 
  } else { 
    seq1 = getSequence(id = ENS, 
                       type = type, 
                       seqType = seqType,
                       upstream = upstream,
                       mart = genome)
    seq_p1 = c(seq1[1]);attributes(seq_p1) = NULL
    seq2 = getSequence(id = ENS, 
                       type = type, 
                       seqType = seqType,
                       downstream = downstream,
                       mart = genome)
    seq_p2 = unlist(seq2[1]);attributes(seq_p2) = NULL
    seq_all = paste(seq_p1,
                    substring(seq_p2,nchar(seq_p2) - downstream + 1,nchar(seq_p2)),
                    sep = '')
    # nchar(seq_all)
    Find_loc = data.frame(str_locate_all(seq_all,find))
    total_c = nrow(Find_loc)
    return(total_c)
  }
}

Gene_all = unique(rownames(networkData_normalized))
Gene_net = unique(rownames(networkData_50var_nocrt))
target = 'ENSBTAG00000005635'
Gene_all = Gene_all[-grep(target,Gene_all)]

#Gene_all = Gene_all_total[1:5000]
Genes_C_count_all = data.frame(Gene = c(),
                               total = c())
for (i in seq_along(Gene_all)){
  message('Working on ',Gene_all[i])
  Genes_C_count_all[i,1] = Gene_all[i]
  Count = AGCTcount(Gene_all[i],genome = genome)
  Genes_C_count_all[i,2] = Count
  print(Count)
}

save(Genes_C_count_all,file = 'Genes_C_count_all_2.RData')

# AGCTcount('ENSBTAG00000005925',genome = genome)
# 
# 
# inputs = list(1, 2, 4,'oops',1,2,4)
# for(input in inputs) {
#   x = try(log(input),TRUE)
#   if(isTRUE(class(x)=="try-error")) {
#     print('error')
#     next 
#   } else { 
#     x = log(input)
#     print(x)
#   }
#  }

# for (i in c(1:5)) {
#   
#   skip_to_next <- FALSE
#   test = AGCTcount(Gene_all[i],genome = genome)
#   # Note that print(b) fails since b doesn't exist
#   tryCatch(print(), error = function(e) { skip_to_next <<- TRUE})
# 
#   if(skip_to_next) { next }
# }
# 
# 
# which(Gene_all_total == 'ENSBTAG00000001455')
# 
# AGCTcount(Gene_all[859],genome = genome)
# 
# save(Genes_C_count_all,file = 'Genes_C_count_all_5k.RData')
# 
# head(gene_pos_info_bta)
# #load('Genes_C_count_all_5k.RData')
# #Genes_C_count_5k = Genes_C_count_all
# library(httr)
# library(jsonlite)
# library(xml2)
# 
# species <- "cow"
# chrom <- "17"
# start <- 67669019 - 5000     
# end <- 67758720 + 5000
# region <- paste0(chrom, ":", start, "..", end)
# server <- "http://rest.ensembl.org/sequence/region"
# r <- GET(paste(server, species, region, sep = "/"), content_type("text/plain"))
# #print(content(r))
# test = (content(r))
# 
# Find_loc = data.frame(str_locate_all(test,'CG'))
# Find_loc2 = data.frame(str_locate_all(seq_all,'CG'))
# 