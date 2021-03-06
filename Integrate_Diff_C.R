# dataset pre
options(stringsAsFactors = FALSE)
data_loci = '/blue/mateescu/lihe.liu/Methylation_WGCNA'
setwd(data_loci)
load("methCov08Stat.rda")
library(tidyverse)
library(methylKit)

# select Significant Diff Cs
chr_index = paste(rep('chr',30),c(seq(1,29),'X'),sep = "")
Diff_C_all = getData(methCov08Stat) %>% 
  mutate_at(vars(chr),as.character) %>% 
  dplyr::filter(chr %in% chr_index) # 5136556

library(readxl)
Diff_C_Sig = Diff_C_all %>% 
  dplyr::filter(qvalue <= 0.10,abs(meth.diff) >= 20)
dim(Diff_C_Sig) # 101094

#save(Diff_C_all,Diff_C_Sig,file = 'Meth_Diff_CpG_Universe.rda')

# make bed file
Diff_C_Sig_BED = Diff_C_Sig  %>% 
  dplyr::select(chr,start,end)
Diff_C_Sig_BED[,1] = str_replace(Diff_C_Sig_BED[,1],'chr','') # 101094 
colnames(Diff_C_Sig_BED) = NULL

# 
# # for alt_splicing # 03022021
# # make bed file
# library(readxl)
# Diff_C_Sig_4alt = Diff_C_all %>% 
#   dplyr::filter(pvalue <= 0.01)
# dim(Diff_C_Sig_4alt) # 
# 
# Diff_C_Sig_BED_4alt = Diff_C_Sig_4alt  %>% 
#   dplyr::select(chr,start,end)
# Diff_C_Sig_BED_4alt[,1] = str_replace(Diff_C_Sig_BED_4alt[,1],'chr','') #
# colnames(Diff_C_Sig_BED_4alt) = NULL
# # rgmatch loci
# #rgmatch_loci = '/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net/rgmatch'
# rgmatch_loci = '/blue/mateescu/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net/rgmatch'
# setwd(rgmatch_loci)
# write.table(Diff_C_Sig_BED_4alt, file='Diff_C_Sig_BED_4alt.bed', row.names = F,quote=FALSE, sep='\t')
# setwd(data_loci)
# ## ================================================================================================================== ##
# #   python rgmatch.py -g Bos_taurus.ARS-UCD1.2.99.gtf -b Diff_C_Sig_BED_4alt.bed -r 'gene' -q 3 -o myassociations_sig_4alt.txt  ##
# ## ================================================================================================================== ##



# TEST_Diff_C = read.xlsx('DiffC_Gene.xlsx')
# table(TEST_Diff_C$Chromosome)
# length(unique(TEST_Diff_C$Position))

# rgmatch loci
rgmatch_loci = '/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net/rgmatch'
setwd(rgmatch_loci)
write.table(Diff_C_Sig_BED, file='1_Sig_Cs_BED.bed', row.names = F,quote=FALSE, sep='\t')
setwd(data_loci)
## ================================================================================================================== ##
#   python rgmatch.py -g Bos_taurus.ARS-UCD1.2.99.gtf -b 1_Sig_Cs_BED.bed -r 'exon' -q 4 -o myassociations_sig.txt  ##
## ================================================================================================================== ##

#Count all Cs
setwd('/blue/mateescu/lihe.liu/Methylation_WGCNA/Network_No_Crt/MethEval')
load('All_Cs.rda')
chr_index = paste(rep('chr',30),c(seq(1,29),'X'),sep = "")
SeqC_all = getData(data) %>% 
  mutate_at(vars(chr),as.character) %>% 
  dplyr::filter(chr %in% chr_index) #
dim(SeqC_all) #  5136556
head(SeqC_all)
#
SeqC_all_BED = SeqC_all %>% 
  dplyr::select(chr,start,end)
SeqC_all_BED[,1] = str_replace(SeqC_all_BED[,1],'chr','')

colnames(SeqC_all_BED) = NULL

# rgmatch loci
rgmatch_loci = '/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net/rgmatch'
setwd(rgmatch_loci)
write.table(SeqC_all_BED, file='2_All_Cs_BED.bed', row.names = F,quote=FALSE, sep='\t')
setwd(data_loci)
## ================================================================================================================== ##
#   python rgmatch.py -g Bos_taurus.ARS-UCD1.2.99.gtf -b 2_All_Cs_BED.bed -r 'exon' -q 4 -o myassociations_all.txt  ##
## ================================================================================================================== ##

# # 
# thres = 0.05
# sig_type = 'qvalue'
# lanchor = 5000
# ranchor = 5000
# ext = 'all'
# # ext = 'prmt'
# 
# MythEval = function(gene_pos_info_bta,
#                     Diff_C_all,
#                     qthres = 0.1,
#                     sig_type = 'qvalue',
#                     fthres = 0.2,
#                     diff_cname = 'meth.diff',
#                     lanchor = 5000,
#                     ranchor = 5000,
#                     ext = 'all', # all means left and right; 'prmt' FOR promoter ONLY
#                     nchr = 30){ # number of chr)
#   # function pre
#   # fucntion pre
#   Rbisect_r =  function(lst, value){
#     low=1
#     high=length(lst)
#     mid=length(lst)%/%2
#     if (lst[low]==value) low
#     else if (lst[high]==value) high
#     else{
#       while (lst[mid] != value) {
#         if (value > lst[mid]){
#           low = mid+1
#         } else if (value < lst[mid]) {
#           high = mid - 1
#         } 
#         if(high<low){
#           mid=low-1;break
#         }
#         mid=(low+high)%/%2
#       }
#       mid
#     }
#   }
#   #
#   "/" <- function(x,y) ifelse(y==0,0,base:::"/"(x,y)) # special division
#   # index pre
#   chr_index = paste(rep('chr',30),c(seq(1,29),'X'),sep = "")
#   out_compile = data.frame(ensembl_gene_id=character(),meth=numeric())
#   fold_thres = fthres * max(abs(Diff_C_all[,diff_cname]))
#   #set container
#   for (m in seq_len(nchr)){
#     #m = 1
#     C_index = Diff_C_all %>% dplyr::filter(chr == chr_index[m]) %>% dplyr::pull(start)
#     #ptm <- proc.time()
#     message('Working on chromosome ',m," :")
#     if(ext == 'all'){
#       tmp = dplyr::filter(gene_pos_info_bta,chromosome_name == chr_index[m]) %>% 
#         arrange(start_position) %>% #dplyr::slice(1:30) %>% 
#         mutate(start_position_ext = start_position - lanchor) %>% 
#         mutate(end_position_ext = start_position + ranchor)
#       index_start = tmp$start_position_ext 
#       index_end = tmp$end_position_ext
#     } else if (ext == 'prmt') {
#       tmp = dplyr::filter(gene_pos_info_bta,chromosome_name == chr_index[m]) %>% 
#         arrange(start_position) %>% #dplyr::slice(1:30) %>% 
#         mutate(start_position_ext = start_position - lanchor) %>% 
#         mutate(end_postition_ext = start_position)
#       index_start = tmp$start_position_ext
#       index_end = tmp$end_postition_ext
#     } else {
#       tmp = dplyr::filter(gene_pos_info_bta,chromosome_name == chr_index[m]) %>% 
#         arrange(start_position)
#       index_start = tmp$start_position 
#       index_end = tmp$end_position
#     }
#     #
#     out_total = rep(0,dim(tmp)[1])
#     total_c_count = rep(0,dim(tmp)[1])
#     sig_c_count = rep(0,dim(tmp)[1])
#     #
#     for (i in seq_along(C_index)){
#       #i = 10000
#       if (i%%10000 == 0) {message("tryingd on ",i,"th ", "location ",C_index[i])}
#       tmp_c = C_index[i]
#       out = Rbisect_r(index_start,tmp_c)
#       out_total[i] = length(out)
#       multi_check = which(index_end[1:out] >= tmp_c)
#       #
#       con1 = (length(multi_check)> 0)
#       con2 = (Diff_C_all[i,sig_type] <= qthres)
#       con3 = (Diff_C_all[i,diff_cname]>= fold_thres)
#       #
#       if(con1){
#         total_c_count[multi_check] = total_c_count[multi_check]+1
#         if (con2 & con3){sig_c_count[multi_check] = sig_c_count[multi_check]+1}
#       }
#     }
#     tmp_out = tmp %>% 
#       dplyr::mutate(meth_total = total_c_count) %>% 
#       dplyr::mutate(meth_sig = sig_c_count) %>% 
#       dplyr::mutate(meth = sig_c_count/total_c_count) 
#     # %>% dplyr::select(ensembl_gene_id,meth)
#     #proc.time() - ptm
#     out_compile = rbind(out_compile,tmp_out)
#   }
#   return(out_compile) # list(MythStatus = 
# }
# 
# 
# Meth_all = MythEval(gene_pos_info_bta,
#                     Diff_C_all,
#                     qthres = 0.1,
#                     sig_type = 'qvalue',
#                     fthres = 0.2,
#                     diff_cname = 'meth.diff',
#                     lanchor = 5000,
#                     ranchor = 5000,
#                     ext = 'all', # all means left and right; 'prmt' FOR promoter ONLY
#                     nchr = 30)
# 
# Meth_all = as_tibble(Meth_all) %>% arrange(desc(meth)) #%>% print(n = 100)
# 
# Meth_prmt = MythEval(gene_pos_info_bta,
#                      Diff_C_all,
#                      qthres = 0.1,
#                      sig_type = 'qvalue',
#                      fthres = 0.2,
#                      diff_cname = 'meth.diff',
#                      lanchor = 5000,
#                      ranchor = 5000,
#                      ext = 'prmt', # all means left and right; 'prmt' FOR promoter ONLY
#                      nchr = 30)
# Meth_prmt = as_tibble(Meth_prmt) %>% arrange(desc(meth)) #%>% print(n = 100)
# 
# 
# getwd()
# save(Meth_all,Meth_prmt,
#      file = "MethEval_all.RData")
# 
# #
# load('MethEval_all.RData')
