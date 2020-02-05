# dataset pre
options(stringsAsFactors = FALSE)
data_loci = '/ufrc/penagaricano/lihe.liu/Methylation_WGCNA'
setwd(data_loci)
load("methCov08Stat.rda")
library(tidyverse)
library(methylKit)
chr_index = paste(rep('chr',30),c(seq(1,29),'X'),sep = "")
Diff_C_all = getData(methCov08Stat) %>% 
  mutate_at(vars(chr),as.character) %>% 
  dplyr::filter(chr %in% chr_index)


# dim(Diff_C_all)
#find.package('biomaRt')
# using biomart
#BiocManager::install("biomaRt")
#install.packages('biomaRt_2.42.0.tgz',repos = NULL)
library(biomaRt)
#sessionInfo()
# get postiion info
genome <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",  dataset = "btaurus_gene_ensembl", host = "grch37.ensembl.org")
gene = getBM(c("ensembl_gene_id","external_gene_name","description", "start_position", "end_position", "chromosome_name"), mart = genome)
gene_pos_info_bta = dplyr::select(gene,ensembl_gene_id,start_position,end_position,chromosome_name) %>%  arrange(ensembl_gene_id)
gene_pos_info_bta$chromosome_name = paste('chr',gene_pos_info_bta$chromosome_name,sep = "")

# output pre
gene_bta_diff_c_out = dplyr::select(gene_pos_info_bta,ensembl_gene_id,chromosome_name)


# 
thres = 0.05
sig_type = 'qvalue'
lanchor = 5000
ranchor = 5000
ext = 'all'
# ext = 'prmt'

MythEval = function(gene_pos_info_bta,
                    Diff_C_all,
                    qthres = 0.1,
                    sig_type = 'qvalue',
                    fthres = 0.2,
                    diff_cname = 'meth.diff',
                    lanchor = 5000,
                    ranchor = 5000,
                    ext = 'all', # all means left and right; 'prmt' FOR promoter ONLY
                    nchr = 30){ # number of chr)
  # function pre
  # fucntion pre
  Rbisect_r =  function(lst, value){
    low=1
    high=length(lst)
    mid=length(lst)%/%2
    if (lst[low]==value) low
    else if (lst[high]==value) high
    else{
      while (lst[mid] != value) {
        if (value > lst[mid]){
          low = mid+1
        } else if (value < lst[mid]) {
          high = mid - 1
        } 
        if(high<low){
          mid=low-1;break
        }
        mid=(low+high)%/%2
      }
      mid
    }
  }
  #
  "/" <- function(x,y) ifelse(y==0,0,base:::"/"(x,y)) # special division
  # index pre
  chr_index = paste(rep('chr',30),c(seq(1,29),'X'),sep = "")
  out_compile = data.frame(ensembl_gene_id=character(),meth=numeric())
  fold_thres = fthres * max(abs(Diff_C_all[,diff_cname]))
  #set container
  for (m in seq_len(nchr)){
    #m = 1
    C_index = Diff_C_all %>% dplyr::filter(chr == chr_index[m]) %>% dplyr::pull(start)
    #ptm <- proc.time()
    message('Working on chromosome ',m," :")
    if(ext == 'all'){
      tmp = dplyr::filter(gene_pos_info_bta,chromosome_name == chr_index[m]) %>% 
        arrange(start_position) %>% #dplyr::slice(1:30) %>% 
        mutate(start_position_ext = start_position - lanchor) %>% 
        mutate(end_position_ext = start_position + ranchor)
      index_start = tmp$start_position_ext 
      index_end = tmp$end_position_ext
    } else if (ext == 'prmt') {
      tmp = dplyr::filter(gene_pos_info_bta,chromosome_name == chr_index[m]) %>% 
        arrange(start_position) %>% #dplyr::slice(1:30) %>% 
        mutate(start_position_ext = start_position - lanchor) %>% 
        mutate(end_postition_ext = start_position)
      index_start = tmp$start_position_ext
      index_end = tmp$end_postition_ext
    } else {
      tmp = dplyr::filter(gene_pos_info_bta,chromosome_name == chr_index[m]) %>% 
        arrange(start_position)
      index_start = tmp$start_position 
      index_end = tmp$end_position
    }
    #
    out_total = rep(0,dim(tmp)[1])
    total_c_count = rep(0,dim(tmp)[1])
    sig_c_count = rep(0,dim(tmp)[1])
    #
    for (i in seq_along(C_index)){
      #i = 10000
      if (i%%10000 == 0) {message("tryingd on ",i,"th ", "location ",C_index[i])}
      tmp_c = C_index[i]
      out = Rbisect_r(index_start,tmp_c)
      out_total[i] = length(out)
      multi_check = which(index_end[1:out] >= tmp_c)
      #
      con1 = (length(multi_check)> 0)
      con2 = (Diff_C_all[i,sig_type] <= qthres)
      con3 = (Diff_C_all[i,diff_cname]>= fold_thres)
      #
      if(con1){
        total_c_count[multi_check] = total_c_count[multi_check]+1
        if (con2 & con3){sig_c_count[multi_check] = sig_c_count[multi_check]+1}
      }
    }
    tmp_out = tmp %>% 
      dplyr::mutate(meth_total = total_c_count) %>% 
      dplyr::mutate(meth_sig = sig_c_count) %>% 
      dplyr::mutate(meth = sig_c_count/total_c_count) %>% 
      dplyr::select(ensembl_gene_id,meth)
    #proc.time() - ptm
    out_compile = rbind(out_compile,tmp_out)
  }
  return(out_compile) # list(MythStatus = 
}

Meth_all = MythEval(gene_pos_info_bta,
                    Diff_C_all,
                    qthres = 0.1,
                    sig_type = 'qvalue',
                    fthres = 0.2,
                    diff_cname = 'meth.diff',
                    lanchor = 5000,
                    ranchor = 5000,
                    ext = 'all', # all means left and right; 'prmt' FOR promoter ONLY
                    nchr = 30)

Meth_all = as_tibble(Meth_all) %>% arrange(desc(meth)) #%>% print(n = 100)



Meth_prmt = MythEval(gene_pos_info_bta,
                     Diff_C_all,
                     qthres = 0.1,
                     sig_type = 'qvalue',
                     fthres = 0.2,
                     diff_cname = 'meth.diff',
                     lanchor = 5000,
                     ranchor = 5000,
                     ext = 'prmt', # all means left and right; 'prmt' FOR promoter ONLY
                     nchr = 30)
Meth_prmt = as_tibble(Meth_prmt) %>% arrange(desc(meth)) #%>% print(n = 100)



getwd()
save(Meth_all,Meth_prmt,
     file = "MethEval_all.RData" )
#
load('MethEval_all.RData')

head(Meth_prmt)

# gene_bta_diff_c_out = merge(gene_bta_diff_c_out,
#                                 out_compile,all=T,
#                                 by.x = "ensembl_gene_id", 
#                                 by.y = "ensembl_gene_id",
#                                 suffixes = c("",""))


######=========================##########
##       Retrive CpG info from NCBI    ##
######========================##########
head(gene_pos_info_bta)




