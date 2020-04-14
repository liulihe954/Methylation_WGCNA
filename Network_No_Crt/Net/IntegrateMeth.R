######=========================##########
##           0. Data Pre               ##
######=========================##########
Pwd = getwd()
setwd('/Users/liulihe95/Desktop/Methionine')
source("Function_Source.R")
setwd(Pwd)
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
#load('MethEval_all.RData')

######=========================##########
##         1. counted ALL Cs           ##
######=========================##########
# run rgmatch find the Diff C location/ gene assignments
## ================================================================================================================== ##
#     python rgmatch.py -g Bos_taurus.ARS-UCD1.2.99.gtf -b Diff_C_Sig_BED.bed -r 'gene' -q 4 -o myassociations.txt    ##
## ================================================================================================================== ##
#Data_loci = '/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net'
setwd('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/Net')
#setwd('/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net/rgmatch')
# Associ_out_raw = read.table('myassoci_exon_5.5k_ext.txt',sep = '\t') %>% data.frame()
library(readxl)
library(tidyverse)
# custom function to transpose while preserving names
transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.) %>% 
    `colnames<-`(.[1,]) %>%
    .[-1,] %>%
    `rownames<-`(NULL)
  return(t_df)
}

rgmatch_loci = '/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net/rgmatch'
setwd(rgmatch_loci)

Associ_out_all_raw = read.table('myassociations_all.txt',header = TRUE,stringsAsFactors = FALSE)
Associ_out_sig_raw = read.table('myassociations_sig.txt',header = TRUE,stringsAsFactors = FALSE)

#
Associ_out_all_raw = as_tibble(Associ_out_all_raw)
Associ_sig_sig_raw = as_tibble(Associ_out_sig_raw)

outtt = Associ_sig_sig_raw[Associ_sig_sig_raw$Midpoint == 76077,]

head(Associ_out_sig_raw)

length(unique(Associ_sig_sig_raw$Region))

test = Associ_sig_sig_raw %>% 
  group_by(Midpoint) %>% 
  count(Gene)






#
Associ_out_all = Associ_out_all_raw %>%
  #dplyr::filter(!(Area == 'DOWNSTREAM')) %>%
  group_by(Gene) %>%
  dplyr::count(Area) %>%
  tidyr::spread(key = Gene, value = n) %>%
  transpose_df() %>%
  mutate_at(vars(-Area), as.numeric) %>%
  replace(is.na(.), 0) %>%
  dplyr::rename(Gene = Area) %>% # 24344
  mutate(Count_Prpt =  PROMOTER + TSS + UPSTREAM) %>%
  mutate(Count_Body =  `1st_EXON`+ GENE_BODY + INTRON) %>%
  mutate(Count_All = Count_Prpt + Count_Body) %>%
  dplyr::select(Gene,Count_Prpt,Count_Body,Count_All)


Associ_out_sig = Associ_out_sig_raw %>%
  #dplyr::filter(!(Area == 'DOWNSTREAM')) %>%
  group_by(Gene) %>%
  dplyr::count(Area) %>%
  tidyr::spread(key = Gene, value = n) %>%
  transpose_df() %>%
  mutate_at(vars(-Area), as.numeric) %>%
  replace(is.na(.), 0) %>%
  dplyr::rename(Gene = Area) %>% # 10247
  mutate(Count_Prpt =  PROMOTER + TSS + UPSTREAM) %>%
  mutate(Count_Body =  `1st_EXON`+ GENE_BODY + INTRON) %>%
  mutate(Count_All = Count_Prpt + Count_Body) %>%
  dplyr::select(Gene,Count_Prpt,Count_Body,Count_All)
#
#setwd('/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net/rgmatch')

Associ_out_sig_out = Diff_C_Sig %>%
  dplyr::select(chr,start,qvalue,meth.diff) %>%
  left_join(dplyr::select(Associ_out_sig_raw,Midpoint,Gene,Area),
            by = c('start'='Midpoint')) %>%
  replace_na(list(Gene = '-', Area = "INTERGENIC")) %>%
  as_tibble()

x1 = Associ_out_sig_raw %>% dplyr::select(Midpoint) %>% unlist(use.names = F)
x2 = Diff_C_Sig %>% dplyr::select(start) %>% unlist(use.names = F)
head(Associ_out_sig_out[duplicated(Associ_out_sig_out$start),])
Associ_out_sig_out[(which(table(Associ_out_sig_out$start) == 12)),]


Associ_out_all_count = Associ_out_all_raw %>%
  #dplyr::filter(!(Area == 'DOWNSTREAM')) %>%
  group_by(Gene) %>%
  dplyr::count(Area) %>%
  tidyr::spread(key = Gene, value = n) %>%
  transpose_df() %>%
  mutate_at(vars(-Area), as.numeric) %>%
  replace(is.na(.), 0) %>%
  dplyr::rename(Gene = Area) %>%
  as_tibble()


require(openxlsx)
list_of_datasets <-
  list("Associations of Diff Meth CpGs" = Associ_out_sig_out,
       "Diff Meth CpGs Counts by Region" = Associ_out_all_count)
write.xlsx(list_of_datasets, file = "CpG_asso_count.xlsx")

# 
# 
# "/" <- function(x,y) ifelse(y==0,0,base:::"/"(x,y))
# Meth_Prop_Univ = Associ_out_all %>%
#   left_join(Associ_out_sig, by = c('Gene' = 'Gene')) %>%
#   replace_na(list(Count_Prpt.y = 0,Count_Body.y = 0,Count_All.y = 0)) %>%
#   mutate(Prop_Prpt = Count_Prpt.y/Count_Prpt.x,
#          Prop_Body = Count_Body.y/Count_Body.x,
#          Prop_All = Count_All.y/Count_All.x) %>%
#   dplyr::select(Gene,Count_Prpt.y,Count_Body.y,Count_All.y,Prop_Prpt,Prop_Body,Prop_All)
# #
# 
# Meth_Prop_Univ %>%
#   dplyr::filter(Prop_Prpt ==0,
#                 Prop_Body==0,
#                 Prop_All ==0) %>% dim() # 8884
# #
# save(Associ_out_all,Associ_out_sig,Meth_Prop_Univ,file = 'Meth_Prop_Univ.rda')

# #head(Associ_out,100) %>% print(n = Inf)
# Associ_out =  Associ_out_raw %>% 
#   #dplyr::select(-Distance,-Transcript,-`Exon/Intron`,-TSSDistance,-PercRegion,-PercArea) %>% 
#   group_by(Gene) %>% 
#   dplyr::count(Area) %>%
#   tidyr::spread(key = Gene, value = n) %>% 
#   transpose_df() %>% 
#   replace(is.na(.), 0) %>% 
#   mutate_at(vars(-Area), as.numeric) %>% 
#   dplyr::rename(Gene = Area)
# 
# Associ_out_count = Associ_out %>% 
#   mutate(Count1 = `1st_EXON`+ GENE_BODY + INTRON) %>% 
#   mutate(Count2 = PROMOTER + TSS + UPSTREAM) %>% 
#   dplyr::select(Gene,Count1,Count2)
# 
# #
# load('Total_C_count_raw.rda')
# Total_C_count = Total_C_count_raw %>% 
#   mutate(Count1_all = `1st_EXON`+ GENE_BODY + INTRON) %>% 
#   mutate(Count2_all = PROMOTER + TSS + UPSTREAM)
# Total_C_count_2join = Total_C_count %>% 
#   dplyr::select(Gene,Count1_all,Count2_all)
# 
# head(Total_C_count_raw)
# dim(Total_C_count_raw)
# 
# "/" <- function(x,y) ifelse(y==0,0,base:::"/"(x,y)) # special division
# Associ_out_count_final = Associ_out_count %>% 
#   dplyr::left_join(Total_C_count_2join,by = c('Gene'= 'Gene')) %>% 
#   dplyr::select(Gene,Count1,Count2,Count1_all,Count2_all) %>% 
#   mutate(Prop_Body = Count1/Count1_all) %>% 
#   mutate(Prop_Prpt = Count2/Count2_all) %>% 
#   mutate(Prop_All = (Count1+Count2)/(Count1_all+Count2_all))
# 
# 
# # Count1 = 1st_EXON+ GENE_BODY + INTRON
# # Count2 = PROMOTER + TSS + UPSTREAM
# Genes_meth_prop = Associ_out_count_final %>% 
#   dplyr::filter(!is.na(Prop_Body) & Prop_Body <=1) %>% 
#   dplyr::filter(!is.na(Prop_Prpt) & Prop_Prpt <=1) %>% 
#   dplyr::filter(!is.na(Prop_All) & Prop_All <=1)
# 
# Genes_meth_select = Genes_meth_prop %>% 
#   dplyr::filter((Count1 >= 20 | Count2 >= 3)) %>% 
#   dplyr::filter(Prop_Prpt >= quantile(Prop_Prpt,0.4)) %>%
#   dplyr::filter(Prop_All >= quantile(Prop_All,0.4)) %>%
#   dplyr::filter(Prop_Body >= quantile(Prop_Body,0.4)) %>%
#   dplyr::filter(Prop_Prpt != 0) %>%
#   dplyr::filter(Prop_All != 0) %>%
#   dplyr::filter(Prop_Body != 0) %>%
#   arrange(Prop_All)
# 
# Diff_Meth_Gene_index = unique(Genes_meth_select$Gene)
# length(Diff_Meth_Gene_index)
# 
# #save(Genes_meth_prop,file = 'Genes_meth_prop.txt')
# save(Genes_meth_prop,
#      Genes_meth_select,
#      Diff_Meth_Gene_index,
#      file = 'Genes_meth_prop.rda')

######=========================##########
##         2. in module inves         ##
######=========================##########
# Gather Info: KME and Meth
load('Meth_Prop_Univ.rda')
Genes_meth_prop = Meth_Prop_Univ %>% 
  rename(Count_Prpt = Count_Prpt.y,
         Count_Body = Count_Body.y,
         Count_All = Count_All.y)
library(WGCNA)
table(Gene_net %in% Meth_Prop_Univ$Gene) # lost 299
table(Gene_all %in% Meth_Prop_Univ$Gene) # lost 1722
#length(Gene_all) # 20479
inMod_con = intramodularConnectivity(adjacency_control,
                                     moduleColors_control,
                                     scaleByMax = T) %>% 
  dplyr::select(-kOut,-kDiff)
inMod_con$Gene = rownames(inMod_con)
head(inMod_con)

datKME_tmp = signedKME(datExpr_control,MEs_control)

summary(datKME$Prop_All)
summary(datKME$Prop_Prpt)

'/' <- base:::"/"
datKME = datKME_tmp %>% 
  dplyr::mutate(Gene = rownames(datKME_tmp)) %>% 
  dplyr::mutate(MdouleAssign = moduleColors_control) %>% 
  dplyr::left_join(Genes_meth_prop, by= c("Gene" = "Gene")) %>% 
  dplyr::left_join(inMod_con,by = c('Gene' = 'Gene'))

#
ref=1; test = 2
Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
Zsummary=Z.PreservationStats$Zsummary.pres
#
nonpres_index_b = (which(Zsummary < 2))
nonpres_modulenames_b = rownames(Z.PreservationStats)[nonpres_index_b]
Mod_Index_NonPre  = nonpres_modulenames_b[-grep("grey",nonpres_modulenames_b)]
Mod_Index_Pre = rownames(Z.PreservationStats)[-nonpres_index_b]
Mod_Index_Pre = Mod_Index_Pre[-grep('gold',Mod_Index_Pre)]

######=========================##########
##        diff C prop vs M.M           ##
######=========================##########
#
library(DCGL)
dcg_raw = DCGL::WGCNA((t(datExpr_control)),(t(datExpr_treatment)),
                   power = 24, variant = "DCp")
Diff_Coexp = data.frame(Gene = names(dcg_raw),
                        DCG = dcg_raw,row.names = NULL) %>% 
  left_join(Genes_meth_prop, by = c('Gene' = 'Gene')) %>% 
  dplyr::select(Gene,DCG,Prop_Prpt,Prop_All) %>% 
  drop_na() %>% as_tibble() %>% 
  left_join(inMod_con,by = c('Gene'='Gene')) %>% 
  left_join(Overal_match_color,by = c('Gene'='Gene'))

for (i in seq_along(Diff_Coexp$Gene)){
  target = Diff_Coexp$Gene[i]
  text = unlist(Diff_Coexp[i,7],use.names = F)
  sub = paste('kME',text,sep = '')
  #
  row_loc = which(datKME$Gene == target)
  col_loc = which(names(datKME) == sub)
  # assign
  Diff_Coexp[i,9] = datKME[row_loc,col_loc]
}
names(Diff_Coexp)[9] = 'Module Membership'
Diff_Coexp$`Module Membership` = as.numeric(Diff_Coexp$`Module Membership`)


print(
  ggplot(Diff_Coexp_plot, aes(x = DCG ,y = Prop, colour = Cat)) + # y = Scale_prop
    geom_point() + 
    geom_jitter(width = 0.0001, height = 0.0001,alpha = 1) + 
    geom_smooth(method=lm)+
    ggtitle("Methylation VS DCG Score") +
    xlab("InModule_Con") + ylab("MethC_Prop") +
    guides(color=guide_legend(title="Genomic Regions")) +
    scale_color_manual(labels = c("All","Upper","Body"),
                       values = c("Green","blue", "turquoise"))+
    theme(legend.key = element_rect(fill = "transparent", colour = "transparent"),
          #legend.key = element_rect(colour = 'black', fill = 'blank', size = 0.5, linetype='dashed'),
          legend.position='top', 
          plot.title = element_text(hjust = 0.5)))
dev.off()

print(
  ggplot(Diff_Coexp_plot, aes(x = DCG ,y = Prop, colour = Cat)) + # y = Scale_prop
    geom_point() + 
    geom_jitter(width = 0.0001, height = 0.0001,alpha = 1) + 
    geom_smooth(method=lm)+
    ggtitle("Methylation VS DCG Score") +
    xlab("InModule_Con") + ylab("MethC_Prop") +
    guides(color=guide_legend(title="Genomic Regions")) +
    scale_color_manual(labels = c("All","Upper","Body"),
                       values = c("Green","blue", "turquoise"))+
    theme(legend.key = element_rect(fill = "transparent", colour = "transparent"),
          #legend.key = element_rect(colour = 'black', fill = 'blank', size = 0.5, linetype='dashed'),
          legend.position='top', 
          plot.title = element_text(hjust = 0.5)))
dev.off()



#i = 3
pdf('PDF_Results_NonP_test.pdf')
for (i in seq_along(Mod_Index_NonPre)){
  text = Mod_Index_NonPre[i]
  sub = paste('kME',text,sep = '')
  # subset input data
  data_tmp = subset(datKME,MdouleAssign == text) %>% drop_na() %>% 
    dplyr::mutate( KME = get(sub)) %>%
    # dplyr::mutate(Scale_prop_all = Diff_Prop_all/max(Diff_Prop_all)) %>% 
    # dplyr::mutate(Scale_prop_prmt = Diff_Prop_prom/max(Diff_Prop_prom)) %>% 
    # dplyr::mutate(Scale_prop_body = Diff_Prop_body/max(Diff_Prop_body)) %>% 
    dplyr::select(Gene,KME,kTotal,kWithin,
                  Prop_Prpt,Prop_Body,Prop_All)
  # pivot rotate
  data_module = data_tmp %>%
    pivot_longer(cols = c(Prop_Prpt,Prop_Body,Prop_All),
                 names_to = "Cat", values_to = "Prop")
  print(
    ggplot(data_module, aes(x = kWithin ,y = Prop, colour = Cat)) + # y = Scale_prop
    geom_point() + 
    geom_jitter(width = 0.0001, height = 0.0001,alpha = 1) + 
    geom_smooth(method=lm)+
    ggtitle(paste("Methylation VS ModuleMembership",text,sep='-')) +
    xlab("InModule_Con") + ylab("MethC_Prop") +
    guides(color=guide_legend(title="Genomic Regions")) +
    scale_color_manual(labels = c("All","Upper","Body"),
                       values = c("Green","blue", "turquoise"))+
    theme(legend.key = element_rect(fill = "transparent", colour = "transparent"),
          #legend.key = element_rect(colour = 'black', fill = 'blank', size = 0.5, linetype='dashed'),
          legend.position='top', 
          plot.title = element_text(hjust = 0.5)))
}
dev.off()

pdf('PDF_Results_Pre_test.pdf')
for (i in seq_along(Mod_Index_Pre)){
  text = Mod_Index_Pre[i]
  sub = paste('kME',text,sep = '')
  # subset input data
  data_tmp = subset(datKME,MdouleAssign == text) %>% drop_na() %>% 
    dplyr::mutate( KME = get(sub)) %>%
    # dplyr::mutate(Scale_prop_all = Diff_Prop_all/max(Diff_Prop_all)) %>% 
    # dplyr::mutate(Scale_prop_prmt = Diff_Prop_prom/max(Diff_Prop_prom)) %>% 
    # dplyr::mutate(Scale_prop_body = Diff_Prop_body/max(Diff_Prop_body)) %>% 
    dplyr::select(Gene,KME,kTotal,kWithin,
                  Prop_Prpt,Prop_Body,Prop_All)
  # pivot rotate
  data_module = data_tmp %>%
    pivot_longer(cols = c(Prop_Prpt,Prop_Body,Prop_All),
                 names_to = "Cat", values_to = "Prop")
  print(
    ggplot(data_module, aes(x = kWithin ,y = Prop, colour = Cat)) + # y = Scale_prop
      geom_point() + 
      geom_jitter(width = 0.0001, height = 0.0001,alpha = 1) + 
      geom_smooth(method=lm)+
      ggtitle(paste("Methylation VS ModuleMembership",text,sep='-')) +
      xlab("InModule_Con") + ylab("MethC_Prop") +
      guides(color=guide_legend(title="Genomic Regions")) +
      scale_color_manual(labels = c("All","Upper","Body"),
                         values = c("Green","blue", "turquoise"))+
      theme(legend.key = element_rect(fill = "transparent", colour = "transparent"),
            #legend.key = element_rect(colour = 'black', fill = 'blank', size = 0.5, linetype='dashed'),
            legend.position='top', 
            plot.title = element_text(hjust = 0.5)))
}
dev.off()

######=========================##########
##         3 Hper G test              ##
######=========================##########
setwd('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/Net/')

# library(biomaRt)
# genome <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",  dataset = "btaurus_gene_ensembl", host = "grch37.ensembl.org")
# gene = getBM(c("ensembl_gene_id","external_gene_name", "start_position", "end_position", "chromosome_name"), mart = genome)
# Gene_genome = unique(gene$ensembl_gene_id)
# 
# table(Gene_all %in% Gene_genome)
# '%!in%' <- function(x,y)!('%in%'(x,y))
# table(Gene_net %in% Gene_genome)
# test = Gene_net[which(Gene_net %!in% Gene_genome)]
# head(test)
# 
# Diff_Meth_Gene_index[which(Diff_Meth_Gene_index%in%test)]

# gene index pre
Diff_Meth_Gene_index = unique(Genes_meth_select$Gene)
length(Diff_Meth_Gene_index)

# Analysis and display of module preservation results
load("modulePreservation_methionine.RData")
# test = mp$preservation$Z$ref.control
# stats = mp$preservation$Z$ref.control$inColumnsAlsoPresentIn.treatment
# Results_mp = stats[order(-stats[,2]),c(1:2)]
# head(stats)
# head(Results_mp)

## get module index
ref=1; test = 2;
Z.PreservationStats=mp$preservation$Z[[ref]][[test]];
Zsummary=Z.PreservationStats$Zsummary.pres;nonpres_index_b = (which(Zsummary < 2));
nonpres_modulenames_b = rownames(Z.PreservationStats)[nonpres_index_b]
#
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
#
MP_Stats_medianRank = 
  cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) %>% 
  dplyr::select(medianRank.pres)
MP_Stats_tmp = mp$preservation$Z[[ref]][[test]] 
MP_Stats = merge(MP_Stats_tmp,MP_Stats_medianRank,by="row.names",all.x=TRUE)
MP_Stats_final = MP_Stats[,c(1,2,ncol(MP_Stats),3:(ncol(MP_Stats)-1))] %>% 
  dplyr::filter(!(Row.names %in% c('gold','grey')))
MP_Stats_final[6,2] = 3440
MP_Stats_final


# coral1 3440
# 
#devtools::install_github("slowkow/ggrepel")
library(ggrepel)
set.seed(415)
names(MP_Stats_final)
out = MP_Stats_final[,1:4]
write.csv(out,file = 'MP_Results_Out.csv')
# median rank
MP_Stats_nobig = MP_Stats_final %>% dplyr::filter(!(Row.names %in% c('coral1')))

# 
ggplot(MP_Stats_nobig, aes(moduleSize,medianRank.pres,label = Row.names)) +
  geom_point(size = 5,aes(colour = Row.names))+
  geom_text_repel(vjust = -2)+
  ggtitle("Preservation Median Rank") +
  xlab("ModuleSize") + ylab("Preservation Median Rank") +
  theme(plot.title = element_text(hjust = 0.5,size=16, face="bold.italic"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))

# Preservation Z summary
ggplot(MP_Stats_nobig, aes(moduleSize,Zsummary.pres,label = Row.names)) +
  geom_point(size = 5,aes(colour = Row.names))+
  geom_text_repel(vjust = -2)+
  ggtitle("Preservation Zsummary") +
  xlab("ModuleSize") + ylab("Zsummary") +
  theme(plot.title = element_text(hjust = 0.5,size=16, face="bold.italic"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))+
  geom_abline(intercept = 0,slope=0,colour="black",linetype="dashed")+
  geom_abline(intercept = 2,slope=0,colour="blue",linetype="dashed")+
  geom_abline(intercept = 10,slope=0,colour="darkgreen",linetype="dashed")

# 
Mod_Index_NonPre  = nonpres_modulenames_b[-grep("grey",nonpres_modulenames_b)]
Mod_Index_Pre = rownames(Z.PreservationStats)[-nonpres_index_b]

# massage gene list
Module_assign_all = data.frame(Gene = Gene_net,
                               Assign = moduleColors_control)
UnPreserved_Gene_list = list()
for ( i in seq_along(Mod_Index_NonPre)){
  tmp = as.character(subset(Module_assign_all,Assign == Mod_Index_NonPre[i])$Gene)
  UnPreserved_Gene_list[[i]] = tmp #
  names(UnPreserved_Gene_list)[i]= Mod_Index_NonPre[i]
}
Preserved_Gene_list = list()
Mod_Index_Pre = Mod_Index_Pre[-grep('gold',Mod_Index_Pre)]
for ( i in seq_along(Mod_Index_Pre)){
  tmp = as.character(subset(Module_assign_all,Assign == Mod_Index_Pre[i])$Gene)
  Preserved_Gene_list[[i]] = tmp #
  names(Preserved_Gene_list)[i]= Mod_Index_Pre[i]
}
#
save(UnPreserved_Gene_list,
     Preserved_Gene_list,
     file = 'Gene_list_by_module.rda')
##
load('Gene_list_by_module.rda')

# genes captured in all category
Gene_all = unique(rownames(networkData_normalized))
Gene_net = unique(rownames(networkData_50var_nocrt))
# 
Gene_grey = datKME %>% 
  dplyr::filter(MdouleAssign == 'grey') %>% 
  dplyr::select(Gene) %>% 
  unlist(use.names = F)


# Define differentially methylated genes
Genes_meth_select = Genes_meth_prop %>% 
  dplyr::filter((Count1 >= 20 | Count2 >= 3)) %>% 
  # dplyr::filter(Prop_All>= 0.05) %>%
  # dplyr::filter(Prop_Prpt>= 0.05) %>%
  # dplyr::filter(Prop_Body>= 0.05) %>%
  dplyr::filter(Prop_Prpt >= quantile(Prop_Prpt,0.2)) %>%
  dplyr::filter(Prop_All >= quantile(Prop_All,0.2)) %>%
  dplyr::filter(Prop_Body >= quantile(Prop_Body,0.2)) %>%
  dplyr::filter(Prop_Prpt != 0) %>%
  dplyr::filter(Prop_All != 0) %>%
  dplyr::filter(Prop_Body != 0) %>%
  arrange(Prop_All)

Diff_Meth_Gene_index = unique(Genes_meth_select$Gene)
length(Diff_Meth_Gene_index)
#
Module_Overrep= data.frame(ModuleName = c(),
                            PreservationStatus = c(),
                            #ModuleSize = c(),
                            No.Overlap = c(),
                            Pvalue = c())
Gene_list_all = append(Preserved_Gene_list,UnPreserved_Gene_list)
Grey = list(Grey = Gene_grey)
Gene_list_all = append(Gene_list_all,Grey)
length(Gene_list_all)

for (i in seq_along(names(Gene_list_all))){
  tmp_md = names(Gene_list_all)[i]
  Module_Overrep[i,1] =  tmp_md
  Module_Overrep[i,2] = ifelse(tmp_md %in% Mod_Index_NonPre, 'Sig','Not')
  total.genes =  Gene_all
  sig.genes = Diff_Meth_Gene_index
  gENEs = unlist(Gene_list_all[i]);attributes(gENEs) = NULL
  N = length(total.genes)
  S = length(sig.genes) 
  m = length(total.genes[total.genes %in% gENEs]) # genes from target interpro and in our dataset
  findG = sig.genes[sig.genes %in% gENEs]
  s = length(findG) # # genes from target interpro also in the non-preserved module
  PastefindG = paste(findG, collapse="/")
  M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
  Pval = round(fisher.test(M, alternative ="g")$p.value,100)
  Module_Overrep[i,3] = s
  Module_Overrep[i,4] = Pval
}
Module_Overrep

#
Diff_Coexp_selected = Diff_Coexp %>% 
  dplyr::filter(DCG >= quantile(DCG,0.8)) %>% 
  dplyr::select(Gene) %>% unlist(.,use.names = F)
length(Diff_Coexp_selected)

Module_Overrep_DCG= data.frame(ModuleName = c(), PreservationStatus = c(),
                               #ModuleSize = c(),
                               No.Overlap = c(),
                               Pvalue = c())
sig.genes = Diff_Coexp_selected
for (i in seq_along(names(Gene_list_all))){
  tmp_md = names(Gene_list_all)[i]
  Module_Overrep_DCG[i,1] =  tmp_md
  Module_Overrep_DCG[i,2] = ifelse(tmp_md %in% Mod_Index_NonPre, 'Sig','Not')
  total.genes =  Gene_all
  sig.genes = sig.genes
  gENEs = unlist(Gene_list_all[i]);attributes(gENEs) = NULL
  N = length(total.genes)
  S = length(sig.genes) 
  m = length(total.genes[total.genes %in% gENEs]) # genes from target interpro and in our dataset
  findG = sig.genes[sig.genes %in% gENEs]
  s = length(findG) # # genes from target interpro also in the non-preserved module
  PastefindG = paste(findG, collapse="/")
  M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
  Pval = round(fisher.test(M, alternative ="g")$p.value,100)
  Module_Overrep_DCG[i,3] = s
  Module_Overrep_DCG[i,4] = Pval
}
Module_Overrep_DCG

######=========================##########
##        Plot Enrichment Results     ##
######========================##########
# gene pre
setwd('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/Net/')
library(readxl)
library(tidyverse)
#
load('MESH_Enrichment_0115.RData')
Mesh_results = Parse_Results(Mesh_results_b) %>% 
  dplyr::select(-findG)

ExampleID_Mesh = c("D012269","D025261","D046988")
Mesh_Select = Mesh_results %>% 
  dplyr::filter(MeshID %in% ExampleID_Mesh)
  
# D012269 D025261 D046988

#
load('Interpro_Enrichment_0113.RData')
Interpro_results = Parse_Results(Interpro_results_b)%>% 
  dplyr::select(-findG)

ExampleID_Interpro = c("IPR001353","IPR023333","IPR014722")
Interpro_Select = Interpro_results %>% 
  dplyr::filter(InterproID %in% ExampleID_Interpro)
# IPR001353 IPR023333 IPR014722
names(Interpro_results)

# 
load('Kegg_Enrichment_0113.RData')
Kegg_results = Parse_Results(KEGG_results_b)%>% 
  dplyr::select(-findG)

ExampleID_kegg = c("bta03010","bta03050")
Kegg_Select = Kegg_results %>% 
  dplyr::filter(KEGGID %in% ExampleID_kegg)

# bta03010 bta03050

#
load('Reactome_Enrich_all_path_0113.RData')
Reactome_results = Parse_Results(Reactome_results_b)%>% 
  dplyr::select(-findG)

ExampleID_Reactome = c("R-BTA-72706","R-BTA-72702")
Reactome_Select = Reactome_results %>% 
  dplyr::filter(ReactomeID %in% ExampleID_Reactome)

# R-BTA-72706 R-BTA-72702

# 
load('Msig_Enrichment_0124.RData')
Msig_results = Parse_Results(Results_b)%>% 
  dplyr::select(-findG)
names(Msig_results)

ExampleID_Msig = c("M10085","M17748")
Msig_Select = Msig_results %>% 
  dplyr::filter(MsigID %in% ExampleID_Msig)
# M10085 M17748 


AllID = c(as.character(Reactome_Select[,1]),
          as.character(Msig_Select[,1]),
          as.character(Interpro_Select[,1]),
          as.character(Kegg_Select[,1]),
          as.character(Mesh_Select[,1]))
          
AllPerc = c(as.numeric(Reactome_Select[,8]),
            as.numeric(Msig_Select[,8]),
            as.numeric(Interpro_Select[,8]),
            as.numeric(Kegg_Select[,8]),
            as.numeric(Mesh_Select[,8]))

AllP = c(as.character(Reactome_Select[,5]),
          as.character(Msig_Select[,5]),
          as.character(Interpro_Select[,5]),
          as.character(Kegg_Select[,5]),
          as.character(Mesh_Select[,5]))

AllCount = c(as.character(Reactome_Select[,4]),
          as.character(Msig_Select[,4]),
          as.character(Interpro_Select[,4]),
          as.character(Kegg_Select[,4]),
          as.character(Mesh_Select[,4]))

AllDef = c(as.character(Reactome_Select[,2]),
             as.character(Msig_Select[,2]),
             as.character(Interpro_Select[,2]),
             as.character(Kegg_Select[,2]),
             as.character(Mesh_Select[,2]))
AllDef[3] = 'Release of a polypeptide chain from the ribosome in a mitochondrion'

AllCate = c(rep('Reactome',2),
            rep('Msig',2),
            rep('Interpro',3),
            rep('Kegg',2),
            rep('Mesh',3))

final = data.frame(cbind(AllID,
                         AllPerc,
                         AllP,
                         AllCount,
                         AllCate,
                         AllDef))
final[,2] = as.numeric(as.character(final[,2]))
final[,3] = formatC(as.numeric(as.character(final[,3])),format = "e", digits = 2)
final[,4] = as.numeric(as.character(final[,4]))
final[,6] = as.character(final[,6])
final = final %>% mutate(AllDef_plot = paste(AllDef,AllP,sep = ' - ')) %>% 
  dplyr::rename(DataBase = AllCate, Count =  AllCount)


library(ggrepel)
set.seed(415)
pdf('EnrichOut.pdf')
print(ggplot(final, aes(x = AllPerc,
                  y = AllID,
                  colour = DataBase,
                  size = Count)) +
  geom_point()+
  #geom_text_repel(size = 8)+
  geom_label_repel(data = final[-c(1:2,11),],
                   show.legend = F,
                   size = 3.2,
                   direction = "y",
                   vjust = 0,
                   hjust = - 0.3,
                   aes(label = AllDef_plot))+
  geom_label_repel(data = final[1:2,],
                   #nudge_x = 2,
                   show.legend = F,
                   size = 3.5,
                   direction = "y",
                   hjust = 1.2,
                   vjust = -0.2,
                   aes(label = AllDef_plot))+
  geom_label_repel(data = final[11,],
                   #nudge_x = 2,
                   show.legend = F,
                   size = 3.5,
                   direction = "y",
                   hjust = 1.1,
                   vjust = 0,
                   aes(label = AllDef_plot))+
  ggtitle("Overrepresentation in Module antiquewhite2") +
  xlab("Hit Percentage") + 
  ylab("Biological Pathway ID") +
  theme_gray()+
  theme(plot.title = element_text(hjust = 0.5,size=16, face="bold.italic"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold")))
dev.off()

##
# TransFactor_raw = read.csv('Bos_taurus_TF.txt',sep = "\t")
# TransFactor_co_raw = read.csv('Bos_taurus_TF_cofactors.txt',sep = "\t")
# head(TransFactor_co_raw)

ov1 = test[which(test %in% Gene_list_all[[14]])]
ov2 = Diff_Meth_Gene_index[which(Diff_Meth_Gene_index %in% test)]
table(ov1 %in% ov2)

for (i in 1:14){
  tmp = Gene_list_all[[i]]
  x = table(test %in% tmp)
  print(length(tmp))
  print(x)
}

######=========================##########
##        4. Transcription factor     ##
######=========================##########

# devtools::install_github("slowkow/tftargets")
library(tftargets)
#detach("package:tftargets", unload=TRUE)
length(TRED)#etz
length(ENCODE)#etz
length(TRRUST)
length(ITFP)
length(Neph2012)#etz
length(Marbach2016)

check_symbol = function(notsure){
  tmp = checkGeneSymbols(notsure) %>% 
    mutate(Suggested.Symbol.Merge = ifelse(is.na(Suggested.Symbol),x,Suggested.Symbol))
  out = unlist(tmp$Suggested.Symbol.Merge,use.names = F)
  return(out)
}
# out = check_symbol(symbol_test_human$Symbol)


# retrive human gene information
library(org.Hs.eg.db);library(org.Bt.eg.db);library(tidyverse)
entrezUniverse_human = AnnotationDbi::select(org.Hs.eg.db, 
                                             as.character(AnnotationDbi::keys(org.Hs.eg.db,keytype = c("ENSEMBL"))), 
                                             columns = c("ENTREZID","SYMBOL"),
                                             keytype = "ENSEMBL") %>% 
  dplyr::distinct(ENSEMBL,.keep_all= TRUE)
save(entrezUniverse_human,file = 'Human_gene.rda')
load('Human_gene.rda')

# find if symbols are official in human database
library(HGNChelper);library(tidyverse)
# for human
symbol_test_human = 
  unique(entrezUniverse_human$SYMBOL) %>% 
  .[.!= ""] %>% 
  checkGeneSymbols(.) %>% as_tibble() %>% 
  dplyr::filter(Approved != 'TRUE',!(is.na(Suggested.Symbol))) %>% 
  rename(Symbol = x) %>% 
  dplyr::select(Symbol,Suggested.Symbol) # 207 records
head(entrezUniverse_human)
#
entrezUniverse_human_corrected = entrezUniverse_human %>% 
  mutate(Suggested.Symbol = check_symbol(SYMBOL)) %>% 
  dplyr::select(-SYMBOL)

# for bovine
load('gene_symbols_genome.rda')
symbol_test_bovine = 
  unique(gene_symbols_genome$external_gene_name) %>% 
  .[.!= ""] %>% 
  checkGeneSymbols(.) %>% as_tibble() %>% 
  dplyr::filter(Approved != 'TRUE',!(is.na(Suggested.Symbol))) %>% 
  rename(Symbol = x) %>% 
  dplyr::select(Symbol,Suggested.Symbol)

symbol_test_bovine_corrected = gene_symbols_genome %>% 
  mutate(Suggested.Symbol = check_symbol(external_gene_name)) %>% 
  dplyr::select(-external_gene_name)

# for database
# symbol - needs trans
length(TRRUST)
length(ITFP)
length(Marbach2016)
#
symbollist = c('TRRUST','ITFP','Marbach2016')
symbollist_out = c('TRRUST_cr','ITFP_cr','Marbach2016_cr')
#
for (i in seq_along(symbollist)){
  tmp = get(symbollist[i])
  name_new = check_symbol(names(tmp))
  for (p in seq_along(name_new)){
    names(tmp)[p] = name_new[p]
    just4print = check_symbol(tmp[[p]])
    #print(length(which(!(tmp[[p]]%in%just4print))))
    print(c(i,p))
    tmp[[p]] = just4print
  }
  assign(symbollist_out[i],tmp)
}
# needs matching
length(TRED)#etz
length(ENCODE)#etz
length(Neph2012)#etz
SHMM = Neph2012[[29]]
entrezlist = c('TRED','ENCODE','SHMM')
entrezlist_out = c('TRED_cr','ENCODE_cr','SHMM_cr')
#
#count = c()
#i = 1; p = 1
for (i in seq_along(entrezlist)){
  tmp = get(entrezlist[i])
  name_new = check_symbol(names(tmp))
  diff = c()
  for (p in seq_along(name_new)){
    names(tmp)[p] = name_new[p]
    entrez = tmp[[p]]
    entrez2symbol_raw = entrezUniverse_human_corrected %>% 
      dplyr::filter(ENTREZID %in% entrez) %>% 
      dplyr::select(Suggested.Symbol) %>% 
      unlist(use.names = F)
    entrez2symbol = unique(entrez2symbol_raw)
    #diff[p] = length(entrez2symbol_raw) - length(entrez2symbol)
    tmp[[p]] = entrez2symbol
    print(c(i,p))
  }
  count = append(count,diff)
  assign(entrezlist_out[i],tmp)
}
#test = entrezUniverse_human_corrected %>% dplyr::filter(ENTREZID %in% TRED[[12]]) %>% dplyr::select(ENTREZID,Suggested.Symbol.Merge)
save(TRED_cr,
     ENCODE_cr,
     SHMM_cr,
     TRRUST_cr,
     ITFP_cr,
     Marbach2016_cr,file = 'TF_data_cr.rda')

# 
tf_database_index  = c('TRED_cr','ENCODE_cr','TRRUST_cr','ITFP_cr','SHMM_cr','Marbach2016_cr')
# 

# container: given a module(gene set); search for overlap for each TF in all databases, record the overlap
library(tidyverse)
# Data Structure Preperation
Gene_list_all = append(UnPreserved_Gene_list,Preserved_Gene_list)
ModuleName_Unpreserved = names(UnPreserved_Gene_list)
ModuleName_Preserves = names(Preserved_Gene_list)
# gene convert output container
GeneConvert = data.frame(ModuleName = c(),
                         Preseravtion = c(),
                         GeneCount = c(),
                         DupCount = c(), # only rm na
                         FinalCount = c()) # rm dup
#
special_intersect = function(test_list,
                             database_list,
                             pattern = ' /// '){
  Try_Multi = function(test_list,database_list,pattern = ' /// '){
    tmp = unlist(str_split(database_list, pattern = pattern))
    out = intersect(tmp,test_list)
    return(out)
  }
  out_list = intersect(test_list,database_list)
  Multi_index = grepl(pattern,database_list)
  Multi_out = database_list[Multi_index]
  for (i in seq_along(Multi_out)){
    findone = Try_Multi(test_list,Multi_out[i])
    out_list = append(out_list,findone,after = length(out_list))
  }
  count = length(out_list)
  return(out_list)
}

# output container
OverlapOut = data.frame(DataBase = c(),
                        TF_Name = c(),
                        OverNum_Ensl = c(),
                        OverGene = c(),
                        Module = c())
#p = 1;i = 1 ;j = 1
for (p in seq_along(Gene_list_all)){
  OUT =  data.frame(DataBase = c(),
                    TF_Name = c(),
                    OverNum_Ensl = c(),
                    OverGene = c(),
                    Module = c())
  # convert gene id to symbols
  tmp = (Gene_list_all)[[p]]
  md_name_tmp = names(Gene_list_all)[p]
  GeneConvert[p,1] = md_name_tmp
  GeneConvert[p,2] = ifelse(md_name_tmp %in% ModuleName_Unpreserved,'Unpreserved','Preserved')
  GeneConvert[p,3] = length(tmp)
  gene_collection_tmp_raw = symbol_test_bovine_corrected %>% 
    dplyr::filter(ensembl_gene_id %in% tmp) %>% 
    dplyr::select(Suggested.Symbol) %>% 
    drop_na() %>% unlist(use.names = F) %>% 
    .[nchar(.)>0]
  gene_collection_tmp = gene_collection_tmp_raw %>% unique()
  GeneConvert[p,4] = length(gene_collection_tmp_raw)
  GeneConvert[p,5] = length(gene_collection_tmp)
  #
  Module.Gene.Input = gene_collection_tmp
  for (i in seq_along(tf_database_index)){
    out_tmp = data.frame(DataBase = c(),
                         TF_Name = c(),
                         OverNum = c(),
                         OverGene = c(),
                         Module = c())
    DB_name = tf_database_index[i]
    tmp_TF_index = names(get(DB_name))
    for (j in seq_along(tmp_TF_index)){
      tmp_TF_gene_list = get(DB_name)[[j]]
      theoverlap = special_intersect(gene_collection_tmp,tmp_TF_gene_list)
      # out_tmp[i+j-1,1] = DB_name
      # out_tmp[i+j-1,2] = tmp_TF_index[j]
      # out_tmp[i+j-1,3] = length(theoverlap)
      # out_tmp[i+j-1,4] = ifelse(length(theoverlap) == 0,0,length(theoverlap))
      # out_tmp[i+j-1,5] = md_name_tmp
      SingleRecord = data.frame(DataBase = DB_name,
                                TF_Name = tmp_TF_index[j],
                                OverNum = length(theoverlap),
                                OverGene = ifelse(length(theoverlap) == 0,0,length(theoverlap)),
                                Module = md_name_tmp)
      out_tmp = rbind(out_tmp,SingleRecord)
      print(c(p,i,j))
    }
    OUT = rbind(OUT,out_tmp)
  }
  OverlapOut = rbind(OverlapOut,OUT)
}
dim(OverlapOut);head(OverlapOut)
table(OverlapOut$OverNum)

#
library(readxl)
library(openxlsx)
# select TF that are (i) high confid (class a/b) (ii) Expressed in Muscle (iii) not necessary to have a official symbol
Bta_TF_raw = read.xlsx('Bta_TF_list.xlsx',sheet = 4,startRow = 6) %>% 
  dplyr::select(Bovine_Class,Ensembl_ID,TF_Symbol,Tissue_expression) %>% 
  dplyr::filter(Bovine_Class %in% c('a','b')) %>% 
  dplyr::filter(str_detect(Tissue_expression,'Muscle',negate = F)) %>% 
  dplyr::mutate(Suggested.Symbol = check_symbol(TF_Symbol)) %>% 
  dplyr::mutate(SymbExist = ifelse(TF_Symbol == c('.'),'No','Yes')) %>% 
  #dplyr::filter(!(TF_Symbol %in% c('.'))) %>% 
  dplyr::select(-Bovine_Class,-Tissue_expression,-TF_Symbol)
head(Bta_TF_raw) # resulting 444 itemsBta_TF_raw

Bta_TF_list = Bta_TF_raw %>% 
  dplyr::select(Ensembl_ID) %>% unlist(use.names = F) %>% unique()
length(Bta_TF_list)


# TcoF
Bta_TcoF_raw = read.xlsx('Bta_TF_list.xlsx',sheet = 5,startRow = 6) %>% 
  dplyr::select(TcoF_ID,TF_ID,TcoF_Class,TcoF_tissue_expression) %>% 
  dplyr::filter(TcoF_Class %in% c('High')) %>% 
  dplyr::filter(str_detect(TcoF_tissue_expression,'Muscle',negate = F)) %>% 
  dplyr::filter(TF_ID %in% Bta_TF_raw$Ensembl_ID) %>% 
  dplyr::select(-TcoF_tissue_expression,-TcoF_Class)
Bta_TcoF_list = Bta_TcoF_raw %>% 
  dplyr::select(TcoF_ID) %>% unlist(use.names = F) %>% unique()
length(Bta_TcoF_list)


Bta_TF_Mstatus_raw = Bta_TF_raw %>% 
  left_join(Bta_TcoF_raw, by = c('Ensembl_ID' = 'TF_ID')) %>% 
  left_join(dplyr::select(symbol_test_bovine_corrected,ensembl_gene_id,Suggested.Symbol),
            by = c('TcoF_ID' = 'ensembl_gene_id')) %>% 
  left_join(Genes_meth_prop,by = c('Ensembl_ID'='Gene')) %>% 
  left_join(Genes_meth_prop,by = c('TcoF_ID'='Gene'))

#
Bta_TF_Mstatus_final = Bta_TF_Mstatus_raw %>%  
  dplyr::select(-c('Count_Prpt.x','Count_Body.x','Count_All.x',
                   'Count_Prpt.y','Count_Body.y','Count_All.y',)) %>% 
  group_by(Suggested.Symbol.x) %>% 
  mutate(TcoF_meth = mean(Prop_All.y)) %>% 
  dplyr::select(-Prop_Body.y,-Prop_Prpt.y,-Prop_All.y,-SymbExist,-TcoF_ID,-Suggested.Symbol.y) %>% 
  sample_n(1) %>% 
  mutate(Prop_All_wTcoF = ifelse(is.na(Prop_All.x),TcoF_meth,TcoF_meth + Prop_All.x))

#view(Bta_TF_Mstatus_final)

#
ModuleName_Unpreserved = names(UnPreserved_Gene_list)
ModuleName_Preserves = names(Preserved_Gene_list)
ModuleSize = data.frame(Module = names(Gene_list_all),
                        Size = sapply(Gene_list_all, length))
Overal_match_color = data.frame(Gene = colnames(datExpr_control),
                                Module = moduleColors_control) %>% 
  mutate(Cate = ifelse(Module %in% ModuleName_Unpreserved,'Unp','Pre'))

all_Bta_TF_muscle = as.character(unique(Bta_TF_Mstatus_final$Suggested.Symbol.x))

Genes_meth_prop_withSymbol = 
  Genes_meth_prop %>% 
  left_join(symbol_test_bovine_corrected,by=c('Gene'='ensembl_gene_id'))

Bta_TF_OverlapMatch = OverlapOut %>% 
  dplyr::filter(OverNum != 0) %>%
  dplyr::filter(TF_Name %in% all_Bta_TF_muscle) %>% 
  dplyr::filter(TF_Name %in% Genes_meth_prop_withSymbol$Suggested.Symbol) %>%
  dplyr::filter(!(Module %in% 'Grey')) %>% 
  left_join(ModuleSize, by = c('Module' = 'Module')) %>% 
  mutate(Pres = ifelse(Module %in% ModuleName_Unpreserved,'Unpre','Pre')) %>% 
  mutate(OverPerc = OverNum/Size) %>% as_tibble() %>% 
  group_by(TF_Name) %>% 
  #dplyr::filter(OverNum == max(OverNum)) %>% sample_n(1) %>% 
  dplyr::filter(OverPerc == max(OverPerc)) %>% sample_n(1) #%>% 
  dplyr::filter(DataBase != 'Marbach2016_cr')

#
Bta_TF_OverlapMatch_plot = Bta_TF_OverlapMatch %>% 
  left_join(Bta_TF_Mstatus_final,by = c('TF_Name'='Suggested.Symbol.x')) %>% 
  mutate(Prop_All_Final = ifelse(is.na(Prop_All_wTcoF),Prop_All.x,Prop_All_wTcoF))
# reorder
reindex = c(which(Bta_TF_OverlapMatch_plot$Pres == 'Pre'),which(Bta_TF_OverlapMatch_plot$Pres == 'Unpre'))
Bta_TF_OverlapMatch_plot =Bta_TF_OverlapMatch_plot[reindex,]
Bta_TF_OverlapMatch_plot$Module = factor(Bta_TF_OverlapMatch_plot$Module,
                                         levels = c(ModuleName_Unpreserved,ModuleName_Preserves))

# melt
library(reshape2)
Bta_TF_OverlapMatch_plot.m =
  melt(Bta_TF_OverlapMatch_plot,
       id.vars='TF_Name', 
       measure.vars=c('Count1.x','Count2.x')) %>% 
  left_join(Bta_TF_OverlapMatch_plot,by = c('TF_Name'='TF_Name'))

dim(Bta_TF_OverlapMatch_plot.m)

class(Bta_TF_OverlapMatch_plot$Module)
view(Bta_TF_OverlapMatch_plot)
names(Bta_TF_OverlapMatch_plot)
library(ggplot2)

#Overall CpG Counts of TF(TcoF) in every module
ggplot(Bta_TF_OverlapMatch_plot,
       aes(x = Module,
           y=(Count1.x + Count2.x),
           fill=Pres)) + 
  geom_boxplot(alpha =.4,width = .2)+
  geom_violin(alpha =.8,width = .8) +
  facet_wrap(~Pres,scales = 'free',
             labeller = labeller(Pres = c(Pre = "Preserved", Unpre = "Unpreserved")))+
  theme(axis.text.x = element_text(face = "bold",size = 7,angle = 45),
        legend.position = "none")+
  coord_flip()+
  xlab("Modules") + ylab("Significant CpG Counts (Total)")

# CpG Counts by large categories

Bta_TF_OverlapMatch_plot_2 = Bta_TF_OverlapMatch_plot %>% 
  group_by(TF_Name) %>% 
  dplyr::filter(OverPerc == max(OverPerc)) %>% 
  dplyr::select(TF_Name,Count1.x,Count2.x,Prop_Body.x,Prop_Prpt.x,Prop_All_Final,Pres)

#Gene_Traced_all
Bta_TF_Mstatus_final_2 = Bta_TF_Mstatus_final %>% 
  mutate(Prop_All_Final = ifelse(is.na(Prop_All_wTcoF),Prop_All.x,Prop_All_wTcoF)) %>% 
  dplyr::select(Suggested.Symbol.x,Count1.x,Count2.x,Prop_Body.x,Prop_Prpt.x,Prop_All_Final) %>% 
  drop_na() %>% 
  mutate(Pres = 'All') %>% 
  rename(TF_Name = Suggested.Symbol.x)
#
LargeCate_plot1 = rbind(Bta_TF_OverlapMatch_plot_2,
                       Bta_TF_Mstatus_final_2)
#
ggplot(LargeCate_plot1,
       aes(x = Pres,y =(Count1.x + Count2.x),fill=Pres))+
  geom_violin(alpha =.8,width = .8) +
  xlab("Module Preservation") + ylab("Prop_ALL (with TcoF)")

ggplot(LargeCate_plot1,
       aes(x = Pres,y =Prop_All_Final,fill=Pres))+
  geom_violin(alpha =.8,width = .8) +
  xlab("Module Preservation") + ylab("Prop_ALL (with TcoF)")

#
# all genes and pre/unpre genes with DNA methy measure
# surprisingly, variance relates to methylation counts
Genes_with_Meth = Genes_meth_prop %>% 
  dplyr::filter(!(Count1 == 0 & Count1 == 0)) %>% 
  mutate(Cate = 'All')

LargeCate_plot2_raw = Genes_meth_prop %>% 
  dplyr::filter(!(Count1 == 0 & Count1 == 0),Gene %in% Gene_net) %>% 
  left_join(dplyr::select(Overal_match_color,Gene,Cate),
                          by = c('Gene'='Gene'))

LargeCate_plot2 = rbind(LargeCate_plot2_raw,Genes_with_Meth)


ggplot()+
  geom_violin(data = LargeCate_plot2,
              aes(x = Cate,y =Prop_Prpt,fill=Cate))+
  #geom_flat_violin(position=position_nudge(x=.2)) +
  #geom_violin(alpha =.8,width = .8) +
  #geom_boxplot(outlier.alpha = 0.2,outlier.size = 0.3) +
  xlab("Module Preservation") + ylab("Prop_ALL (with TcoF)")
  #coord_flip()

#
require(ggpubr)
ggarrange(a,b + rremove("x.text"), 
          labels = c("A", "B"),
          ncol = 2)

#Bta_TF_OverlapMatch_plot,aes(x = Module,y=(Count1.x + Count2.x),fill=Pres)



test = Genes_meth_prop %>% 
  left_join(Overal_match_color, by = c('Gene' ='Gene')) %>% 
  mutate(Pres = ifelse(Module %in% ModuleName_Unpreserved,'Unpre','Pre')) %>% 
  filter(Count1 == 0 | Count2 == 0) %>% 
  filter(Gene %in% Gene_net) %>% 
  dplyr::select(Gene) %>% unlist(use.names = F)

outtestunpre = Genes_meth_prop_final[which(Genes_meth_prop_final$Pres =='Unpre'),]
outtestpre = Genes_meth_prop_final[which(Genes_meth_prop_final$Pres =='Pre'),]
boxplot(outtestunpre$Prop_All,
        outtestpre$Prop_All)


Bta_TF_Meth_plot = Bta_TF_OverlapMatch %>% 
  group_by(TF_Name) %>% 
  dplyr::filter(OverNum == max(OverNum))

table(Bta_TF_Meth_plot$Pres)
table(Bta_TF_Meth_plot$Module)

left_join(Bta_TF_Mstatus_final,by = c('TF_Name' = 'Suggested.Symbol.x')) %>% 
  group_by(TF_Name) %>% 
  replace_na(list(Prop_All.x = 0, 
                  Prop_Body.x = 0,
                  Prop_Prpt.x = 0,
                  Prop_All.y = 0, 
                  Prop_Body.y = 0,
                  Prop_Prpt.y = 0)) %>% 
  mutate(Prop_All_comb = Prop_All.x + mean(Prop_All.y)) %>% 
  mutate(Prop_Body_comb = Prop_Body.x + mean(Prop_Body.y)) %>% 
  mutate(Prop_Prpt_comb = Prop_Prpt.x + mean(Prop_Prpt.y)) %>% 
  dplyr::select(-c('Prop_All.x','Prop_Body.x','Prop_Prpt.x',
                   'Prop_All.y','Prop_Body.y','Prop_Prpt.y'))
Bta_TF_Meth_plot_final = Bta_TF_Meth_plot %>% 
  group_by(TF_Name) %>%
  slice(1)

view(Bta_TF_Meth_plot_final)
table(Bta_TF_Meth_plot_final$Pres)
table(Bta_TF_Meth_plot_final$Module)

# sig
Bta_TF_match_sig = Bta_TF_match %>% 
  dplyr::filter(OverGene >= 0.4 * Size) #%>% filter(DataBase != "Marbach2016_cr")

#
Gene_list_all = append(UnPreserved_Gene_list,Preserved_Gene_list)
ModuleName_Unpreserved = names(UnPreserved_Gene_list)
ModuleName_Preserves = names(Preserved_Gene_list)
Gene_list_Unpre = c()
Gene_list_Pre = c()
for (i in names(UnPreserved_Gene_list)){
  tmp = unlist(UnPreserved_Gene_list[[i]],use.names = F)
  Gene_list_Unpre = append(Gene_list_Unpre,tmp,length(Gene_list_Unpre))
}
for (i in names(Preserved_Gene_list)){
  tmp = unlist(Preserved_Gene_list[[i]],use.names = F)
  Gene_list_Pre = append(Gene_list_Pre,tmp,length(Gene_list_Pre))
}


######=========================##########
##        Hyper G test                ##
######========================##########
setwd('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/Net/')

# library(biomaRt)
# genome <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",  dataset = "btaurus_gene_ensembl", host = "grch37.ensembl.org")
# gene = getBM(c("ensembl_gene_id","external_gene_name", "start_position", "end_position", "chromosome_name"), mart = genome)
# Gene_genome = unique(gene$ensembl_gene_id)
# 
# table(Gene_all %in% Gene_genome)
# '%!in%' <- function(x,y)!('%in%'(x,y))
# table(Gene_net %in% Gene_genome)
# test = Gene_net[which(Gene_net %!in% Gene_genome)]
# head(test)
# 
# Diff_Meth_Gene_index[which(Diff_Meth_Gene_index%in%test)]

# gene index pre
Diff_Meth_Gene_index = unique(Genes_meth_select$Gene)
length(Diff_Meth_Gene_index)

# Analysis and display of module preservation results
load("modulePreservation_methionine.RData")
# test = mp$preservation$Z$ref.control
# stats = mp$preservation$Z$ref.control$inColumnsAlsoPresentIn.treatment
# Results_mp = stats[order(-stats[,2]),c(1:2)]
# head(stats)
# head(Results_mp)

Mod_Index_NonPre  = nonpres_modulenames_b[-grep("grey",nonpres_modulenames_b)]
Mod_Index_Pre = rownames(Z.PreservationStats)[-nonpres_index_b]

# massage gene list
Module_assign_all = data.frame(Gene = Gene_net,
                               Assign = moduleColors_control)
UnPreserved_Gene_list = list()
for ( i in seq_along(Mod_Index_NonPre)){
  tmp = as.character(subset(Module_assign_all,Assign == Mod_Index_NonPre[i])$Gene)
  UnPreserved_Gene_list[[i]] = tmp #
  names(UnPreserved_Gene_list)[i]= Mod_Index_NonPre[i]
}
Preserved_Gene_list = list()
Mod_Index_Pre = Mod_Index_Pre[-grep('gold',Mod_Index_Pre)]
for ( i in seq_along(Mod_Index_Pre)){
  tmp = as.character(subset(Module_assign_all,Assign == Mod_Index_Pre[i])$Gene)
  Preserved_Gene_list[[i]] = tmp #
  names(Preserved_Gene_list)[i]= Mod_Index_Pre[i]
}
#
save(UnPreserved_Gene_list,
     Preserved_Gene_list,
     file = 'Gene_list_by_module.rda')
##
load('Gene_list_by_module.rda')

names(Preserved_Gene_list)
# genes captured in all category
Gene_all = unique(rownames(networkData_normalized))
Gene_net = unique(rownames(networkData_50var_nocrt))
# 
Gene_grey = datKME %>% 
  dplyr::filter(MdouleAssign == 'grey') %>% 
  dplyr::select(Gene) %>% 
  unlist(use.names = F)


# Define differentially methylated genes

Genes_meth_select = Genes_meth_prop %>% 
  dplyr::filter((Count1 >= 20 | Count2 >= 3)) %>% 
  # dplyr::filter(Prop_All>= 0.05) %>%
  # dplyr::filter(Prop_Prpt>= 0.05) %>%
  # dplyr::filter(Prop_Body>= 0.05) %>%
  dplyr::filter(Prop_Prpt >= quantile(Prop_Prpt,0.2)) %>%
  dplyr::filter(Prop_All >= quantile(Prop_All,0.2)) %>%
  dplyr::filter(Prop_Body >= quantile(Prop_Body,0.2)) %>%
  dplyr::filter(Prop_Prpt != 0) %>%
  dplyr::filter(Prop_All != 0) %>%
  dplyr::filter(Prop_Body != 0) %>%
  arrange(Prop_All)

Diff_Meth_Gene_index = unique(Genes_meth_select$Gene)
length(Diff_Meth_Gene_index)
#
Module_Overrep= data.frame(ModuleName = c(),
                           PreservationStatus = c(),
                           #ModuleSize = c(),
                           No.Overlap = c(),
                           Pvalue = c())
Gene_list_all = append(Preserved_Gene_list,UnPreserved_Gene_list)
Grey = list(Grey = Gene_grey)
Gene_list_all = append(Gene_list_all,Grey)
length(Gene_list_all)

for (i in seq_along(names(Gene_list_all))){
  tmp_md = names(Gene_list_all)[i]
  Module_Overrep[i,1] =  tmp_md
  Module_Overrep[i,2] = ifelse(tmp_md %in% Mod_Index_NonPre, 'Sig','Not')
  total.genes =  Gene_all
  sig.genes = Diff_Meth_Gene_index
  gENEs = unlist(Gene_list_all[i]);attributes(gENEs) = NULL
  N = length(total.genes)
  S = length(sig.genes) 
  m = length(total.genes[total.genes %in% gENEs]) # genes from target interpro and in our dataset
  findG = sig.genes[sig.genes %in% gENEs]
  s = length(findG) # # genes from target interpro also in the non-preserved module
  PastefindG = paste(findG, collapse="/")
  M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
  Pval = round(fisher.test(M, alternative ="g")$p.value,100)
  Module_Overrep[i,3] = s
  Module_Overrep[i,4] = Pval
}
Module_Overrep

#
Diff_Coexp_selected = Diff_Coexp %>% 
  dplyr::filter(DCG >= quantile(DCG,0.8)) %>% 
  dplyr::select(Gene) %>% unlist(.,use.names = F)
length(Diff_Coexp_selected)

Module_Overrep_DCG= data.frame(ModuleName = c(), PreservationStatus = c(),
                               #ModuleSize = c(),
                               No.Overlap = c(),
                               Pvalue = c())
sig.genes = Diff_Coexp_selected
for (i in seq_along(names(Gene_list_all))){
  tmp_md = names(Gene_list_all)[i]
  Module_Overrep_DCG[i,1] =  tmp_md
  Module_Overrep_DCG[i,2] = ifelse(tmp_md %in% Mod_Index_NonPre, 'Sig','Not')
  total.genes =  Gene_all
  sig.genes = sig.genes
  gENEs = unlist(Gene_list_all[i]);attributes(gENEs) = NULL
  N = length(total.genes)
  S = length(sig.genes) 
  m = length(total.genes[total.genes %in% gENEs]) # genes from target interpro and in our dataset
  findG = sig.genes[sig.genes %in% gENEs]
  s = length(findG) # # genes from target interpro also in the non-preserved module
  PastefindG = paste(findG, collapse="/")
  M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
  Pval = round(fisher.test(M, alternative ="g")$p.value,100)
  Module_Overrep_DCG[i,3] = s
  Module_Overrep_DCG[i,4] = Pval
}
Module_Overrep_DCG

# 
# sig
Bta_TF_match_sig = Bta_TF_match %>% 
  dplyr::filter(OverGene >= 0.2 * Size) #%>% filter(DataBase != "Marbach2016_cr")
#
Bta_TF_match_sig %>% group_by(Module) %>% count()
#names(Bta_TF_match_sig)
#
Bta_TF_match_sig_select = Bta_TF_match_sig %>% 
  dplyr::select(Ensembl_ID) %>% unlist(use.names = F) %>% unique()
length(Bta_TF_match_sig_select)
#
Module_Overrep_TF= data.frame(ModuleName = c(), PreservationStatus = c(),
                               #ModuleSize = c(),
                               No.Overlap = c(),
                               Pvalue = c())
sig.genes = Bta_TF_list
for (i in seq_along(names(Gene_list_all))){
  tmp_md = names(Gene_list_all)[i]
  Module_Overrep_TF[i,1] =  tmp_md
  Module_Overrep_TF[i,2] = ifelse(tmp_md %in% Mod_Index_NonPre, 'Sig','Not')
  total.genes =  Gene_all
  sig.genes = sig.genes
  gENEs = unlist(Gene_list_all[i]);attributes(gENEs) = NULL
  N = length(total.genes)
  S = length(sig.genes) 
  m = length(total.genes[total.genes %in% gENEs]) # genes from target interpro and in our dataset
  findG = sig.genes[sig.genes %in% gENEs]
  s = length(findG) # # genes from target interpro also in the non-preserved module
  PastefindG = paste(findG, collapse="/")
  M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
  Pval = round(fisher.test(M, alternative ="g")$p.value,100)
  Module_Overrep_TF[i,3] = s
  Module_Overrep_TF[i,4] = Pval
}
Module_Overrep_TF
