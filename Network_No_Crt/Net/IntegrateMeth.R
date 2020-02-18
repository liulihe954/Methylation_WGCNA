######=========================##########
##            Data Pre                 ##
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


# counted ALL Cs 
# measure meth level
# run rgmatch find the Diff C location/ gene assignments
## ================================================================================================================== ##
#     python rgmatch.py -g Bos_taurus.ARS-UCD1.2.99.gtf -b Diff_C_Sig_BED.bed -r 'gene' -q 5 -o myassociations.txt    ##
## ================================================================================================================== ##
Data_loci = '/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net'
#setwd('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/Net')
setwd('/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net/rgmatch')
#Associ_out_raw = read.table('myassociations_gene_new.txt',sep = '\t') %>% data.frame()
Associ_out_raw = read.table('myassociations_exon_new.txt',sep = '\t') %>% data.frame()
setwd(Data_loci)

### massage diff C associated gene
cname =as.character(unlist(Associ_out_raw[1,]));attributes(cname) = NULL
colnames(Associ_out_raw) = cname
Associ_out_raw= Associ_out_raw[-1,]
####
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
#head(Associ_out,100) %>% print(n = Inf)
Associ_out =  Associ_out_raw %>% 
  dplyr::select(-Distance,-Transcript,-`Exon/Intron`,-TSSDistance,-PercRegion,-PercArea) %>% 
  group_by(Gene) %>% 
  dplyr::count(Area) %>%
  tidyr::spread(key = Gene, value = n) %>% 
  transpose_df() %>% 
  replace(is.na(.), 0) %>% 
  mutate_at(vars(-Area), as.numeric) %>% 
  dplyr::rename(Gene = Area)

Associ_out_count = Associ_out %>% 
  mutate(Count1 = `1st_EXON`+ GENE_BODY + INTRON) %>% 
  mutate(SigG1 = ifelse(Count1 > 30, "Sig", "Not")) %>% 
  mutate(Count2 = PROMOTER + TSS + UPSTREAM) %>% 
  mutate(SigG2 = ifelse(Count2 > 5, "Sig", "Not"))

DiffC2Gene_sig = Associ_out_count %>% 
  dplyr::filter(SigG1 == 'Sig' | SigG2 == 'Sig' )

# sig meth gene whole genome
Gene_DiffC_index = unique(DiffC2Gene_sig$Gene)

# count all Cs using self-wraped function: AGCTcount (getSEQ from biomart mainly)
load('Genes_C_count_all_Final.RData')
Genes_C_count_all = Genes_C_count_all %>% 
  mutate_at(vars(Total_CG),as.numeric) %>% 
  replace(is.na(.), 0)

## calculate proption of diff Cs
"/" <- function(x,y) ifelse(y==0,0,base:::"/"(x,y)) # special division
Genes_meth_prop = Genes_C_count_all %>% 
  dplyr::left_join(Associ_out_count, by = c('Gene' = 'Gene')) %>% 
  dplyr::select(Gene,Total_CG,Count1,Count2) %>% 
  replace(is.na(.), 0) %>% 
  mutate(Count_all = Count1 + Count2) %>% 
  mutate(Diff_Prop = Count_all/Total_CG)
#save(Genes_meth_prop,file = 'Genes_meth_prop.txt')
save(Genes_meth_prop,file = 'Genes_meth_prop.rda')

#dim(Genes_meth_prop)

# gene index pre
Gene_all_rna = unique(rownames(networkData_normalized))
Gene_all_net = unique(rownames(networkData_50var_nocrt)) 

# Gather Info: KME and Meth
load('Genes_meth_prop.rda')
datKME_tmp = signedKME(datExpr_control, MEs_control)
datKME = datKME_tmp %>% 
  dplyr::mutate(Gene = rownames(datKME_tmp)) %>% 
  dplyr::mutate(MdouleAssign = moduleColors_control) %>% 
  dplyr::left_join(Genes_meth_prop, by= c("Gene" = "Gene")) %>% 
  dplyr::mutate(Scale_prop = Diff_Prop/max(Diff_Prop)) #%>% 
  d%>% plyr::filter(Scale_prop != 1)
head(datKME)
summary(datKME$Scale_prop)

######=========================##########
##            Index Pre                ##
######=========================##########
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
pdf('PDF_Results_NonP.pdf')
for (i in seq_along(Mod_Index_NonPre)){
  text = Mod_Index_NonPre[i]
  sub = paste('kME',text,sep = '')
  print(ggplot(datKME,aes(x=get(sub),y = Scale_prop)) + 
    geom_point(colour="grey") +
    geom_point(data = subset(datKME, MdouleAssign == text), 
               aes(x=get(sub),y = Scale_prop),colour="red", size=1)+
    ggtitle(paste("Methylation VS ModuleMembership",text,sep='-')) +
    xlab("InModule_Con") + ylab("MethC_Prop")+
    theme(plot.title = element_text(hjust = 0.5)))
}
dev.off()

pdf('PDF_Results_PreS.pdf')
for (i in seq_along(Mod_Index_Pre)){
  text = Mod_Index_Pre[i]
  sub = paste('kME',text,sep = '')
  print(ggplot(datKME,aes(x=get(sub),y = Scale_prop)) + 
          geom_point(colour="grey") +
          geom_point(data = subset(datKME, MdouleAssign == text), 
                     aes(x=get(sub),y = Scale_prop),colour="red", size=1)+
          ggtitle(paste("Methylation VS ModuleMembership",text,sep='-')) +
          xlab("InModule_Con") + ylab("MethC_Prop")+
          theme(plot.title = element_text(hjust = 0.5)))
}
dev.off()

######=========================##########
##        Hyper G test                ##
######========================##########

setwd('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/Net/')
library(readxl)
library(tidyverse)
DiffC2Gene_raw = read_xlsx('DiffC_Gene.xlsx')
dim(DiffC2Gene_raw)
# data pre
DiffC2Gene = DiffC2Gene_raw %>% 
  dplyr::filter(Gene != '-') %>% 
  group_by(Gene) %>% dplyr::count(Region) %>% 
  tidyr::spread(key = Gene, value = n) %>% 
  transpose_df()%>% 
  replace(is.na(.), 0) %>% 
  mutate_at(vars(-Region), as.numeric) %>% 
  rename(Gene = Region)

DiffC2Gene_count = DiffC2Gene %>% 
  mutate(Count1 = `1st_EXON`+ GENE_BODY + INTRON) %>% 
  mutate(SigG1 = ifelse(Count1 > 30, "Sig", "Not")) %>% 
  mutate(Count2 = PROMOTER + TSS + UPSTREAM) %>% 
  mutate(SigG2 = ifelse(Count2 > 5, "Sig", "Not"))

DiffC2Gene_sig = DiffC2Gene_count %>% 
  dplyr::filter(SigG1 == 'Sig' | SigG2 == 'Sig' )

dim(DiffC2Gene_sig)

save(DiffC2Gene,
     DiffC2Gene_count,
     DiffC2Gene_sig,file = 'DiffC2Gene_Sig_Genes.RData')

# get sig diff C gene
Sig_gene_index = unique(DiffC2Gene_sig$Gene)

## get module index 
ref=1; test = 2
Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
Zsummary=Z.PreservationStats$Zsummary.pres
nonpres_index_b = (which(Zsummary < 2))
nonpres_modulenames_b = rownames(Z.PreservationStats)[nonpres_index_b]
Mod_Index_NonPre  = nonpres_modulenames_b[-grep("grey",nonpres_modulenames_b)]
Mod_Index_Pre = rownames(Z.PreservationStats)[-nonpres_index_b]


# get module assign and seperate them
Gene_all = unique(rownames(networkData_normalized))
Gene_net = unique(rownames(networkData_50var_nocrt))
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

######=========================##########
##        Myth_extent porp            ##
######========================##########
# gene pre
setwd('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/Net/')
library(readxl)
library(tidyverse)
# DiffC2Gene_raw = read_xlsx('DiffC_Gene.xlsx')
# DiffC2Gene.extend = DiffC2Gene_raw %>% 
#   dplyr::filter(Gene != '-')
head(DiffC2Gene.extend)
# genome pre
genome <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",  dataset = "btaurus_gene_ensembl", host = "grch37.ensembl.org")
gene = getBM(c("ensembl_gene_id","external_gene_name","description", "start_position", "end_position", "chromosome_name"), mart = genome)
gene_pos_info_bta = dplyr::select(gene,ensembl_gene_id,start_position,end_position,chromosome_name) %>% 
  arrange(ensembl_gene_id) %>% 
  dplyr::mutate_at(vars(chromosome_name),add)

