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
load('MethEval_all.RData')

# Gather Info: KME and Meth
datKME_tmp = signedKME(datExpr_control, MEs_control)
datKME = datKME_tmp %>% 
  dplyr::mutate(ensembl_gene_id = rownames(datKME_tmp)) %>% 
  dplyr::mutate(MdouleAssign = moduleColors_control) %>% 
  dplyr::left_join(Meth_prmt, by= c("ensembl_gene_id" = "ensembl_gene_id"))
  #%>% dplyr::filter(meth != 1)
head(datKME)
summary(datKME$meth)
table(datKME$ensembl_gene_id%in% Meth_prmt$ensembl_gene_id)

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

######=========================##########
##        diff C prop vs M.M           ##
######=========================##########
# 
pdf('PDF_Results.pdf')
for (i in seq_along(Mod_Index_NonPre)){
  text = Mod_Index_NonPre[i]
  sub = paste('kME',text,sep = '')
  print(ggplot(datKME,aes(x=get(sub),y=meth)) + 
    geom_point(colour="grey") +
    geom_point(data = subset(datKME,MdouleAssign == text), 
               aes(x=get(sub),y=meth),colour="red", size=1)+
    ggtitle(paste("Methylation VS ModuleMembership",text,sep='-')) +
    xlab("InModule_Con") + ylab("MethC_Prop")+
    theme(plot.title = element_text(hjust = 0.5)))
}
dev.off()

######=========================##########
##        Hyper G test                ##
######========================##########
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
  transpose_df() %>% 
  replace(is.na(.), 0) %>% 
  mutate_at(vars(-Region), as.numeric) %>% 
  rename(Gene = Region)

DiffC2Gene_count = DiffC2Gene %>% 
  mutate(Count1 = `1st_EXON`+ GENE_BODY + INTRON) %>% 
  mutate(SigG1 = ifelse(Count1 > 40, "Sig", "Not")) %>% 
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


table(Sig_gene_index%in%Gene_all)
table(Sig_gene_index%in%Gene_net)

getOption("max.print")
test = listAttributes(genome)
test[334:dim(test)[1],]
test[666:dim(test)[1],]
test[998:dim(test)[1],]

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
genome <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",  dataset = "btaurus_gene_ensembl", host = "grch37.ensembl.org")
gene = getBM(c("ensembl_gene_id","external_gene_name","description", "start_position", "end_position", "chromosome_name"), mart = genome)
gene_pos_info_bta = dplyr::select(gene,ensembl_gene_id,start_position,end_position,chromosome_name) %>% 
  arrange(ensembl_gene_id) %>% 
  dplyr::mutate_at(vars(chromosome_name),add)

# function pre
library(biomaRt)
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
                  substring(seq_p2,nchar(seq_p2)-downstream + 1,nchar(seq_p2)),
                  sep = '')
  nchar(seq_all)
  Find_loc = data.frame(str_locate_all(seq_all,find))
  total_c = nrow(Find_loc)
  return(total_c)
}

Genes_C_count_all = data.frame(Gene = c(),
                               total = c())
for (i in seq_along(Gene_all)){
  Genes_C_count_all[i,1] = Gene_all[i]
  Genes_C_count_all[i,2] = AGCTcount(Gene_all[i],genome = genome)
  message('Working on ',Gene_all[i])
}

#
Gene_all = unique(rownames(networkData_normalized))
Gene_net = unique(rownames(networkData_50var_nocrt)) 
#
Gene_DiffC_index = unique(DiffC2Gene.extend$Gene)

#length(unique(Gene_DiffC_index))
#table(Gene_DiffC_index %in% Gene_all)

### massage diff C associated gene
Associ_out_raw = read.table('myassociations2.txt',sep = '\t') %>% data.frame()
cname =as.character(unlist(Associ_out_raw[1,]));attributes(cname) = NULL
colnames(Associ_out_raw) = cname
Associ_out_raw= Associ_out_raw[-1,]
####
length(table(Associ_out_raw$Gene))
table(Associ_out_raw$Area)

Associ_within 

Associ_up


