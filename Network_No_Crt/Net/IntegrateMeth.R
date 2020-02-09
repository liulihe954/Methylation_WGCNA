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
library(readxl)
library(tidyverse)
DiffC2Gene_raw = read_xlsx('DiffC_Gene.xlsx')
DiffC2Gene = DiffC2Gene_raw %>% 
  dplyr::filter(Gene != '-') %>% 
  group_by(Gene) %>% 
  count(Region)

selected_gene = c()
gene_index = unique(DiffC2Gene$Gene)
for (i in seq_along(gene_index)){
  anchor = gene_index[i]
  tmp %>% dplyr::filter(Gene == anchor)
}

names(DiffC2Gene)

head(DiffC2Gene,500)
table(DiffC2Gene_raw$Region)

names(DiffC2Gene)

mtcars %>% group_by(cyl) %>% tally()
mtcars %>% group_by(gear) %>% count(carb)


Con1 <- "Region %in% c('1st_EXON','INTRON','GENE_BODY')"
Con2 <- "Region %in% c('PROMOTER','UPSTREAM','TSS')"
DiffC2Gene1 = DiffC2Gene_raw %>% dplyr::filter(Region %in% c('1st_EXON','INTRON','GENE_BODY'))
DiffC2Gene2 = DiffC2Gene_raw %>% dplyr::filter(Region %in% c('PROMOTER','UPSTREAM','TSS'))





# DiffC2Gene = DiffC2Gene_raw %>% dplyr::filter(Gene != '-')
# DiffC2Gene_p1 = dplyr::filter(`Meth Change %` >= 0.2*max(abs(`Meth Change %`))｜`q-value` <= 0.1)
# DiffC2Gene_p2 = dplyr::filter(`Meth Change %` >= 0.2*max(abs(`Meth Change %`))｜`q-value` <= 0.1)
# Region %in% c('1st_EXON','PROMOTER','TSS'))
length(unique(DiffC2Gene$Gene))
table(DiffC2Gene_raw$Region)
head(DiffC2Gene)












