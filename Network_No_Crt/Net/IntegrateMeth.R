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
#     python rgmatch.py -g Bos_taurus.ARS-UCD1.2.99.gtf -b Diff_C_Sig_BED.bed -r 'gene' -q 4 -o myassociations.txt    ##
## ================================================================================================================== ##
#Data_loci = '/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net'
setwd('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/Net')
#setwd('/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Network_No_Crt/Net/rgmatch')
# Associ_out_raw = read.table('myassoci_exon_5.5k_ext.txt',sep = '\t') %>% data.frame()
library(readxl)
Associ_out_raw = read_xlsx('DiffC_Gene.xlsx') %>% dplyr::filter(Gene != "-") %>% 
  rename(Area=Region)
#%>% dplyr::filter(Region != 'DOWNSTREAM')
#setwd(Data_loci)

### massage diff C associated gene
# cname =as.character(unlist(Associ_out_raw[1,]));attributes(cname) = NULL
# colnames(Associ_out_raw) = cname
# Associ_out_raw = Associ_out_raw[-1,]

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
  #dplyr::select(-Distance,-Transcript,-`Exon/Intron`,-TSSDistance,-PercRegion,-PercArea) %>% 
  group_by(Gene) %>% 
  dplyr::count(Area) %>%
  tidyr::spread(key = Gene, value = n) %>% 
  transpose_df() %>% 
  replace(is.na(.), 0) %>% 
  mutate_at(vars(-Area), as.numeric) %>% 
  dplyr::rename(Gene = Area)

Associ_out_count = Associ_out %>% 
  mutate(Count1 = `1st_EXON`+ GENE_BODY + INTRON) %>% 
  mutate(Count2 = PROMOTER + TSS + UPSTREAM)


# count all Cs using self-wraped function: AGCTcount (getSEQ from biomart mainly)


# load('Genes_C_count_all_Final_api.RData')
# head(Genes_C_count_all_api)
# Genes_C_count_all = Genes_C_count_all_api %>% 
#   rename(Gene = V1, Total_CG = V2,Prom_CG = V3) %>% 
#   mutate_at(vars(Total_CG,Prom_CG),as.numeric) %>% 
#   replace(is.na(.), 0)

load('Total_C_count_raw.rda')
Total_C_count = Total_C_count_raw %>% 
  mutate(Count1_all = `1st_EXON`+ GENE_BODY + INTRON) %>% 
  mutate(Count2_all = PROMOTER + TSS + UPSTREAM)
Total_C_count_2join = Total_C_count %>% 
  dplyr::select(Gene,Count1_all,Count2_all)

"/" <- function(x,y) ifelse(y==0,0,base:::"/"(x,y)) # special division
Associ_out_count_final = Associ_out_count %>% 
  dplyr::left_join(Total_C_count_2join,by = c('Gene'= 'Gene')) %>% 
  dplyr::select(Gene,Count1,Count2,Count1_all,Count2_all) %>% 
  mutate(Prop_Body = Count1/Count1_all) %>% 
  mutate(Prop_Prpt = Count2/Count2_all) %>% 
  mutate(Prop_All = (Count1+Count2)/(Count1_all+Count2_all))


# Count1 = 1st_EXON+ GENE_BODY + INTRON
# Count2 = PROMOTER + TSS + UPSTREAM
Genes_meth_prop = Associ_out_count_final %>% 
  dplyr::filter(!is.na(Prop_Body) & Prop_Body <=1) %>% 
  dplyr::filter(!is.na(Prop_Prpt) & Prop_Prpt <=1) %>% 
  dplyr::filter(!is.na(Prop_All) & Prop_All <=1)

sort(Genes_meth_prop$Prop_Prpt,decreasing = T)

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

#save(Genes_meth_prop,file = 'Genes_meth_prop.txt')
save(Genes_meth_prop,
     Genes_meth_select,
     Diff_Meth_Gene_index,
     file = 'Genes_meth_prop.rda')

# Gather Info: KME and Meth
load('Genes_meth_prop.rda')
datKME_tmp = signedKME(datExpr_control, MEs_control)
'/' <- base:::"/"
datKME = datKME_tmp %>% 
  dplyr::mutate(Gene = rownames(datKME_tmp)) %>% 
  dplyr::mutate(MdouleAssign = moduleColors_control) %>% 
  dplyr::left_join(Genes_meth_prop, by= c("Gene" = "Gene"))


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

names(datKME)

plot(test)
test = scale(data_tmp$Scale_prop_all)
summary(test)
names(datKME)

#i = 3
pdf('PDF_Results_NonP.pdf')
for (i in seq_along(Mod_Index_NonPre)){
  text = Mod_Index_NonPre[i]
  sub = paste('kME',text,sep = '')
  # subset input data
  data_tmp = subset(datKME,MdouleAssign == text) %>% 
    dplyr::mutate(Scale_KME = scale(get(sub))) %>%
    dplyr::mutate(Scale_prop_all = scale(Diff_Prop_all)) %>% 
    dplyr::mutate(Scale_prop_prmt = scale(Diff_Prop_prom)) %>% 
    dplyr::mutate(Scale_prop_body = scale(Diff_Prop_body)) %>% 
    # dplyr::mutate(Scale_prop_all = Diff_Prop_all/max(Diff_Prop_all)) %>% 
    # dplyr::mutate(Scale_prop_prmt = Diff_Prop_prom/max(Diff_Prop_prom)) %>% 
    # dplyr::mutate(Scale_prop_body = Diff_Prop_body/max(Diff_Prop_body)) %>% 
    dplyr::select(Gene,Scale_KME,Scale_prop_all,Scale_prop_prmt,Scale_prop_body)
  # pivot rotate
  data_module = data_tmp %>%
    pivot_longer(cols = c(Scale_prop_all,Scale_prop_prmt,Scale_prop_body),
                 names_to = "Cat", values_to = "Prop")
  print(
    ggplot(data_module, aes(x = Scale_KME ,y = Prop, colour = Cat)) + # y = Scale_prop
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


pdf('PDF_Results_Pre.pdf')
for (i in seq_along(Mod_Index_Pre)){
  text = Mod_Index_Pre[i]
  sub = paste('kME',text,sep = '')
  # subset input data
  data_tmp = subset(datKME,MdouleAssign == text) %>% 
    dplyr::mutate(Scale_KME = scale(get(sub))) %>%
    dplyr::mutate(Scale_prop_all = scale(Diff_Prop_all)) %>% 
    dplyr::mutate(Scale_prop_prmt = scale(Diff_Prop_prom)) %>% 
    dplyr::mutate(Scale_prop_body = scale(Diff_Prop_body)) %>% 
    # dplyr::mutate(Scale_prop_all = Diff_Prop_all/max(Diff_Prop_all)) %>% 
    # dplyr::mutate(Scale_prop_prmt = Diff_Prop_prom/max(Diff_Prop_prom)) %>% 
    # dplyr::mutate(Scale_prop_body = Diff_Prop_body/max(Diff_Prop_body)) %>% 
    dplyr::select(Gene,Scale_KME,Scale_prop_all,Scale_prop_prmt,Scale_prop_body)
  # pivot rotate
  data_module = data_tmp %>%
    pivot_longer(cols = c(Scale_prop_all,Scale_prop_prmt,Scale_prop_body),
                 names_to = "Cat", values_to = "Prop")
  print(
    ggplot(data_module, aes(x = Scale_KME ,y = Prop, colour = Cat)) + # y = Scale_prop
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

names(Preserved_Gene_list)
# genes captured in all category
Gene_all = unique(rownames(networkData_normalized))
Gene_net = unique(rownames(networkData_50var_nocrt))

Gene_grey = datKME %>% 
  dplyr::filter(MdouleAssign == 'grey') %>% 
  dplyr::select(Gene) %>% 
  unlist(use.names = F)

table(Diff_Meth_Gene_index %in% Gene_all)
table(Diff_Meth_Gene_index %in% Gene_net)
table(Diff_Meth_Gene_index %in% Gene_grey)

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


# devtools::install_github("slowkow/tftargets")
library(tftargets)
length(TRED)#etz
length(ENCODE)#etz
length(TRRUST)
length(ITFP)
length(Neph2012)#etz
length(Marbach2016) 


#
tf_database_index  = c('TRED','ENCODE','TRRUST','ITFP','Neph2012','Marbach2016')
# 


# get all genes _ name index
load('gene_symbols_genome.rda')
head(gene_symbols_genome)

Trun2Bt = function(list,species = 'Bt'){
  library(tidyverse)
  list_new = list()
  for (i in seq_along(names(list))){
    tmp_name = names(list)[i]
    tmp_index = list[[i]] %>% unlist(use.names = F)
    for(p in seq_along(tmp_index)){
      sub_tmp = alias2Symbol(tmp_index[p],species = species, expand.symbols = F)
      tmp_index[p] = ifelse(length(sub_tmp)==0,tmp_index[p],sub_tmp)
    }
    list_new[[i]] = tmp_index
    nwe_name_tmp = alias2Symbol(tmp_name,species = species,expand.symbols = F)
    names(list_new)[i] = ifelse(length(nwe_name_tmp)==0,tmp_name,nwe_name_tmp)
  }
  return(list_new)
}

ITFP_test = ITFP[1:3]

ITFP_new = Trun2Bt(ITFP_test)



# container: given a module(gene set); search for overlap for each TF in all databases, record the overlap
library(tidyverse)
OUT = data.frame(DataBase = c(),
                 TF_Name = c(),
                 OverNum_Entrz = c(),
                 OverNum_Ensl = c(),
                 OverGene = c(),
                 Module = c())
for (p in seq_along(UnPreserved_Gene_list)){
  tmp = (UnPreserved_Gene_list)[[p]]
  module.name = names(UnPreserved_Gene_list)[p]
  gene_collection_tmp = gene_symbols_genome %>% 
    dplyr::filter(ensembl_gene_id %in% tmp) %>% 
    rename(Ens = ensembl_gene_id, 
           Symbol = external_gene_name,
           Entrez = ENTREZID)
  out = data.frame(DataBase = c(),
                   TF_Name = c(),
                   OverNum_Entrz = c(),
                   OverNum_Ensl = c(),
                   OverGene = c(),
                   Module = c())
  Module.Gene.Input = gene_collection_tmp
  for (i in seq_along(tf_database_index)){
    out_tmp = data.frame(DataBase = c(),
                         TF_Name = c(),
                         OverNum_Entrz = c(),
                         OverNum_Ensl = c(),
                         OverGene = c(),
                         Module = c())
    DB_name = tf_database_index[i]
    tmp_TF_index = names(get(DB_name))
    for (j in seq_along(tmp_TF_index)){
      tmp_TF_gene_list = get(DB_name)[[j]]
      out_tmp[i+j-1,1] = DB_name
      out_tmp[i+j-1,2] = tmp_TF_index[j]
      module.gene.entrz = Module.Gene.Input[,3]
      module.gene.ensl = Module.Gene.Input[,2]
      Ov1_Entrz = intersect(module.gene.entrz,tmp_TF_gene_list)
      Ov2_Ensl  = intersect(module.gene.ensl,tmp_TF_gene_list)
      out_tmp[i+j-1,3] = length(Ov1_Entrz) 
      out_tmp[i+j-1,4] = length(Ov2_Ensl)
      findG1_Entrz = paste(Ov1_Entrz,collapse = '/')
      findG2_Ensl = paste(Ov2_Ensl,collapse = '/')
      out_tmp[i+j-1,5] = ifelse(length(Ov1_Entrz) == 0,findG2_Ensl,findG1_Entrz)
      out_tmp[i+j-1,6] = module.name
    }
    out = rbind(out,out_tmp)
  }
  OUT = rbind(OUT,out)
}
#
ModuleSize = data.frame(Module = names(UnPreserved_Gene_list),
                        Size = sapply(UnPreserved_Gene_list, length))
OUT_final = OUT %>%
  dplyr::filter(!(V3 == 0) | !(V4 == 0)) %>% 
  drop_na() %>% 
  rename(Database = V1,
         TF_Name = V2,
         Overlap_Entrz = V3,
         Overlep_Ensl = V4,
         FindG = V5,
         Module = V6) %>% 
  group_by(Module) %>% 
  left_join(ModuleSize, by = c('Module' = 'Module')) %>% 
  dplyr::filter(Overlep_Ensl >= 0.4 * Size)

OUT_Final_subname = OUT_final %>% 
  mutate(subname = TF_Name)
for (i in seq_len(dim(OUT_Final_subname)[1])){
  tmp_name = OUT_Final_subname[i,2]
  sub = alias2Symbol(tmp_name,species = "Bt", expand.symbols = F)
  re = ifelse(length(sub) == 0,tmp_name,sub)
  OUT_Final_subname[i,9] = re
}

OUT_Final = OUT_Final_subname %>% dplyr::select(-FindG) %>% dplyr::select(-subname) %>% 
  rename(Sub_name = V9) %>% 
  left_join(gene_symbols_genome, by = c('Sub_name' = 'external_gene_name')) %>% 
  filter(!(Database == "Marbach2016"))
OUT_Final[OUT_Final$Sub_name=='VARS',8] = 'ENSBTAG00000005631'
OUT_Final[OUT_Final$Sub_name=='VARS',9] = '505556'
#
view(OUT_Final)

names(gene_symbols_genome)
filter(gene_symbols_genome,external_gene_name == 'VARS1')
filter(gene_symbols_genome,ensembl_gene_id == 'ENSBTAG00000054322')
filter(gene_symbols_genome,ENTREZID == '518045')




