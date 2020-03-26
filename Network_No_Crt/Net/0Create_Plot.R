#install.packages('customLayout')
library(ggplot2)
library(ggplotify)
library(cowplot)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)

sizeGrWindow(12,36)
pdf('test.pdf')
#par(mfrow=c(6,1),mar = c(1,1,1,1))

par(mar = c(2,4.1,4.1,2))
plot(METree_control, 
     main = "Clustering of Initial MEs (Control Diet)",
     xlab = '', sub = "", ylab = "Height",
     hang = -1,cex = 0.35);MEDissThres = 0.2;abline(h = MEDissThres, col = "red")

par(mar = c(2,4.4,4.1,2))
plotDendroAndColors(geneTree_control, 
                    cbind(dynamicColors_control,moduleColors_control), 
                    c("Dynamic Tree Cut", "Merged Modules"), 
                    dendroLabels=F, hang=0.03, addGuide=TRUE,
                    guideHang=0.05,
                    main="Gene Cluster Dendrogram (Control Diet)")

plotDendroAndColors(geneTree_treatment,
                    moduleColors_control,
                    "Modules", 
                    dendroLabels=F,hang=0.03, addGuide=TRUE,
                    guideHang=0.05,setLayout = F,
                    main="Gene Cluster Dendrogram (Methionine Diet)")
dev.off()



######################  Fig1 Net-Dengdro #############################

DengroControl = as.ggplot(~plotDendroAndColors(geneTree_control, 
                                     cbind(dynamicColors_control,moduleColors_control), 
                                     c("Dynamic Tree Cut", "Merged Modules"), 
                                     dendroLabels=F,
                                     hang = 0.03,addGuide=TRUE,
                                     guideHang=0.05, cex.colorLabels = 1,
                                     main="Gene Cluster Dendrogram (Control Diet)"))
par(mar = c(2,4.4,4.1,2))
DengroTreatment = as.ggplot(~plotDendroAndColors(geneTree_treatment,moduleColors_control,
                                                  "Modules", 
                                                  dendroLabels=F,hang=0.03, addGuide=TRUE,
                                                  guideHang=0.05, cex.colorLabels = 1,
                                                  main="Gene Cluster Dendrogram (Methionine Diet)"))


library(ggpubr)
tiff("Fig1-Net-Dengdro.tiff", width = 12, height = 20, units = 'in', res = 300)
#grid.arrange(DengroControl,DengroTreatment,ncol=1,nrow =2)
ggarrange(DengroControl,DengroTreatment, labels = c("A", "B"),ncol=1,nrow =2)
# plot_grid(DengroControl,DengroTreatment,align = c("v"),
#           labels = c("A","B"), label_size= 20,label_colour = "black")
dev.off()


######################  Fig 2 Net-Presv-stats #####################
library(ggrepel)
# Preservation Z summary
PresZsum =
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
# Preservation MedianR
PresMedianR = 
ggplot(MP_Stats_nobig, aes(moduleSize,medianRank.pres,label = Row.names)) +
  geom_point(size = 5,aes(colour = Row.names))+
  geom_text_repel(vjust = -2)+
  ggtitle("Preservation Median Rank") +
  xlab("ModuleSize") + ylab("Preservation Median Rank") +
  theme(plot.title = element_text(hjust = 0.5,size=16, face="bold.italic"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))


tiff("Fig2-Net-Presv-stats.tiff", width = 16, height = 10, units = 'in', res = 300)
#PresZsum / PresMedianR
ggarrange(PresZsum,PresMedianR,labels = c("A", "B"),ncol=2,nrow =1)
# plot_grid(PresZsum,PresMedianR,
#           align = c("h"),
#           labels = c("A","B"), label_size= 20,label_colour = "black")
dev.off()


####################  Fig3 Functional Enrich ###################
library(openxlsx);library(tidyverse)
require(openxlsx)
setwd('/Users/liulihe95/Desktop/Methionine/Network_No_Crt/enrich_results')
emptyentry = data.frame(ID = c(),
                        Name = c(),
                        Total_Genes = c(),
                        Significant_Genes = c(),
                        pvalue_r = c(),hitsPerc = c(),Database = c())
enrich_terms_out <- list()
for (i in seq_along(ModuleName_Unpreserved)){
  module = ModuleName_Unpreserved[i]
  target_sheet = i;pthres = 0.001
  #
  module_go_raw = read.xlsx('GO_Results_all_0113.xlsx',sheet = target_sheet)
  if (is.null(module_go_raw)) {
    module_go = emptyentry
  } else {module_go = module_go_raw %>% dplyr::select(-findG,-namespace_1003) %>% 
    dplyr::filter(pvalue<=pthres) %>% 
    rename(pvalue_r = pvalue,ID = go_id,Name = GO_Name) %>% 
    mutate(Database = 'GO')}
  #
  module_interpro_raw = read.xlsx('Interpro_Results_all_0113.xlsx',sheet = target_sheet)
  if (is.null(module_interpro_raw)){
    module_interpro = emptyentry
  } else {
    module_interpro = module_interpro_raw %>% dplyr::select(-findG) %>% 
                               dplyr::filter(pvalue_r<=pthres) %>% 
                               rename(ID = InterproID,Name = Interpro_Name) %>% 
                               mutate(Database = 'Interpro')
  }
  #
  module_kegg_raw = read.xlsx('KEGG_Results_all_0113.xlsx',sheet = target_sheet) 
  if (is.null(module_kegg_raw)){
    module_kegg = emptyentry
  } else {
    module_kegg = module_kegg_raw %>% dplyr::select(-findG) %>% 
      dplyr::filter(pvalue_r<=pthres) %>% 
      rename(ID = KEGGID,Name = KEGGTERM) %>% 
      mutate(Database = 'KEGG')
  }
  #
  module_mesh_raw = read.xlsx('Mesh_Results_all_0115.xlsx',sheet = target_sheet)
  if (is.null(module_mesh_raw)){
    module_mesh = emptyentry
  } else {
    module_mesh = module_mesh_raw %>% dplyr::select(-findG) %>% 
      dplyr::filter(pvalue_r<=pthres) %>% 
      rename(ID = MeshID,Name = MeshTerm) %>% 
      mutate(Database = 'MeSH')
  }
  #
  module_reactome_raw = read.xlsx('Reactome_all_path_Results_all_0113.xlsx',sheet = target_sheet)
  if (is.null(module_reactome_raw)){
    module_reactome = emptyentry
  } else {
    module_reactome = module_reactome_raw %>% dplyr::select(-findG) %>%
      dplyr::filter(pvalue_r<=pthres) %>% 
      rename(ID = ReactomeID,Name = ReactomeName) %>% 
      mutate(Database = 'reactome')
  }
  #
  module_msig_raw = read.xlsx('Msig_Results_all_0124.xlsx',sheet = target_sheet)
  if (is.null(module_msig_raw)){
    module_msig = emptyentry
  } else {
    module_msig = module_msig_raw %>% dplyr::select(-findG) %>% 
      dplyr::filter(pvalue_r<=pthres) %>% 
      rename(ID = MsigID,Name = MsigTerm) %>% 
      mutate(Database = 'Msig')
  }
  #
  module_summary = rbind(module_go,module_interpro,module_kegg,
                         module_mesh,module_reactome,module_msig)
  enrich_terms_out[[i]] = module_summary
  names(enrich_terms_out)[i] = module
  write.xlsx(enrich_terms_out, file = 'Enrich_results_summary.xlsx')
}

#
item_index = c(1,4,30,56,70,212,461,8,9,11,48,235,28,74,18,65,47,50)
item_cate_func = c(rep('ribosomal',7),rep('Mitochondrial',5),rep('ATP',2),rep('rRNA',2),rep('NAD(P)H',2))
length(item_index)

antiquewhile_enrichsummary_plot = 
  read.xlsx('Enrich_results_summary.xlsx',sheet = 1) %>% 
  slice(item_index) %>% mutate(Func_cate = item_cate_func) %>% 
  mutate(Ylable = paste0('[',ID,']',':',Name,'(',Total_Genes,')')) %>% 
  mutate(ID_Y = paste0('[',ID,']',' (',Total_Genes,')')) %>% 
  mutate(log10pvalue = -log10(pvalue_r)) %>% 
  mutate_at('Func_cate',fct_infreq)
antiquewhile_enrichsummary_plot$ID =
  factor(antiquewhile_enrichsummary_plot$ID,
         levels = c('R-BTA-72706','R-BTA-72702','M17739','M14381','IPR011332','GO:0022627','GO:0003735',
                    'R-BTA-6791226','GO:0019843',
                    'R-BTA-163200','GO:0015986',
                    'M11099','GO:0032981','GO:0005753','GO:0005743','D025261',
                    'D016660','D009247'))
#
ggEnrich = ggplot()+
  geom_bar(antiquewhile_enrichsummary_plot,
           mapping = aes(x = ID_Y,y = hitsPerc,fill = Func_cate),
           stat = "identity", width=0.2) +
  scale_y_continuous(name = expression(bold("Hits Percentage")),
                     sec.axis = sec_axis(~.* 30/100, name = expression(bold("-log10(pvalue)"))),
                     limits = c(0,100), breaks = seq(0,100,by=5))+
  geom_point(antiquewhile_enrichsummary_plot,
             mapping = aes(x = ID_Y,y = log10pvalue*100/30,fill = Func_cate)) +
  theme(legend.position="none",
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=8,face="bold",color = "black"),
        axis.text.x = element_text(size=10,face="bold",color = "black")) + 
  #facet_wrap(~Func_cate,scales = 'free',dir = "v")+
  coord_flip()

dev.off()
 
tiff("Fig3-Enrichment-Bar.tiff", width = 10, height = 6, units = 'in', res = 300)
plot(ggEnrich)
dev.off()

####################  Fig4 Methylation Prop ###################

#all -c-associ-gene 25491(24344)
# 1147 all zeros in any region

# sig -c-associ-gene 10247

# gene-rna = 20479 # no two low
# gene-net (50%) = 7035

length(unique(Associ_out_count_final$Gene))

names(Total_C_count_raw)

test = Total_C_count_raw %>% dplyr::select(-DOWNSTREAM) %>% 
  mutate(rowsum =rowSums(.[2:7])) %>% 
  dplyr::filter(!(rowsum == 0))


length(which(test$rowsum ==0))


sigflaggene = (unique(Associ_out_raw$Gene))
allflaggene = (unique(Total_C_count_raw$Gene))

Gene_all = unique(rownames(networkData_normalized))
Gene_net = unique(rownames(networkData_50var_nocrt))

length(Gene_net)

table(Gene_all %in% allflaggene)

table(Gene_net %in% allflaggene)


#
save(UnPreserved_Gene_list,
     Preserved_Gene_list,
     file = 'Gene_list_by_module.rda')

ModuleSize = data.frame(Module = names(Gene_list_all),
                        Size = sapply(Gene_list_all, length))
Overal_match_color = data.frame(Gene = colnames(datExpr_control),
                                Module = moduleColors_control) %>% 
  mutate(Cate = ifelse(Module %in% ModuleName_Unpreserved,'Unp','Pre'))

#Bta_TF_list;Bta_TcoF_list

Gene_Meth_Viol = Meth_Prop_Univ %>% 
  left_join(Overal_match_color,by = c('Gene'='Gene')) %>% drop_na() %>% 
  mutate(TFinfo = ifelse(Gene %in% Bta_TF_list,'TF',
                         ifelse(Gene %in% Bta_TcoF_list,'TcoF','NOT')))

table(Gene_Meth_Viol$Cate)


#
Gene_Meth_Viol2 = Gene_Meth_Viol
Gene_Meth_Viol2$Cate = 'All'
#
rbind1 = rbind(Gene_Meth_Viol,Gene_Meth_Viol2) %>% 
  dplyr::select(Gene,Prop_Prpt,Prop_All,Module,Cate)
rbind1_tmp = rbind1;rbind1_tmp$Cate = 'All'
rbind1_final = rbind(rbind1,rbind1_tmp) %>% 
  dplyr::select(-Prop_Prpt) %>% 
  mutate(Region = 'whole') %>% 
  rename(Prop = Prop_All)
#
rbind2 = rbind(Gene_Meth_Viol,Gene_Meth_Viol2) %>% 
  dplyr::select(Gene,Prop_Prpt,Prop_All,Module,Cate)
rbind2_tmp = rbind2;rbind2_tmp$Cate = 'All'
rbind2_final = rbind(rbind2,rbind2_tmp) %>%
  dplyr::select(-Prop_All) %>% 
  mutate(Region = 'promoter') %>% 
  rename(Prop = Prop_Prpt)

rbind_final = rbind(rbind1_final,rbind2_final) %>% 
  mutate_at('Region',as.factor) %>% 
  mutate_at('Cate',as.factor)
str(rbind_final)
#labels <- labeller(`promter` = expression('Methprop'['REG']), `whole` = expression('Methprop'['GENE']))

levels(rbind_final$Region)[levels(rbind_final$Region)=="promoter"] <- 'Methprop_REG'
levels(rbind_final$Region)[levels(rbind_final$Region)=="whole"] <-  'Methprop_GENE'
levels(rbind_final$Cate)[levels(rbind_final$Cate)=="All"] <-  'All(6735)'
levels(rbind_final$Cate)[levels(rbind_final$Cate)=="Pre"] <-  'Preserved(5589)'
levels(rbind_final$Cate)[levels(rbind_final$Cate)=="Unp"] <-  'Unpreserved(1146)'

P1 = ggplot() +
  geom_violin(data = rbind_final, aes(x = Cate,y = Prop,fill=Cate))+
  #geom_violin(alpha =.8,width = .8) +
  #geom_boxplot(outlier.alpha = 0.2,outlier.size = 0.3) +
  xlab("Gene Preservation Status") + ylab('') +
  theme(legend.position="none",
        axis.title.x = element_text(size=12, face="bold"),
        axis.text.x = element_text(size=9,color ='black',face="bold"),
        axis.text.y = element_text(size=9,color ='black',face="bold"),
        strip.text = element_text(size=12,color ='black',face="bold")) +
  facet_grid(~Region)
P1

####################  Fig5 Methylation Prop ###################
#
table(Gene_Meth_Viol_TF$Cate)

Gene_Meth_Viol_TF = Gene_Meth_Viol %>% 
  dplyr::filter(TFinfo != 'NOT')

Gene_Meth_Viol2 = Gene_Meth_Viol_TF
Gene_Meth_Viol2$Cate = 'All'
#
rbind1 = rbind(Gene_Meth_Viol_TF,Gene_Meth_Viol2) %>% 
  dplyr::select(Gene,Prop_Prpt,Prop_All,Module,Cate)
rbind1_tmp = rbind1;rbind1_tmp$Cate = 'All'
rbind1_final = rbind(rbind1,rbind1_tmp) %>% 
  dplyr::select(-Prop_Prpt) %>% 
  mutate(Region = 'whole') %>% 
  rename(Prop = Prop_All)
#
rbind2 = rbind(Gene_Meth_Viol_TF,Gene_Meth_Viol2) %>% 
  dplyr::select(Gene,Prop_Prpt,Prop_All,Module,Cate)
rbind2_tmp = rbind2;rbind2_tmp$Cate = 'All'
rbind2_final = rbind(rbind2,rbind2_tmp) %>%
  dplyr::select(-Prop_All) %>% 
  mutate(Region = 'promoter') %>% 
  rename(Prop = Prop_Prpt)

rbind_final = rbind(rbind1_final,rbind2_final) %>% 
  mutate_at('Region',as.factor) %>% 
  mutate_at('Cate',as.factor)
str(rbind_final)

levels(rbind_final$Region)[levels(rbind_final$Region)=="promoter"] <- 'Methprop_REG'
levels(rbind_final$Region)[levels(rbind_final$Region)=="whole"] <-  'Methprop_GENE'
levels(rbind_final$Cate)[levels(rbind_final$Cate)=="All"] <-  'All(292)'
levels(rbind_final$Cate)[levels(rbind_final$Cate)=="Pre"] <-  'Preserved(253)'
levels(rbind_final$Cate)[levels(rbind_final$Cate)=="Unp"] <-  'Unpreserved(39)'

P2 = ggplot() +
  geom_violin(data = rbind_final, aes(x = Cate,y = Prop,fill=Cate))+
  #geom_violin(alpha =.8,width = .8) +
  #geom_boxplot(outlier.alpha = 0.2,outlier.size = 0.3) +
  xlab("TF/TcoF Detected") + ylab('') +
  theme(legend.position="none",
        axis.title.x = element_text(size=12, face="bold"),
        axis.text.x = element_text(size=9,color ='black',face="bold"),
        axis.text.y = element_text(size=9,color ='black',face="bold"),
        strip.text = element_text(size=12,color ='black',face="bold")) +
  facet_grid(~Region)
P2

library(gridExtra)

tiff("Fig4-Gene-TF-Methy-by-Cate.tiff",
     width = 14, height = 8, units = 'in', res = 300)
ggarrange(P1,P2,labels = c("A", "B"),ncol=1,nrow =2)
#grid.arrange(P1,P2,ncol=1,nrow =2)
# plot_grid(P1,P2,align = c("v"),
#            labels = c("A","B"), label_size= 12,label_colour = "black")
dev.off()









