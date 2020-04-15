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
dev.off()



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
tiff("Fig1-Net-Dengdro.tiff", width = 12, height = 20, units = 'in', res = 500)
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
  #ggtitle("Preservation Zsummary") +
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
  #ggtitle("Preservation Median Rank") +
  xlab("ModuleSize") + ylab("Preservation Median Rank") +
  theme(plot.title = element_text(hjust = 0.5,size=16, face="bold.italic"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))

library(gridExtra)
tiff("Fig2-Net-Presv-stats.tiff", width = 16, height = 10, units = 'in', res = 500)
#PresZsum / PresMedianR
# grid.arrange(PresZsum,PresMedianR,
#              top = "Title",
#              layout_matrix = matrix(c(1,2), ncol=2, byrow=TRUE))
figure = ggarrange(PresZsum,PresMedianR,
          labels = c("A", "B"),ncol=2,nrow =1)
annotate_figure(figure,
                top = text_grob("Module Preservation Statistics", color = "black", face = "bold", size = 20))
# plot_grid(PresZsum,PresMedianR,
#           align = c("h"),
#           labels = c("A","B"), label_size= 20,label_colour = "black")
dev.off()

####################  Fig3 Functional Enrich ###################
library(openxlsx);library(tidyverse)
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
library(gplots)
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
  mutate_at('Func_cate',fct_infreq) %>% 
  mutate(fill_plot = as.character(recode(Func_cate,
                            ribosomal = col2hex('red'),
                            Mitochondrial = col2hex('blue'),
                            ATP = col2hex('green'),
                            `NAD(P)H` = col2hex('orange'),
                            rRNA = col2hex('purple'))))
antiquewhile_enrichsummary_plot$ID_Y =
  factor(antiquewhile_enrichsummary_plot$ID_Y,
         levels = rev(c('[R-BTA-72706] (66)','[R-BTA-72702] (42)','[M17739] (585)','[M14381] (310)','[IPR011332] (8)','[GO:0022627] (45)','[GO:0003735] (157)',
                    '[R-BTA-6791226] (99)','[GO:0019843] (30)',
                    '[R-BTA-163200] (52)','[GO:0015986] (17)',
                    '[M11099] (68)','[GO:0032981] (34)','[GO:0005753] (18)','[GO:0005743] (224)','[D025261] (19)',
                    '[D016660] (39)','[D009247] (38)')))
colorindex = as.character(antiquewhile_enrichsummary_plot$fill_plot)

# 
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# 
# n = 18
# cols = gg_color_hue(n)
# 
# tst = colourName(cols)
# ggplotColours <- function(n = 6, h = c(0, 360) + 15){
#   if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
#   hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
# }
# ggplotColours(n = 3)
# colourName(ggplotColours(n = 3))
# 
# color = data.frame(ggplot_build(ggEnrich)$data)
# table(color$fill)

#fill_plot = as.vector(antiquewhile_enrichsummary_plot$fill_plot)
ggEnrich = ggplot()+
  geom_bar(antiquewhile_enrichsummary_plot,
           mapping = aes(x = ID_Y, y = hitsPerc),
           fill = colorindex, ## colour = fill_plot,fill = rep('#619CFF',18)
           stat = "identity", width=0.2) +
  #scale_colour_manual(values= colorindex) +
  #scale_fill_manual(values = fill_plot)+
  scale_y_continuous(name = expression(bold("Percentage of Significant Genes")),
                     sec.axis = sec_axis(~.* 30/100, name = expression(bold("-log10(P-value)"))),
                     limits = c(0,100), breaks = seq(0,100,by=5))+
  geom_point(antiquewhile_enrichsummary_plot,
             mapping = aes(x = ID_Y,y = log10pvalue*100/30,fill = Func_cate)) +
  theme(legend.position="none",
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=8,face="bold",color = "black"),
        axis.text.x = element_text(size=10,face="bold",color = "black")) + 
  coord_flip()
ggEnrich


#F8766D #00BFC4
tiff("Fig3-Enrichment-Bar.tiff", width = 10, height = 6, units = 'in', res = 500)
print(ggEnrich)
dev.off()

####################  Fig4 Methylation Prop ###################

#all -c-associ-gene 25491(24344)
# 1147 all zeros in any region

# sig -c-associ-gene 10247

# gene-rna = 20479 # no two low
# gene-net (50%) = 7035

alltf=c(Bta_TF_list,Bta_TcoF_list)



table(alltf %in% Gene_net)

table(alltf %in% Gene_list_Pre)
table(alltf %in% Gene_list_Unpre)
table(alltf %in% Gene_grey)





test = Associ_out_all %>% dplyr::select(Gene) %>% unlist(use.names = F)

lost = Gene_net[!(Gene_net %in% test)]

table(lost %in% Gene_list_Pre)
table(lost %in% Gene_list_Unpre)
table(lost %in% Gene_grey)

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
  geom_violin(data = rbind_final, aes(x = Cate,y = Prop,fill=Region))+
  #geom_violin(alpha =.8,width = .8) +
  #geom_boxplot(outlier.alpha = 0.2,outlier.size = 0.3) +
  xlab("Gene Preservation Status") + ylab('') +
  theme(#legend.position="none",
        axis.title.x = element_text(size=12, face="bold"),
        axis.text.x = element_text(size=9,color ='black',face="bold"),
        axis.text.y = element_text(size=9,color ='black',face="bold"),
        strip.text = element_text(size=12,color ='black',face="bold"))+
  facet_wrap(~Region,labeller = variable_labeller)
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
library(tidyverse)
# tiff("Fig4-Gene-TF-Methy-by-Cate.tiff",
#      width = 14, height = 8, units = 'in', res = 300)
# ggarrange(P1,P2,labels = c("A", "B"),ncol=1,nrow =2)
# #grid.arrange(P1,P2,ncol=1,nrow =2)
# # plot_grid(P1,P2,align = c("v"),
# #            labels = c("A","B"), label_size= 12,label_colour = "black")
# dev.off()

test = Gene_Meth_Viol %>% 
  dplyr::filter(TFinfo =='TF') %>% 
  group_by(Gene) %>% sample_n(1)
length(unique(test$Gene))
table(test$Cate)

Gene_Meth_Viol_tmp_Gene = Gene_Meth_Viol %>% 
  dplyr::select(Gene,Prop_Prpt,Prop_All,Cate) %>% 
  mutate(Cate=recode(Cate, `Pre`="Preserved(5589)",`Unp`="Unpreserved(1146)"))

Gene_Meth_Viol_tmp_Gene_dup = Gene_Meth_Viol_tmp_Gene %>% 
  mutate(Cate=recode(Cate, `Preserved(5589)`="All(6735)",`Unpreserved(1146)`="All(6735)")) %>% 
  rbind(Gene_Meth_Viol_tmp_Gene) %>% 
  pivot_longer(cols = c(Prop_All,Prop_Prpt),
               names_to = "Cat", values_to = "Prop") %>% 
  mutate(Source = 'Genes')

Gene_Meth_Viol_tmp_TF = Gene_Meth_Viol %>% 
  dplyr::filter(TFinfo != 'NOT') %>%
  dplyr::select(Gene,Prop_Prpt,Prop_All,Cate) %>% 
  mutate(Cate=recode(Cate, `Pre`="Preserved(253)",`Unp`="Unpreserved(39)"))

Gene_Meth_Viol_tmp_TF_dup = Gene_Meth_Viol_tmp_TF %>% 
  mutate(Cate=recode(Cate, `Preserved(253)`="All(292)",`Unpreserved(39)`="All(292)")) %>% 
  rbind(Gene_Meth_Viol_tmp_TF) %>% 
  pivot_longer(cols = c(Prop_All,Prop_Prpt),
               names_to = "Cat", values_to = "Prop") %>% 
  mutate(Source = 'Transcription Factors')


Meth_Viol_plot = rbind(Gene_Meth_Viol_tmp_Gene_dup,
                       Gene_Meth_Viol_tmp_TF_dup)

##### fig 4 new ####
Gene_Meth_Viol_tmp_TF_new = Gene_Meth_Viol %>% 
  dplyr::filter(TFinfo != 'NOT') %>%
  dplyr::select(Gene,Prop_Prpt,Prop_Body,Cate) %>% 
  pivot_longer(cols = c(Prop_Prpt,Prop_Body),
               names_to = "Region", values_to = "Prop") %>% 
  mutate(Source = 'Transcription Factor')

library(gplots)
P4_Meth_bycate_new = Gene_Meth_Viol %>% 
  dplyr::select(-c(TFinfo,Module,Count_Prpt.y,Count_Body.y,Count_All.y,Prop_All)) %>% 
  pivot_longer(cols = c(Prop_Prpt,Prop_Body),
               names_to = "Region", values_to = "Prop") %>% 
  mutate(Source = 'Gene') %>% 
  rbind(Gene_Meth_Viol_tmp_TF_new) %>% 
  dplyr::filter(Prop<= .2) %>% 
  mutate(fill_plot = as.character(recode(Cate,Pre = col2hex('red'),Unp = col2hex('blue'))))
colorindex_viol = as.character(P4_Meth_bycate_new$fill_plot)


variable_names <- list(
  'Prop_Body' = 'Gene Body',
  'Prop_Prpt' = 'Regulatory Region',
  'Gene' = 'Genes',
  'Transcription Factor' = 'Transcription Factors')

variable_labeller <- function(variable,value){ return(variable_names[value]) } 


col2hex(c('red','blue'))

P4_Meth_bycate = 
ggplot() +
  geom_violin(data = P4_Meth_bycate_new,
              aes(x = Region,y = Prop,fill = Cate))+
    xlab('') + ylab('Methylation Level') + 
    facet_grid(Source ~ Region,scales = "free",labeller=variable_labeller)+
  theme(legend.position='bottom',
        # Change legend key size and key width
        legend.key.size = unit(0.6, "cm"),
        legend.key.width = unit(0.6,"cm"),
    #legend.position="none",
    #legend.position = c(0.7, 0.2)
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=6,color ='black'),
    legend.title = element_text(colour="black", size=10,face="bold"),
    legend.text = element_text(colour="black", size=14,face="bold"),
    legend.text.align = 0,
    strip.text = element_text(size=12,color ='black',face="bold"),
    strip.background = element_rect(colour="black", size=1)) + 
  scale_fill_manual(values = c("red","blue"),labels = c('Preserved Modules','Unpreserved Modules'),
                    name= " ", guide = guide_legend(reverse = F))
P4_Meth_bycate

#scale_fill_discrete(name = "", labels = c('Preserved','Unpreserved'))
  #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

tiff("Fig4-Gene-TF-Methy-by-Cate.tiff",
     width = 14, height = 8, units = 'in', res = 500)
print(P4_Meth_bycate)
dev.off()

################  Figure 0 - background plots ####################
#mat <- matrix(c(1,1,1,2,2,2,rep(3,6)), nrow = 2, byrow = TRUE)
tiff("Fig0-Net-Constr-Background.tiff",width = 14, height = 8, units = 'in', res = 500)
mat <- matrix(c(1,1,1,2,2,2,rep(3,6)), nrow = 2, byrow = TRUE);layout(mat);plot(sft_b_cl$fitIndices[,1], -sign(sft_b_cl$fitIndices[,3])*sft_b_cl$fitIndices[,2],
     cex.main = 1.5,cex.lab=1.2,font.lab=2,
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit",type="n",
     main = paste("Scale Independence"));text(sft_b_cl$fitIndices[,1], -sign(sft_b_cl$fitIndices[,3])*sft_b_cl$fitIndices[,2],
                                                      labels=powers,cex=0.9,col="red");abline(h=0.80,col="red");mtext(side=3, line=-2.4, at=0.05,text="A.", cex = 1,adj=0, outer=T);plot(sft_b_cl$fitIndices[,1], sft_b_cl$fitIndices[,5],
     cex.main = 1.5,cex.lab=1.2,font.lab=2,
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",type="n",
     main = paste("Mean connectivity"));text(sft_b_cl$fitIndices[,1], sft_b_cl$fitIndices[,5], labels=powers, cex=.9,col="red");abline(h=MeanK_b,col="red");par(mar = c(1,4.1,4.1,2));mtext(side=3,line=-2.4,at=0.5,cex = 1,text="B.", adj=0, outer=T);plot(METree_control, 
     cex.main = 1.5,cex.lab=1.2,font.lab=2,
     main = "Clustering of Initial MEs (Control Diet)",
     xlab = '', sub = "", ylab = "Height",
     hang = -1,cex = 0.6);MEDissThres = 0.2;abline(h = MEDissThres, col = "red");mtext(side=3, line=-32, at=0.05,cex = 1,text="C.", adj=0, outer=T);
#axis(side = 2, at = seq(0,1,.2), col = "#F38630",labels = FALSE, lwd = 2)
dev.off()
# 
# plot(rnorm(10))
# mtext(side=3, line=-3,at=0.05,
#       text="Here again?", adj=0, outer=T)
getwd()

################  Files -  ####################
################  1. gene measure ments
################  2. associations related

# 1. meth prop vs DCG and inCon
dcg_raw = DCGL::WGCNA((t(datExpr_control)),(t(datExpr_treatment)),power = 24, variant = "DCp")

Gene_deg_ME_count_out = Diff_Coexp[,c(1,17:19,2:16,20:25)]
head(Gene_deg_ME_count_out)

#2 module preservation 
head(MP_Stats_final)

# require(openxlsx)
# list_of_datasets <-
#   list("Measurements by Genes" = Diff_Coexp,
#        "Module Preservation Statistics" = MP_Stats_final)
# write.xlsx(list_of_datasets, file = "Gene_Measure_and_MP_Stats.xlsx")

#3 loci matching

# done!
# output module preservation statistics (excel)

################  Figure 5 - linear relationship seperate ####################
Diff_Coexp_plot1 = Diff_Coexp %>%
  dplyr::select(-Module,-Prop_Prpt) %>% 
  pivot_longer(cols = c(DCG,kWithin,kTotal,`Module Membership`), #,Prop_Body
               names_to = "Kind", values_to = "Gene Measurement") %>% 
  mutate(Region = 'All') %>% 
  rename(Proportion = Prop_All)

Diff_Coexp_plot2 = Diff_Coexp %>%
  dplyr::select(-Module,-Prop_All) %>% 
  pivot_longer(cols = c(DCG,kWithin,kTotal,`Module Membership`), #,Prop_Body
               names_to = "Kind", values_to = "Gene Measurement") %>% 
  mutate(Region = 'Promoter')%>% 
  rename(Proportion = Prop_Prpt)

Diff_Coexp_plot = rbind(Diff_Coexp_plot1,Diff_Coexp_plot2) %>% 
  rename(`Methylation Proportion` = Proportion) %>% 
  mutate(alpha = ifelse(Cate == 'Pre','1','.9'))

#

Diff_Coexp_plot_new = Diff_Coexp_plot %>% 
  dplyr::filter(Kind == 'kWithin')

library(ggpmisc)
library(ggpubr)
library(broom)
#my.format ="b[0]~`=`~%.3g*\", \"*b[1]~`=`~%.3g*\""

my.format <-
  "b[0]~`=`~%.3g*\",\"*b[1]~`=`~%.3g"

variable_names_lm <- list(
  'All' = 'Gene Body',
  'Promoter' = 'Regulatory Region')

variable_labeller_lm <- function(variable,value){ return(variable_names_lm[value]) } 

scale_fill_discrete2 <- function(...) {
  scale_colour_manual(..., values = col2hex(c('red','blue')))
}


head(Diff_Coexp_plot_new)

genemeasure_vs_prop =
  ggplot(Diff_Coexp_plot_new,aes(x = `Methylation Proportion`,
                             y = `Gene Measurement`,color = Cate),) + 
  geom_point(aes(alpha = alpha),size = 1) +
  facet_wrap(~Region,ncol = 2,scales = 'free',labeller=variable_labeller_lm)+
  coord_flip() +
  geom_smooth(method = "lm",linetype="dashed",
              aes(group = factor(Cate),colour = factor(Cate)),
              formula =y ~ x,se = F)+
  ylab('Intramodular Connectivity')+
  xlab('Methylation Level')+
  # stat_fit_tidy(method = "lm",
  #               label.x = .8,
  #               vstep = .4,
  #               hstep = .05,
  #               method.args = list(formula = y ~ x),
  #               mapping = aes(size = 3,
  #                             label = sprintf("slope = %.1e\np-val = %.1e",
  #                                             stat(x_estimate),
  #                                             stat(x_p.value))))+
  theme(legend.position='bottom',
        # Change legend key size and key width
        legend.key.size = unit(.6, "cm"),
        legend.key.width = unit(.6,"cm"),
        legend.title = element_text(colour="black", size=10,face="bold"),
        legend.text = element_text(colour="black", size=14,face="bold"),
        legend.text.align = 0,
        axis.title.x = element_text(size=17,face="bold"),
        axis.title.y = element_text(size=17,face="bold"),
        axis.text.x = element_text(size=9,color ='black',face="bold"),
        axis.text.y = element_text(size=9,color ='black',face="bold"),
        strip.text = element_text(size=12,color ='black',face="bold")) +
  scale_alpha_manual(values=c(1,0.3),guide=F)+
  #scale_color_manual(values=col2hex(c('red','blue')))+
  scale_fill_discrete2(name = "",labels = c('Preserved Modules','Unpreserved Modules'))

genemeasure_vs_prop


col2hex(c('red','blue'))
cols2 <- c("Pre" = "red", "Unp" = "blue")
genemeasure_vs_prop

#
tiff("Fig5-Gene-Measure-vs-MethyProp.tiff",
     width = 16, height = 10, units = 'in', res = 500)
print(genemeasure_vs_prop)
dev.off()

##
DCG_Pre = Diff_Coexp %>% dplyr::filter(Cate == 'Pre') %>% dplyr::select(DCG) %>% unlist(use.names = F)
kTotal_Pre = Diff_Coexp %>% dplyr::filter(Cate == 'Pre') %>% dplyr::select(kTotal) %>% unlist(use.names = F)
kWithin_Pre = Diff_Coexp %>% dplyr::filter(Cate == 'Pre') %>% dplyr::select(kWithin) %>% unlist(use.names = F)
MM_Pre = Diff_Coexp %>% dplyr::filter(Cate == 'Pre') %>% dplyr::select(`Module Membership`) %>% unlist(use.names = F)
#
PropPrmp_Pre = Diff_Coexp %>% dplyr::filter(Cate == 'Pre') %>% dplyr::select(Prop_Prpt) %>% unlist(use.names = F)
PropAll_Pre = Diff_Coexp %>% dplyr::filter(Cate == 'Pre') %>% dplyr::select(Prop_All) %>% unlist(use.names = F)

##
DCG_Unp = Diff_Coexp %>% dplyr::filter(Cate == 'Unp') %>% dplyr::select(DCG) %>% unlist(use.names = F)
kTotal_Un = Diff_Coexp %>% dplyr::filter(Cate == 'Unp') %>% dplyr::select(kTotal) %>% unlist(use.names = F)
kWithin_Unp = Diff_Coexp %>% dplyr::filter(Cate == 'Unp') %>% dplyr::select(kWithin) %>% unlist(use.names = F)
MM_Unp = Diff_Coexp %>% dplyr::filter(Cate == 'Unp') %>% dplyr::select(`Module Membership`) %>% unlist(use.names = F)
#
PropPrmp_Unp = Diff_Coexp %>% dplyr::filter(Cate == 'Unp') %>% dplyr::select(Prop_Prpt) %>% unlist(use.names = F)
PropAll_Unp = Diff_Coexp %>% dplyr::filter(Cate == 'Unp') %>% dplyr::select(Prop_All) %>% unlist(use.names = F)

library(broom)
v1 = get('kWithin_Pre')
v2 = get('kWithin_Unp')
v3 = get('PropPrmp_Pre')
v4 = get('PropAll_Pre')
v5 = get('PropPrmp_Unp')
v6 = get('PropAll_Unp')

df = data.frame(v1 = get('kWithin_Unp'),
                v2 = get('PropAll_Unp'))
lmod   = lm(v1~v2,df)
summary(lmod)




################  Figure 6 - select TF by controling percentage ####################
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

table(Bta_TF_OverlapMatch$Pres)

Bta_TF_Mstatus_final = Bta_TF_Mstatus_raw %>%  
  dplyr::select(-c('Count_Body.x',#'Count_All.x
                   'Count_Prpt.y','Count_Body.y','Count_All.y',)) %>% 
  group_by(Suggested.Symbol.x) %>% 
  mutate(TcoF_meth = mean(Prop_All.y)) %>% 
  dplyr::select(-Prop_Body.y,-Prop_Prpt.y,-Prop_All.y,-SymbExist,-TcoF_ID,-Suggested.Symbol.y) %>% 
  sample_n(1) %>% 
  mutate(Prop_All_wTcoF = ifelse(is.na(Prop_All.x),TcoF_meth,TcoF_meth + Prop_All.x))


Bta_TF_OverlapMatch_plot = Bta_TF_OverlapMatch %>% 
  left_join(Bta_TF_Mstatus_final,by = c('TF_Name'='Suggested.Symbol.x')) %>% 
  replace_na(list(Prop_All_wTcoF = 0)) %>% 
  pivot_longer(cols = c(Prop_Prpt.x,Prop_All.x),
               names_to = "Kind", values_to = "Prop") %>% 
  mutate(Pres=recode(Pres, `Pre`="Preserved(113)",`Unpre`="Unpreserved(199)"))


kstest  = Bta_TF_OverlapMatch %>% 
  left_join(Bta_TF_Mstatus_final,by = c('TF_Name'='Suggested.Symbol.x')) %>% 
  replace_na(list(Prop_All_wTcoF = 0)) %>% 
  mutate(Pres=recode(Pres, `Pre`="Preserved(113)",`Unpre`="Unpreserved(199)"))


#
BT_Unp_all = kstest %>% 
  ungroup() %>% 
  dplyr::filter(Pres =='Unpreserved(199)') %>% 
  dplyr::select(Prop_All_wTcoF) %>% 
  unlist(use.names = F)


BT_Pre_all = kstest %>% 
  ungroup() %>% 
  dplyr::filter(Pres =='Preserved(113)') %>% 
  dplyr::select(Prop_All_wTcoF) %>% 
  unlist(use.names = F)


BT_Unp_prmp = kstest %>% 
  ungroup() %>% 
  dplyr::filter(Pres =='Unpreserved(199)') %>% 
  dplyr::select(Prop_Prpt.x) %>% 
  unlist(use.names = F)


BT_Pre_prmp = kstest %>% 
  ungroup() %>% 
  dplyr::filter(Pres =='Preserved(113)') %>% 
  dplyr::select(Prop_Prpt.x) %>% 
  unlist(use.names = F)

aa = ks.test(BT_Unp_all,BT_Pre_all)
bb = ks.test(BT_Unp_prmp,BT_Unp_prmp)
aa$p.value
bb$p.value

#####
anno3 <- data.frame(xstar = c(1.5), 
                    ystar = c(.2),
                    lab = c("KS-test p-val: 0.775"),
                    Source = c("Genes", "Transcription Factors"))

anno4 <- data.frame(xstar = c(1.5), 
                    ystar = c(.15),
                    lab = c("KS-test p-val: 1"),
                    Source = c("Genes", "Transcription Factors"))

Force_assign = 
ggplot() + 
  geom_violin(data = Bta_TF_OverlapMatch_plot,
              aes(x = Pres,y =Prop,fill = Kind))+
  #geom_violin(alpha =.8,width = .8) +
  #geom_boxplot(outlier.alpha = 0.2,outlier.size = 0.3) +
  xlab('Preservation Category') + ylab('Proportion of DMCs') +
  theme(legend.position='bottom',
        # Change legend key size and key width
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"),
        #legend.position="none",
        #legend.position = c(0.7, 0.2)
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(size=12,color ='black'),
        axis.text.y = element_text(size=9,color ='black'),
        strip.text = element_text(size=12,color ='black',face="bold"),
        legend.title = element_text(colour="black", size=10,face="bold"),
        legend.text = element_text(colour="black", size=8,face="bold"),
        legend.text.align = 0) +
  labs(fill = "Region") +
  scale_fill_discrete(name = "Region", labels = c(expression('Methprop'['GENE']),expression('Methprop'['REG'])))+
  geom_text(data = anno3,
            aes(x = xstar,  y = ystar, label = lab),
            size=5,
            color = '#F8766D', fontface = "bold")+
  geom_text(data = anno4,
            aes(x = xstar,  y = ystar, label = lab),
            size=5,
            color = '#00BFC4', fontface = "bold")
Force_assign


tiff("Fig6-forcely-assgin-tf.tiff",
     width = 14, height = 8, units = 'in', res = 500)
print(Force_assign)
dev.off()

#### 
ks.test(male,female)
wilcox.test(male,female) 



Gene_Unp_All = Gene_Meth_Viol %>% 
  dplyr::filter(Cate == 'Unp',TFinfo == 'NOT') %>% 
  dplyr::select(Prop_Body) %>% 
  unlist(use.names = F)

Gene_Unp_Prmp = Gene_Meth_Viol %>% 
  dplyr::filter(Cate == 'Unp',TFinfo == 'NOT') %>% 
  dplyr::select(Prop_Prpt) %>% 
  unlist(use.names = F)

Gene_Pre_All = Gene_Meth_Viol %>% 
  dplyr::filter(Cate == 'Pre',TFinfo == 'NOT') %>% 
  dplyr::select(Prop_Body) %>% 
  unlist(use.names = F)

Gene_Pre_Prmp = Gene_Meth_Viol %>% 
  dplyr::filter(Cate == 'Pre',TFinfo == 'NOT') %>% 
  dplyr::select(Prop_Prpt) %>% 
  unlist(use.names = F)

TF_Unp_All  = Gene_Meth_Viol %>% 
  dplyr::filter(Cate == 'Unp',TFinfo != 'NOT') %>% 
  dplyr::select(Prop_Body) %>% 
  unlist(use.names = F)

TF_Unp_Prmp = Gene_Meth_Viol %>% 
  dplyr::filter(Cate == 'Unp',TFinfo != 'NOT') %>% 
  dplyr::select(Prop_Prpt) %>% 
  unlist(use.names = F)

TF_Pre_All = Gene_Meth_Viol %>% 
  dplyr::filter(Cate == 'Pre',TFinfo != 'NOT') %>% 
  dplyr::select(Prop_Body) %>% 
  unlist(use.names = F)

TF_Pre_Prmp = Gene_Meth_Viol %>% 
  dplyr::filter(Cate == 'Pre',TFinfo != 'NOT') %>% 
  dplyr::select(Prop_Prpt) %>% 
  unlist(use.names = F)

# KS test
a = ks.test(Gene_Pre_All,Gene_Unp_All)
b = ks.test(Gene_Pre_Prmp,Gene_Unp_Prmp)
a$p.value
b$p.value

c = ks.test(TF_Pre_All,TF_Unp_All)
d = ks.test(TF_Pre_Prmp,TF_Unp_Prmp)
c$p.value
d$p.value


library(grid)
# Create a text
grob <- grobTree(textGrob("Scatter plot", x=0.1,  y=0.95, hjust=0,
                          gp=gpar(col="red", fontsize=13, fontface="italic")))
# Plot
# anno1 <- data.frame(x1 = c(1.75, 0.75), 
#                    x2 = c(2.25, 1.25), 
#                    y1 = c(1, 1), 
#                    y2 = c(1, 1), 
#                    xstar = c(2.5,2.58), 
#                    ystar = c(.7, .7),
#                    lab = c("KS-test p-val: 1.1e-16", "KS-test p-val: p-val: 1.5e-3"),
#                    Source = c("Genes", "Transcription Factors"))
# 
# anno2 <- data.frame(x1 = c(1.75, 0.75), 
#                    x2 = c(2.25, 1.25), 
#                    y1 = c(1, 1), 
#                    y2 = c(1, 1), 
#                    xstar = c(2.5,2.55), 
#                    ystar = c(.4, .4),
#                    lab = c("KS-test p-val: 0.735", "KS-test p-val: p-val: 1"),
#                    Source = c("Genes", "Transcription Factors"))





