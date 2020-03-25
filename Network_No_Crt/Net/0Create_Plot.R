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

plotDendroAndColors(geneTree_control,
                    cbind(dynamicColors_control,moduleColors_control), 
                    c("Dynamic Tree Cut", "Merged Modules"), 
                    dendroLabels=F, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, setLayout = F,
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
                                     guideHang=0.05, 
                                     main="Gene Cluster Dendrogram (Control Diet)"))
par(mar = c(2,4.4,4.1,2))
DengroTreatment = as.ggplot(~plotDendroAndColors(geneTree_treatment,moduleColors_control,
                                                  "Modules", 
                                                  dendroLabels=F,hang=0.03, addGuide=TRUE,
                                                  guideHang=0.05,
                                                  main="Gene Cluster Dendrogram (Methionine Diet)"))


tiff("Fig1-Net-Dengdro.tiff", width = 20, height = 16, units = 'in', res = 300)
DengroControl / DengroTreatment
# plot_grid(DengroControl,DengroTreatment,
#           align = c("v"),
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


tiff("Fig2-Net-Presv-stats.tiff", width = 10, height = 8, units = 'in', res = 300)
PresZsum / PresMedianR
# plot_grid(DengroControl,DengroTreatment,
#           align = c("v"),
#           labels = c("A","B"), label_size= 20,label_colour = "black")
dev.off()


####################  Fig3 Functional Enrich ###################



