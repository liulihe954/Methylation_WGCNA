#install.packages('customLayout')
library(ggplot2);library(customLayout)

# 创建拼图画布
lay1 <- lay_new( 
mat = matrix(1:3, ncol = 1), # 矩阵分布，mat表示指定排版的数字矩阵 
widths = c(4,2),# 设定宽度比例
heights = c(5,1))             
# 设置高度比例)
lay_show(lay1) # 显示拼图画布

# 创建第2个拼图画布，与第1个结构一样，只是比例不一样
lay2 <- lay_new(matrix(1:4, nc = 2), widths = c(3, 5),heights = c(2, 4))
lay_show(lay2)
rm(lay1)
lay1 <- lay_new( 
  mat = matrix(1:3, ncol = 1), # 矩阵分布，mat表示指定排版的数字矩阵 
  widths = c(4,2),# 设定宽度比例
  heights = c(5,1))

par(mfrow=c(3,1))

par(mar = c(2,4.1,4.1,2))
plot(METree_control, 
     main = "Clustering of Initial MEs (Control Diet)",
     xlab = '', sub = "", ylab = "Height",
     hang = -1,cex = 0.35);MEDissThres = 0.2;abline(h = MEDissThres, col = "red")
dev.off()

par(mar = c(2,4.4,4.1,2))
#par(mar = c(1,1,1,1))
#par(mfrow=c(4,2))
plotDendroAndColors(geneTree_control,
                    cbind(dynamicColors_control,moduleColors_control), 
                    c("Dynamic Tree Cut", "Merged Modules"), 
                    dendroLabels=F, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, setLayout = F,
                    main="Gene Cluster Dendrogram (Control Diet)")
dev.off()

plotDendroAndColors(geneTree_treatment,moduleColors_control,
                    "Modules", 
                    dendroLabels=F,hang=0.03, addGuide=TRUE,
                    guideHang=0.05,setLayout = F,
                    main="Gene Cluster Dendrogram (Methionine Diet)")
