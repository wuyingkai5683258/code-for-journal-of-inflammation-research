#gsva
library(ggplot2)
library(limma)
library(ComplexHeatmap) 
library(clusterProfiler) 
library(GSVA) 
library(GSEABase)
library(dplyr) 
library(data.table) 
library(tidyverse)
COAD <- read.csv("expmerge.csv", check.names = FALSE, header = TRUE, row.names = 1)
#COAD <- log2(COAD+1)
gene_set <- getGmt("c5.go.v2023.2.Hs.symbols.gmt")
group_list <-  read.csv("groupmerge.csv", check.names = FALSE, header = TRUE, row.names = 1)
table(group_list)
annotation <- data.frame(group_list)
rownames(annotation) <- colnames(COAD)
head(annotation)
gsva_result<- gsva(as.matrix(COAD), gene_set, method = "gsva",min.sz=1,
                   max.sz=Inf,kcdf="Gaussian")
library(dendextend)
library(circlize)
library(RColorBrewer)
colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
values <- seq(-0.8, 0.8, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)
top_annotation<-HeatmapAnnotation(df=annotation,col=list(group=c("NC"="blue","LDH"="red")))
Heatmap(gsva_result, name = "GSVA", col = col_fun,cluster_rows = T,cluster_columns = F,show_row_names = T,
        show_column_names = F,column_split = annotation$group,)
annotation <- annotation[order(annotation$group == "LDH", decreasing = TRUE), , drop = FALSE]
exp="LDH"
ctr="NC"
design <- model.matrix(~0+factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_result)
contrast.matrix <- makeContrasts(contrasts=paste0(exp,'-',ctr),  #"exp/ctrl"
                                 levels = design)

fit1 <- lmFit(gsva_result,design)                 
fit2 <- contrasts.fit(fit1, contrast.matrix)
efit <- eBayes(fit2)                         

summary(decideTests(efit,lfc=1, p.value=0.05)) 
tempOutput <- topTable(efit, coef=paste0(exp,'-',ctr), n=Inf)
degs <- na.omit(tempOutput) 
keep <- rownames(degs[degs$adj.P.Val<1& degs$P.Value<1, ])
length(keep)
dat <- gsva_result[keep,] 
pheatmap(dat,cluster_rows = T,cluster_cols = F,
         color=colorRampPalette(c("navy","white","firebrick3"))(100),
         show_colnames = F,border_color = NA,scale = "row",show_rownames =T,
         annotation_col = group_list,fontsize=5.0,name = "GSVA")
dat=t(dat)
write.csv(dat,"gsva.csv")

################
library(ggstatsplot)
library(ggplot2)
dt<-read.csv("boxplot of gsva.csv",sep = ",",check.names = T,row.names = 1)
p1 <- ggbetweenstats(
  data  = dt,
  x= group,
  y= Obesity,
  plot.type = "box",#图标类型，不加参数则默认为boxviolin，其它可选为：box、violin
  title = ""
)
p2 <- ggbetweenstats(
  data  = dt,
  x     = group,
  y     = Obesity,
  plot.type = "boxviolin",
  title = "GSVA score of Obesity",
  results.subtitle = TRUE,#决定是否将统计检验的结果显示为副标题（默认TRUE）;如果设置为FALSE,则仅返回绘图
  subtitle = NULL,#副标题,默认显示统计检验结果,自定义则results.subtitle=FALSE
  outlier.tagging = TRUE,#是否标记离群异常值，默认FALSE
  outlier.shape = 19,#异常值形状,可设置为NA将其隐藏（不是删除，因此不会影响统计检验结果）
  outlier.color = "pink",#异常值颜色
  outlier.label.args = list(size = 4),#异常值标签大小
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6),
                    alpha = 0.5, size = 3, stroke = 0),#传递给geom_point的参数设置
  violin.args = list(width = 0.4, alpha = 0.2),#传递给geom_violin的参数设置
  ggtheme = theme_bw(),#主题修改，可直接调用ggplot2的主题，默认主题为ggstatsplot::theme_ggstatsplot()
  package = "ggsci",#提取调色板所需的包
  palette = "uniform_startrek"#选择提取包中的调色板
)
p2





library(WGCNA)
#fpkm = read.csv("exp70362.csv",header=T,sep = ",",check.names=F,row.names = 1)#########file name can be changed#####数据文件名，根据实际修改，如果工作路径不是实际数据路径，需要添加正确的数据路径
#fpkm=log2(fpkm+1)
#write.csv(fpkm,"fpkm.csv")
fpkm = read.csv("expmerge.csv",header=T,sep = ",",check.names=F,row.names = 1)
cibersort="gsva.csv"
meanFPKM=0.5
#
#rownames(fpkm)=fpkm[,1]
#dim(fpkm)
#names(fpkm)
datExpr0 =fpkm# as.data.frame(fpkm[,-1])


dataExpr=t(datExpr0[order(apply(datExpr0, 1, mad), decreasing = TRUE)[1:10000], ])
gsg <- goodSamplesGenes(dataExpr, verbose = 3)
# 如果存在异常样本或基因
if (!gsg$allOK) {
  # 异常的基因
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  # 异常的样本
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # 删除异常的样本和基因
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}



# 再次检查是否所有基因和样本均合格
gsg2 <- goodSamplesGenes(dataExpr, verbose = 3)
gsg$allOK 

datExpr0<-dataExpr
sampleTree = hclust(dist(datExpr0), method = "average")


pdf(file ="1_sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

### Plot a line to show the cut
##abline(h = 15, col = "red")##是否选择剪切
dev.off()


### Determine cluster under the line
##clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
##table(clust)


### clust 1 contains the samples we want to keep.
##keepSamples = (clust==1)
##datExpr0 = datExpr0[keepSamples, ]


####载入性状数据
#Loading clinical trait data
traitData=read.csv(cibersort,header = T,sep = ",",check.names = F,row.names = 1)



for (df in colnames(traitData)) {
  traitData[,df]=traitData[,df]/max(traitData[,df])
  print(sd(traitData[,df]))
}
max(traitData)
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
fpkmSamples = rownames(datExpr0)
traitSamples =rownames(allTraits)
traitRows = match(fpkmSamples, traitSamples)
datTraits = allTraits[traitRows,]
rownames(datTraits) 
collectGarbage()

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.

#sizeGrWindow(12,12)
pdf(file="2_Sample dendrogram and trait heatmap.pdf",width=12,height=11)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap", 
                    marAll = c(1, 12, 3, 1))
dev.off()







#############################network constr########################################
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Choose a set of soft-thresholding powers
powers = c(1:30)

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# Plot the results:
#sizeGrWindow(9, 5)
pdf(file="3_Scale independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


######chose the softPower
softPower =sft$powerEstimate
#显示软阈值
print(softPower)

adjacency = adjacency(datExpr0, power = softPower)

##### Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)

#sizeGrWindow(12,9)
pdf(file="4_Gene clustering on TOM-based dissimilarity.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
pdf(file="5_Dynamic Tree Cut.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
#sizeGrWindow(7, 6)
pdf(file="6_Clustering of module eigengenes.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25######模块剪切高度
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()


# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

#sizeGrWindow(12, 9)
pdf(file="7_merged dynamic.pdf", width = 9, height = 6.5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

# Save module colors and labels for use in subsequent parts
#save(MEs, TOM, dissTOM,  moduleLabels, moduleColors, geneTree, sft, file = "networkConstruction-stepByStep.RData")




##############################relate modules to external clinical triats######################################
# Define numbers of genes and samples
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#sizeGrWindow(10,6)
pdf(file="8_Module-trait relationships.pdf",width=4,height=9)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(12, 8, 3, 3))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


######## Define variable weight containing all column of datTraits

###MM and GS


# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#names of those trait
traitNames=names(datTraits)

geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")


####plot MM vs GS for each trait vs each module


##########example:royalblue and CK
#module="royalblue"
#column = match(module, modNames)
#moduleGenes = moduleColors==module

#trait="CK"
#traitColumn=match(trait,traitNames)

#sizeGrWindow(7, 7)

######

for (trait in traitNames){
  traitColumn=match(trait,traitNames)
  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      
      #sizeGrWindow(7, 7)
      pdf(file=paste("9_", trait, "_", module,"_Module membership vs gene significance.pdf",sep=""),width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}

#####
datExpr0<-as.data.frame(datExpr0)
names(datExpr0)
probes = names(datExpr0)


#################export GS and MM############### 

geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)

for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]

write.table(geneInfo, file = paste0("10_GS_and_MM.xls"),sep="\t",row.names=F)


trait="Obesity"
traitColumn=match(trait,traitNames)  
for (module in modNames){
  column = match(module, modNames)
  moduleGenes = moduleColors==module
  if (nrow(geneModuleMembership[moduleGenes,]) > 1){
    outPdf=paste(trait, "_", module,".pdf",sep="")
    pdf(file=outPdf,width=7,height=7)
    par(mfrow = c(1,1))
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, traitColumn]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for ",trait),
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module,
                       pch = 16,cex = 1.2)
    abline(v = 0.5, col = "red", lty = 2, lwd = 1.5) # MM=0.4
    abline(h = 0.5, col = "red", lty = 2, lwd = 1.5) # GS=0.1
    dev.off()
  }
}

probes = colnames(dataExpr)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.xls",sep="\t",row.names=F)

for (mod in 1:nrow(table(moduleColors))){  
  modules = names(table(moduleColors))[mod]
  probes = colnames(dataExpr)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("module_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

















###########
#############
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(limma)
expFile="expmerge.csv"          #?????ļ?
#??ȡ?????ļ????????????ļ?????
expr_data=read.csv(expFile,sep=",",header=T,row.names = 1,check.names=F)
group<-read.csv("groupmerge.csv",header = T,row.names = 1,sep = ",")
design <- model.matrix(~0+factor(group$group))
colnames(design) <- levels(factor(group$group))
rownames(design) <- colnames(expr_data)
contrast.matrix <- makeContrasts(LDH-NC,levels = design) 
fit <- lmFit(expr_data,design) 
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
DEG$regulate <- ifelse(DEG$adj.P.Val> 0.05, "unchanged",
                       ifelse(DEG$logFC > 1, "up-regulated",
                              ifelse(DEG$logFC < -1, "down-regulated", "unchanged")))
write.table(DEG,"差异基因.csv",row.names=T,col.names=T,sep=",")

DE_1_0.05 <- DEG[DEG$adj.P.Val<0.05&abs(DEG$logFC)>1,]
upGene_1_0.05 <- DE_1_0.05[DE_1_0.05$regulate == "up-regulated",]
downGene_1_0.05 <- DE_1_0.05[DE_1_0.05$regulate == "down-regulated",]
write.csv(DE_1_0.05,"DEGused.csv")
#write.csv(downGene_1_0.05,"downGene_1_005.csv")
############
library(dplyr)
library(ggfun)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(tidyverse)
library(ggrepel)

#读取文件
{
  deg_result <- read.csv(file = "差异基因.csv") %>%
    dplyr::mutate(rank = -1 * rank(logFC, ties.method = "max"))
}

## 手动选择的基因“CTSL”

# ## 上调 top2
top_2 <- deg_result %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::arrange(desc(logFC)) %>%
  dplyr::slice_head(n =5)
# 
# ## 下调2
tail_2 <- deg_result %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::arrange(desc(logFC)) %>%
  dplyr::slice_tail(n = 5)

## 绘制火山图
deg_result %>%
  ggplot() + 
  geom_point(aes(x = rank, y = logFC, color = P.Value, size = abs(logFC))) + 
  
  scale_color_gradient2(low = '#2C7BB6', high = '#D7191C', mid = '#FFFFBF', 
                        midpoint = 0.5, name = "Pvalue") +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +  
  
  ## 添加手动选择的基因“CTSL”标签
  
  geom_text_repel(data = top_2,
                  aes(x = rank + 50, y = logFC, label = SYMBOL),
                  box.padding = 0.5, 
                  nudge_x = 10, 
                  nudge_y = 0.2, 
                  segment.curvature = -0.1,
                  segment.ncp = 3, 
                  segment.angle = 20,
                  direction = "y", 
                  hjust = "left",
                  max.overlaps = 20) +  # 增大值（根据数据调整，如20、30）
  
  # 第二次：tail_2的标签
  geom_text_repel(data = tail_2,
                  aes(x = rank + 10, y = logFC, label = SYMBOL),
                  box.padding = 0.5, 
                  nudge_x = 10, 
                  nudge_y = -0.2,
                  segment.curvature = -0.1, 
                  segment.ncp = 3,
                  segment.angle = 20, 
                  direction = "y", 
                  hjust = "right",
                  max.overlaps = 20) 






############vein
#devtools::install_github("jolars/eulerr")
library(ggvenn)
library(eulerr)
library(scales)
library(dplyr)  

# 读取数据
data <- read.csv("vein input.csv")

# 提取两列基因作为集合
comgene <- list(
  'DEGs' = unique(data$DEGs),
  'WGCNA' = unique(data$WGCNA)
)


# 计算交集基因
com_genes <- intersect(comgene$DEGs, comgene$WGCNA)
#保存为 CSV
{
  write.csv(data.frame(Intersection_Genes = com_genes),
            "交集基因.csv", row.names = FALSE)
}

#普通韦恩图
p1 <- ggvenn(
  comgene,
  show_percentage = F,
  show_elements = FALSE,
  label_sep = ",",
  digits = 2,
  stroke_color = NA,  # 去除边框
  fill_color = c("#C2E0F7", "#A4CBA8"), 
  set_name_color = c("#C2E0F7", "#A4CBA8"),
  text_color = "black",  # 调整文字颜色
  text_size = 5          # 增大文字
)
print(p1)

# 绘制比例图
p2 <- euler(comgene)
plot(
  p2,
  labels = list(col = "black", font = 1, cex = 1.5),  # 调整字体大小和颜色
  edges = NULL,  # 去除边框
  fills = c("#C2E0F7", "#A4CBA8"),  
  alpha = 0.9,
  quantities = list(cex = 1.5, col = 'black')  # 设置数量文字更突出
)








#############

VEIN=read.csv("交集基因.csv",sep=",",header=T) 
colnames(VEIN)
rownames(VEIN)
exp=read.csv("expmerge.csv",sep=",",check.names=F,header=T,row.names = 1)
DEG_gene_exprvein1<-exp[VEIN$Intersection_Genes,]
DEG_gene_exprvein1<-na.omit(DEG_gene_exprvein1)
DEG_gene_exprvein1<-t(DEG_gene_exprvein1)
write.table(DEG_gene_exprvein1,file="machinelearning.csv",sep=",",quote=T,col.names=T)

library(tidyverse)

library(glmnet)

source('msvmRFE.R') 
library(sigFeature)

library(e1071)

library(caret)

library(randomForest) 
library(plotmo)
###lasso
library(tidyverse)    
library(broom)
library(glmnet)



# 读取表达矩阵和分类文件，指定第一列作为行名
train<-read.csv("machinelearning.csv",row.names = 1)

set.seed(111)

x <- as.matrix(train[,-1])

(y <- ifelse(train$group == "NC", 0,1)) #把分组信息换成01

fit = glmnet(x, y, family = "binomial", alpha = 1, nlambda=100)
#10 validation
cvfit = cv.glmnet(x, y,
                  
                  nfold=10, #10-fold cross-validation
                  
                  family = "binomial", alpha=1)

pdf("2cvfit.pdf")

plot(cvfit)

dev.off()

#best lambda

cvfit$lambda.min

# 

myCoefs <- coef(cvfit, s="lambda.min")

lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]

(lasso_fea <- lasso_fea[-1])
write.csv(lasso_fea,"lasso.csv")
# 

pdf("1A_lasso.pdf")

plot(fit, xvar = "dev", label = TRUE)

dev.off()

library(broom)

# 提取数据，就是
# 
library(ggpubr)
par(mfrow=c(1,2))
plot(fit)
plot(cvfit)
###svm
input <- read.csv("machinelearning.csv",row.names = 1)
input$group<-ifelse(input$group == "NC", 0,1)
set.seed(123)
##下面开始正式分析
###采用五折交叉验证
svmRFE(input, k = 5, halve.above = 100)

nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100)

# #查看主要变量
top.features = WriteFeatures(results, input, save=F)
head(top.features)

#把SVM-REF找到的特征保存到文件
write.csv(top.features,"feature_svm.csv")

##变量越多，实际越长
##只用20变量个简单示范一下，大概运行2分钟

featsweep = lapply(1:38, FeatSweep.wrap, results, input) 

save(featsweep,file = "result.RData")
#load("result.RData")
#画图
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
par(mfrow=c(1,2))
#绘制错误率图形
PlotErrors(errors, no.info=no.info) 

#绘制准确率图形  
Plotaccuracy(1-errors,no.info=no.info) 



# 图中红色圆圈所在的位置，即错误率最低点
which.min(errors)

##保存核心基因集
top<-top.features[1:which.min(errors), "FeatureName"]

write.csv(top,"svm筛选的核心基因1.csv")
top.features=top.features[c(1:23),]
p1 = ggplot(top.features, aes(x = reorder(FeatureName, -AvgRank), y = AvgRank, fill = AvgRank)) +  
  geom_bar(stat = "identity", width = 0.7) +  
  coord_flip() +  
  scale_fill_gradient(low = ggsci::pal_npg()(2)[1], high = ggsci::pal_npg()(2)[2]) +  
  labs(x = "Gene", y = "E", title = "Genes listed by AvgRank") +  
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, hjust = 1), 
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  
p1
p2 = ggdotchart(top.features, x = "FeatureName", y = "AvgRank",
                color = "AvgRank", 
                sorting = "ascending",                        # 这里改为从小到大排序
                add = "segments",                             
                add.params = list(color = "lightgray", size = 2), 
                dot.size = 6,                        
                font.label = list(color = "white", size = 9,
                                  vjust = 0.5),               
                ggtheme = theme_bw()         ,               
                rotate = TRUE                                       )
p3 = p2 + geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  gradient_color(palette =c(ggsci::pal_npg()(2)[2],ggsci::pal_npg()(2)[1])      ) 
p3
####RANDOMFOREST
data<- read.csv("machinelearning.csv", header = T, sep=",",row.names =1)
data$group<- ifelse(data$group == "NC", 0,1)
set.seed(123)

y<-factor(data$group)
x<-data[,-39]

rf<-randomForest(x, y,data=data,importance=TRUE)
plot(margin(rf, data$group), main = '观测值被判断正确的概率图')
importance_otu <- data.frame(importance(rf), check.names = FALSE)
head(importance_otu)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
#根据表格输出前30变量
rf
varImpPlot(rf, n.var = min(30, nrow(rf$importance)), main = 'Top 30 - variable importance')
write.table(importance_otu, 'importance_otu.txt', sep = '\t', col.names = NA, quote = FALSE)
otu.cv <- replicate(5, rfcv(data[-ncol(data)], data$group, cv.fold = 10,step = 1.5), simplify = FALSE)
otu.cv
otu.cv <- data.frame(sapply(otu.cv, '[[', 'error.cv'))
otu.cv$otus <- rownames(otu.cv)
otu.cv <- reshape2::melt(otu.cv, id = 'otus')
otu.cv$otus <- as.numeric(as.character(otu.cv$otus))
otu.cv.mean <- aggregate(otu.cv$value, by = list(otu.cv$otus), FUN = mean)
head(otu.cv.mean, 14)
with(otu.cv.mean, plot(Group.1, x, log="x", type="o", lwd=2)) 
library(ggplot2)
library(splines)  #用于在 geom_smooth() 中添加拟合线，或者使用 geom_line() 替代 geom_smooth() 绘制普通折线
p <- ggplot(otu.cv, aes(otus, value)) +
  geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of Variables selected by RandomForest', y = 'Cross-validation error')
p
p + geom_vline(xintercept =5)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
importance_otu.select <- importance_otu[1:5, ]
importance_otu.select
write.table(importance_otu.select, 'randomforest.select.csv', sep = ',', col.names = NA, quote = FALSE)
importance=importance(x=rf)
# 从优化后的随机森林模型中提取各基因的重要性指标
importance=as.data.frame(importance)
# 将重要性指标转换为数据框格式，以便后续操作
importance$size=rownames(importance)
# 将行名（基因名）作为新列添加到数据框中
importance=importance[,c(2,1)]

names(importance)=c("Gene","importance")

#展示前20个基因的重要性
# 按照重要性指标从高到低对基因进行排序
af=importance[order(importance$importance,decreasing = T),]
af=af[c(1:5),]
af$Gene=rownames(af)
library(ggpubr)
p2=ggdotchart(af, x = "Gene", y = "importance",
              color = "importance", 
              sorting = "descending",                       
              add = "segments",                             
              add.params = list(color = "lightgray", size = 2), 
              dot.size = 6,                        
              font.label = list(color = "white", size = 9,
                                vjust = 0.5),               
              ggtheme = theme_bw()         ,               
              rotate=TRUE                                       )
p3=p2+ geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  gradient_color(palette =c(ggsci::pal_npg()(2)[2],ggsci::pal_npg()(2)[1])      ) 
p3
pdf(file="重要性基因2.pdf", width=6, height=6)
print(p3)
dev.off()
#############
##boruta
library(Boruta)
data<- read.csv("machinelearning.csv", header = T, sep=",",row.names =1)
data$group<- ifelse(data$group == "NC", 0,1)
set.seed(12345)
boruta_obj<-Boruta(group~.,
                   data=data,
                   doTrace=0,
                   ntree=500,
                   pValue=0.001)


print(boruta_obj)

ori_plot<-plot(boruta_obj,
               las=3,
               xlab = "",
               ylab = "Importance;Z-score",
               main="Variables selected by Boruta")

legend(8,9,
       c("Shadow","Rejected","Tendensive","Confirmed"),
       fill=c("blue","red","yellow","green"))
write.table(boruta_obj$finalDecision, 'boruta_select1.csv', sep = ',', col.names = NA, quote = FALSE)

#install.packages("xgboost")
#####xgbosst
library(ggplot2)
library(xgboost)
library(caret)

#读取输入文件
data=read.csv("RF.csv",sep=",",row.names = 1)
data$group<- ifelse(data$group == "NC", 0,1)
library(shapviz)
library(xgboost)
library(caret)
library(pROC)
library(tibble)
library(ROCit)
library(ggplot2)
library(data.table)
library(ggprism)
#data=read.csv("xgboost.csv",row.names = 1)
head(data)# 划分训练集和测试集
set.seed(123)#设
grid<-expand.grid(nrounds=c(75,100,150),
                  colsample_bytree=1,
                  min_child_weight=1,
                  eta=c(0.01,0.1,0.3),
                  gamma=c(0.5,0.25,0.1),
                  subsample=0.5,
                  max_depth=c(2,3,4))
cntrl<-trainControl(method="cv",
                    number=5,
                    verboseIter=F,
                    returnData=F,
                    returnResamp="final")
set.seed(123)
traindata=data
head(traindata)
train.xgb<-train(x=traindata[,1:38],y=as.factor(traindata[,39]),trControl=cntrl,tuneGrid=grid,method="xgbTree")
train.xgb#更改参数
model_xgboost=xgboost( data=as.matrix(traindata[,c(1:(ncol(traindata)-1))]), label=traindata$group, 
                       nrounds=100,
                       max_depth=3,
                       eta =0.3,
                       gamma=0.25,
                       colsample_bytree=1,
                       min_child_weight=1,
                       subsample=0.5,
                       objective="binary:logistic")
traindata$pred<-predict(model_xgboost,as.matrix(traindata[,c(1:(ncol(traindata)-1))]))

shap_xgboost=shapviz(model_xgboost,X_pred=as.matrix(traindata[,c(1:(ncol(traindata)-2))]))#瀑布图
sv_importance(shap_xgboost,kind="beeswarm")+theme_bw()
sv_importance(shap_xgboost)+theme_bw()
sv_force(shap_xgboost)
sv_importance(shap_xgboost, kind="beeswarm") + 
  theme_bw() +
  scale_size(range = c(2, 8)) 
sv_importance(shap_xgboost, kind="beeswarm") + 
  theme_bw() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))  # 上、右、下、左的边距均为1cm

shap_xgboost=shapviz(model_xgboost,X_pred=as.matrix(data[,c((ncol(data)-1))]))
sv_waterfall(shap_xgboost,row_id = 3, max_display = 30,fill_colors=c("#FF0000", "#0085FF"))
sv_force(shap, row_id = 2,  max_display = 30,fill_colors=c("#FF0000", "#0085FF"))
sv_importance(shap_xgboost, kind = "beeswarm",                  
              max_display = 48)
sv_importance(shap,fill = "#0085FF",              
              max_display = 48)


##################


###############
library(rms)
library(rms)
rt=read.csv("machinelearning.csv",row.names = 1)
rt$group<- ifelse(rt$group == "NC", 0,1)
# 使用Hmisc包的datadist函数计算变量的分布信息
ddist=datadist(rt)
options(datadist="ddist")
# 构建逻辑回归模型（logistic regression model）
# 'Type'是因变量（通常是二分类变量），后面是自变量（即多个基因）



lrmModel=lrm(group~CTSD+CYP1B1+METRNL, data=rt, x=T, y=T,maxit=1000)
# 创建Nomogram（列线图），显示逻辑回归模型的结果
# 'fun=plogis'指定使用逻辑函数将线性预测转换为概率
# 'fun.at'定义了概率轴上的刻度点
nomo=nomogram(lrmModel, fun=plogis,
              fun.at=c(0.01,0.1,0.5,0.99),
              lp=F, funlabel="Nomogram")


pdf("诺莫图train.pdf", width=8, height=6)
# 绘制Nomogram，调整坐标轴标签的字体大小
plot(nomo, cex.axis=0.8)

dev.off()

# 计算每个样本的Nomogram风险得分，即逻辑回归模型的预测概率
nomoRisk=predict(lrmModel, type="fitted")
# 将风险得分与原始数据表格合并，并添加表头
outTab=cbind(rt, Nomogram=nomoRisk)
outTab=rbind(id=colnames(outTab), outTab)
# 将结果保存为文本文件
write.table(outTab, file="模型评分TRAIN.txt", sep="\t", quote=F, col.names=F)

# 校准图，评估模型的预测概率与实际观察到的结果之间的吻合度
cali=calibrate(lrmModel, method="boot", B=1000)

# 绘制校准图，指定x轴为预测概率，y轴为实际概率
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=F)

{# 计算C-index，用于评估模型的区分能力
  cindex=rcorrcens(group~predict(lrmModel),data=rt)
  # 计算C-index的标准误
  se=cindex[,"SD"]/2
  # 格式化输出C-index及其95%置信区间
  c_index=sprintf("%.03f", cindex[,"C"])
  c_index.ci_low=sprintf("%.03f", cindex[,"C"]-(se*1.96))
  c_index.ci_high=sprintf("%.03f", cindex[,"C"]+(se*1.96))
  # 拼接C-index及其95%置信区间的字符串
  cindexLabel=paste0(c_index, " (95% CI: ", c_index.ci_low, "-", c_index.ci_high, ")")
  # 在图中添加C-index标签
  text(0.1, 0.85, "C-index:")
  text(0.28, 0.78, cindexLabel)
}

##ROC曲线
# 从合并后的表格中提取因变量，即Type列的数据
outTab=outTab[,c("group","CTSD","CYP1B1","METRNL","Nomogram")]
y=outTab[,"group"]
# 提供回归建模、校准图、C-index等函数
library(rmda)
library(Hmisc)
library(glmnet)
library(pROC)
library(ggsci)
# 使用Simpsons配色方案，生成曲线的颜色
bioCol=pal_simpsons(palette=c("springfield"), alpha=1)(length(2:ncol(outTab)))
aucText=c()  # 初始化保存AUC信息的向量
k=0  # 初始化计数器
# 遍历表格中从第2列开始的每个自变量，计算ROC曲线和AUC
for(x in colnames(outTab)[2:ncol(outTab)]){
  k=k+1
  # 计算ROC曲线，预测值转换为数值类型
  roc1=roc(y, as.numeric(outTab[,x]))     
  if(k==1){
    # 绘制第一条ROC曲线，不添加其他曲线
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="")
    # 保存AUC信息
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }else{
    # 绘制后续的ROC曲线，叠加到第一条曲线上
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE)
    # 保存AUC信息
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}

# 在图的右下角添加图例，显示每个变量的AUC值
legend("bottomright", aucText, lwd=1, bty="n", col=bioCol, cex=0.9)
rt=read.csv("machinelearning.csv",row.names = 1)
rt$group<- ifelse(rt$group == "Normal", 0,1)
library(pROC)
dfroc1<- roc(rt$group, rt$MXRA5)
dfroc2 <- roc(rt$group, rt$LRRC15)
dfroc3 <- roc(rt$group, rt$FOXC1)



library(pcutils)
library(ggstatsplot)
library(PMCMRplus) 
library(ggplot2) 
library(cowplot) 

set.seed(123)

#包含基因名和分组这两列
data <- read.csv("machinelearning.csv",sep = ",",check.names = F,stringsAsFactors = F,header = T)

{
  # 请根据您的实际数据结构进行调整
  plist = list()
  for (i in 1:9) {
    plist[[i]] = group_box(tab = data[c("CTSD","CYP1B1","METRNL")], group = "group", metadata = data, 
                           mode = i,p_value1 = "t.test") +#t.test，wilcox.test，anova，kruskal.test等方法进行比较
      ggtitle(paste0("mode", i)) +
      theme_classic() + 
      theme(legend.position = "none")
  }
  
  # 绘制箱式图
  plot_grid(plotlist = plist, ncol = 3)
}


# 单独展示第某一个模型
plot(plist[[4]])


###########
library(rmda)
data("dcaData")
#构架基础DCA模型
rt=read.table("模型评分TRAIN.txt",row.names = 1,header = T)
baseline.model1 <- decision_curve(group~ CTSD,
                                 family=binomial(link='logit'),
                                 data = rt,
                                 thresholds = seq(0, .4, by = .005),
                                 bootstraps = 10)
baseline.model2 <- decision_curve(group~ CYP1B1,
                                 family=binomial(link='logit'),
                                 data = rt,
                                 thresholds = seq(0, .4, by = .005),
                                 bootstraps = 10)
baseline.model3 <- decision_curve(group~ METRNL,
                                 family=binomial(link='logit'),
                                 data = rt,
                                 thresholds = seq(0, .4, by = .005),
                                 bootstraps = 10)

baseline.model4 <- decision_curve(group~ Nomogram,
                                  family=binomial(link='logit'),
                                  data = rt,
                                  thresholds = seq(0, .4, by = .005),
                                  bootstraps = 10)

plot_decision_curve(list(baseline.model1, baseline.model2,baseline.model3,baseline.model4),
                    curve.names = c("CSTD", "CYP1B1","METRNL","Nomogram"),
                    col = c("blue", "red","yellow","pink"),#设置线的颜色
                    lty = c(1,2),  #  设置线形
                    lwd = c(3,3, 2, 2),   #设置线宽，分别为baseling.model、full.model、All和None
                    confidence.intervals = F,
                    legend.position = "topright")


plot_decision_curve(
  list(baseline.model1, baseline.model2, baseline.model3, baseline.model4),
  curve.names = c("CSTD", "CYP1B1","METRNL","Nomogram"),
  col = c("blue", "red","yellow","pink"),
  lty = c(1,2,1,2),
  lwd = c(3,3,2,2),
  confidence.intervals = F,
  legend.position = "topright",
  legend.args = list(cex = 0.1)  # 0.8倍大小，越小越紧凑
)



par(cex = 0.4)  

# 2. 绘制决策曲线（去掉legend.args参数）
plot_decision_curve(
  list(baseline.model1, baseline.model2, baseline.model3, baseline.model4),
  curve.names = c("CSTD", "CYP1B1","METRNL","Nomogram"),
  col = c("blue", "red","yellow","pink"),
  lty = c(1,2,1,2),
  lwd = c(3,3,2,2),
  confidence.intervals = F,
  legend.position = "right"
)

# 3. 恢复全局参数默认值，避免影响后续绘图
par(cex = 1)


plot_decision_curve(
  list(baseline.model1, baseline.model2, baseline.model3, baseline.model4),
  curve.names = c("CSTD", "CYP1B1","METRNL","Nomogram"),
  col = c("blue", "orange", "green", "red"), # 修改为对比度高的颜色
  lty = c(1,2,1,2),
  lwd = c(3,3,2,2),
  confidence.intervals = F,
  legend.position = "right"
)


##绘图
plot(dfroc1,col="red",grid=c(0.2,0.2),grid.col=c("blue","yellow"))
plot(dfroc2,add=TRUE,col="blue")
plot(dfroc3,add=TRUE,col="black")
legend(0.4,0.5,c("MXRA5","LRRC15","FOXC1"),
       col=c("red","blue","black"),lty=1)



VEIN=read.csv("交集基因.csv",sep=",",header=T) 
colnames(VEIN)
rownames(VEIN)
exp=read.csv("exp23130.csv",sep=",",check.names=F,header=T,row.names = 1)
DEG_gene_exprvein1<-exp[VEIN$Intersection_Genes,]
DEG_gene_exprvein1<-na.omit(DEG_gene_exprvein1)
DEG_gene_exprvein1<-t(DEG_gene_exprvein1)
write.table(DEG_gene_exprvein1,file="validation.csv",sep=",",quote=T,col.names=T)






library(rms)
library(rms)
rt=read.csv("validation.csv",row.names = 1)
rt$group<- ifelse(rt$group == "NC", 0,1)
# 使用Hmisc包的datadist函数计算变量的分布信息
ddist=datadist(rt)
options(datadist="ddist")

# 构建逻辑回归模型（logistic regression model）
# 'Type'是因变量（通常是二分类变量），后面是自变量（即多个基因）
lrmModel=lrm(group~IRX3+CPAMD8+DMKN, data=rt, x=T, y=T,maxit=1000)
# 创建Nomogram（列线图），显示逻辑回归模型的结果
# 'fun=plogis'指定使用逻辑函数将线性预测转换为概率
# 'fun.at'定义了概率轴上的刻度点
nomo=nomogram(lrmModel, fun=plogis,
              fun.at=c(0.01,0.1,0.5,0.99),
              lp=F, funlabel="Nomogram")


pdf("诺莫图validation.pdf", width=8, height=6)
# 绘制Nomogram，调整坐标轴标签的字体大小
plot(nomo, cex.axis=0.8)

dev.off()

# 计算每个样本的Nomogram风险得分，即逻辑回归模型的预测概率
nomoRisk=predict(lrmModel, type="fitted")
# 将风险得分与原始数据表格合并，并添加表头
outTab=cbind(rt, Nomogram=nomoRisk)
outTab=rbind(id=colnames(outTab), outTab)
# 将结果保存为文本文件
write.table(outTab, file="模型评分validation.txt", sep="\t", quote=F, col.names=F)

# 校准图，评估模型的预测概率与实际观察到的结果之间的吻合度
cali=calibrate(lrmModel, method="boot", B=1000)

# 绘制校准图，指定x轴为预测概率，y轴为实际概率
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=F)

{# 计算C-index，用于评估模型的区分能力
  cindex=rcorrcens(group~predict(lrmModel),data=rt)
  # 计算C-index的标准误
  se=cindex[,"SD"]/2
  # 格式化输出C-index及其95%置信区间
  c_index=sprintf("%.03f", cindex[,"C"])
  c_index.ci_low=sprintf("%.03f", cindex[,"C"]-(se*1.96))
  c_index.ci_high=sprintf("%.03f", cindex[,"C"]+(se*1.96))
  # 拼接C-index及其95%置信区间的字符串
  cindexLabel=paste0(c_index, " (95% CI: ", c_index.ci_low, "-", c_index.ci_high, ")")
  # 在图中添加C-index标签
  text(0.1, 0.85, "C-index:")
  text(0.28, 0.78, cindexLabel)
}

##ROC曲线


DMKN


# 从合并后的表格中提取因变量，即Type列的数据
outTab=outTab[,c("group","DMKN","IRX3","CPAMD8")]
y=outTab[,"group"]
# 提供回归建模、校准图、C-index等函数
library(rmda)
library(Hmisc)
library(glmnet)
library(pROC)
library(ggsci)
# 使用Simpsons配色方案，生成曲线的颜色
bioCol=pal_simpsons(palette=c("springfield"), alpha=1)(length(2:ncol(outTab)))
aucText=c()  # 初始化保存AUC信息的向量
k=0  # 初始化计数器
# 遍历表格中从第2列开始的每个自变量，计算ROC曲线和AUC
for(x in colnames(outTab)[2:ncol(outTab)]){
  k=k+1
  # 计算ROC曲线，预测值转换为数值类型
  roc1=roc(y, as.numeric(outTab[,x]))     
  if(k==1){
    # 绘制第一条ROC曲线，不添加其他曲线
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="")
    # 保存AUC信息
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }else{
    # 绘制后续的ROC曲线，叠加到第一条曲线上
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE)
    # 保存AUC信息
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}

# 在图的右下角添加图例，显示每个变量的AUC值
legend("bottomright", aucText, lwd=1, bty="n", col=bioCol, cex=0.9)
rt=read.csv("machinelearning.csv",row.names = 1)
rt$group<- ifelse(rt$group == "Normal", 0,1)
library(pROC)
dfroc1<- roc(rt$group, rt$MXRA5)
dfroc2 <- roc(rt$group, rt$LRRC15)
dfroc3 <- roc(rt$group, rt$FOXC1)

##绘图
plot(dfroc1,col="red",grid=c(0.2,0.2),grid.col=c("blue","yellow"))
plot(dfroc2,add=TRUE,col="blue")
plot(dfroc3,add=TRUE,col="black")
legend(0.4,0.5,c("MXRA5","LRRC15","FOXC1"),
       col=c("red","blue","black"),lty=1)






##############

library(pcutils)
library(ggstatsplot)
library(PMCMRplus) 
library(ggplot2) 
library(cowplot) 

set.seed(123)

#包含基因名和分组这两列
data <- read.csv("machinelearning.csv",sep = ",",check.names = F,stringsAsFactors = F,header = T)

{
  # 请根据您的实际数据结构进行调整
  plist = list()
  for (i in 1:9) {
    plist[[i]] = group_box(tab = data[c("DMKN","IRX3","CPAMD8")], group = "group", metadata = data, 
                           mode = i,p_value1 = "t.test") +#t.test，wilcox.test，anova，kruskal.test等方法进行比较
      ggtitle(paste0("mode", i)) +
      theme_classic() + 
      theme(legend.position = "none")
  }
  
  # 绘制箱式图
  plot_grid(plotlist = plist, ncol = 3)
}


# 单独展示第某一个模型
plot(plist[[4]])
















library(Biobase)
library(GEOquery)

#Download GPL file, put it in the current directory, and load it:
gpl <- getGEO('GPL11154', destdir=".")
colnames(Table(gpl)) ## [1] 41108 17
## 重点就是要花时间来摸索这个返回值
head(Table(gpl)[,c(1,10,13)]) ## you need to check this , which column do you need 
probe2symbol=Table(gpl)[,c(1,13)]
library(hgu133a.db)
ids=toTable(hgu133aSYMBOL)
head(ids)
## 或者
platformDB='hugene10sttranscriptcluster.db'
library(platformDB, character.only=TRUE)
probeset <- featureNames(GSE62832[[1]])


########finn-b-M13_LOWBACKPAIN


######
library(tidyverse)
library(GEOquery)
library(tidyverse)
library(GEOquery)
library(limma) 
library(affy)
library(stringr)
gset <- getGEO('GSE23130', destdir=".",
               AnnotGPL = F,     
               getGPL = F)  
exp<-exprs(gset[[1]])
#exp<-log2(exp+1)
cli<-pData(gset[[1]])	
GPL<-read.table("GPL1352-9802 (1).txt",sep = "\t",header = T)
GPL<-GPL[,c(1,2)]
#GPL<-fData(gset[[1]])
gpl<-GPL[,c(1,2)]        
gpl$Gene.Symbol<-data.frame(sapply(gpl$Gene.Symbol,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]
exp<-as.data.frame(exp)
exp$ID<-rownames(exp)
exp_symbol<-merge(exp,gpl,by="ID")
exp_symbol<-na.omit(exp_symbol)
table(duplicated(exp_symbol$Gene.Symbol))
exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$Gene.Symbol)
write.table(exp_unique,"exp23130.csv",sep=",")
write.table(cli,"group23130.csv",sep=",")

exp1<-read.csv("exp70362.csv",row.names = 1)
exp2<-read.csv("exp56081.csv",row.names = 1)
table(rownames(exp1) %in% rownames(exp2))
length(intersect(rownames(exp1),rownames(exp2)))
exp1 <- exp1[intersect(rownames(exp1),rownames(exp2)),]
exp2 <- exp2[intersect(rownames(exp1),rownames(exp2)),]
boxplot(exp1)
boxplot(exp2)


#(4)合并表达矩阵
# exp2的第三个样本有些异常，可以去掉或者用normalizeBetweenArrays标准化，把它拉回正常水平。
exp = cbind(exp1,exp2)
#exp=exp[,-(24:33)]
#exp=exp[,-(34:43)]
boxplot(exp)


#处理批次效应
library(limma)
#?removeBatchEffect()
batch <- c(rep("exp1",24),rep("exp2",10))
exp3 <- removeBatchEffect(exp, batch)
par(mfrow=c(1,2))  # 展示的图片为一行两列
boxplot(as.data.frame(exp),main="Original")
boxplot(as.data.frame(exp3),main="Batch corrected")
write.table(exp3,"expmerge.csv",sep=",",quote = F,row.names =T,col.names = T)

########
dir.create("./DEG")
setwd("./DEG")
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

expFile="exp70362.csv"          #?????ļ?
#??ȡ?????ļ????????????ļ?????
expr_data=read.csv(expFile,sep=",",header=T,row.names = 1,check.names=F)
group<-read.csv("group70362.csv",header = T,row.names = 1,sep = ",")
design <- model.matrix(~0+factor(group$group))
colnames(design) <- levels(factor(group$group))
rownames(design) <- colnames(expr_data)
contrast.matrix <- makeContrasts(LDH-NC,levels = design) 
fit <- lmFit(expr_data,design) 
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
DEG$regulate <- ifelse(DEG$P.Value> 0.05, "unchanged",
                       ifelse(DEG$logFC > 0.5, "up-regulated",
                              ifelse(DEG$logFC < -0.5, "down-regulated", "unchanged")))
write.table(DEG,"DEG.csv",row.names=T,col.names=T,sep=",")

DE_1_0.05 <- DEG[DEG$P.Value<0.05&abs(DEG$logFC)>0.5,]
upGene_1_0.05 <- DE_1_0.05[DE_1_0.05$regulate == "up-regulated",]
downGene_1_0.05 <- DE_1_0.05[DE_1_0.05$regulate == "down-regulated",]
write.csv(upGene_1_0.05,"upGene_1_005.csv")
write.csv(downGene_1_0.05,"downGene_1_005.csv")
#volcano plot of DEGs
library(dplyr) # 
library(gt) 
data<-read.csv("DEG.csv", header = TRUE, sep = ",", dec = ".", quote = "\"", fill = TRUE, comment.char = "")
data<- data%>%
  
  mutate(expression = case_when(logFC>= 1 &P.Value< 0.05 ~ "up", 
                                
                                logFC<= -1 &P.Value< 0.05 ~ "down", 
                                
                                TRUE ~ "Unchanged")) 
data<-na.omit(data)
library(ggplot2)
library(ggrepel)
volc_plot<- ggplot(data, aes(logFC, -log10(P.Value))) +xlim(-5,5)+ylim(0,7.5)+# 将p值进行-log10转化
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1.2,1.2), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = expression),
             size =2.5, 
             alpha =0.5) +
  theme_bw(base_size = 12)+
  ggsci::scale_color_jama() +
  theme(panel.grid = element_blank(),
        legend.position = 'right')

volc_plot


expFile="exp88837.csv"          #?????ļ?
#??ȡ?????ļ????????????ļ?????
expr_data=read.csv(expFile,sep=",",header=T,row.names = 1,check.names=F)

group <- read.csv(
  file = "group.csv",  # 若修复后另存为新文件，需改名为"group88837_fixed.csv"
  header = TRUE,            # 保留表头
  row.names = 1,            # 用第1列作为行名（样本名）
  sep = ",",                # 分隔符为逗号（CSV默认）
  skipNul = TRUE,           # 跳过嵌入空值
  fill = TRUE               # 填充不完整的行
)

design <- model.matrix(~0+factor(group$group))
colnames(design) <- levels(factor(group$group))
rownames(design) <- colnames(expr_data)
contrast.matrix <- makeContrasts(Obesity-NObesity,levels = design) 
fit <- lmFit(expr_data,design) 
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
DEG$regulate <- ifelse(DEG$adj.P.Val> 0.05, "unchanged",
                       ifelse(DEG$logFC > 1, "up-regulated",
                              ifelse(DEG$logFC < -1, "down-regulated", "unchanged")))
write.table(DEG,"DEGobesity.csv",row.names=T,col.names=T,sep=",")

#volcano plot of DEGs
library(dplyr) # 
library(gt) 
data<-read.csv("DEGobesity.csv", header = TRUE, sep = ",", dec = ".", quote = "\"", fill = TRUE, comment.char = "")
data<- data%>%
  
  mutate(expression = case_when(logFC>= 1 &P.Value< 0.05 ~ "up", 
                                
                                logFC<= -1 &P.Value< 0.05 ~ "down", 
                                
                                TRUE ~ "Unchanged")) 
data<-na.omit(data)
library(ggplot2)
library(ggrepel)
volc_plot<- ggplot(data, aes(logFC, -log10(P.Value))) +xlim(-5,5)+ylim(0,7.5)+# 将p值进行-log10转化
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1.2,1.2), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = expression),
             size =2.5, 
             alpha =0.5) +
  theme_bw(base_size = 12)+
  ggsci::scale_color_jama() +
  theme(panel.grid = element_blank(),
        legend.position = 'right')

volc_plot






#

############
library(glmnet)
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)
rt=read.csv("exp70362.csv",row.names = 1,header = T)
gene=read.csv("gene.csv",header = T)
rt=rt[gene$GENE,]
rt=t(rt)
rt=as.data.frame(rt)
group<-read.csv("group70362.csv",row.names = 1)
merged_data <- cbind(group,rt)
#write.csv(rt,"rt.csv")
rt=merged_data
inTrain<-createDataPartition(y=rt[,1], p=0.7, list=F) 
#train=merged_data
train<-rt[inTrain,]
test<-rt[-inTrain,]
set.seed(123) 
dir.create("../9machinelearning/")
setwd("../9machinelearning/")
colnames(train)[1] <- "Type"
colnames(test)[1] <- "Type"
all <- rbind(train,test)

control=trainControl(method="repeatedcv", number=5, savePredictions="final",classProbs = TRUE,summaryFunction = twoClassSummary)
# ranger
rfGrid <- expand.grid(mtry =c(1,2,3,4,5,6), 
                      splitrule = c( "gini", "extratrees", "hellinger"),
                      min.node.size = seq(1,15,2))
mod_rf <- caret::train(Type~.,train,method = "ranger",trControl = control,
                       num.trees = 1000,tuneGrid = rfGrid,
                       metric = "ROC")
## lasso
grid_glmnet <- expand.grid(
  alpha  = seq(0, 1, 0.1),                  # 0=ridge, 1=lasso
  lambda = 10^seq(-4, 1, length = 50)
)
mod_lasso <- train(
  Type ~ ., data = train,
  method = "glmnet",
  trControl = control, metric = "ROC",
  tuneGrid = grid_glmnet
)
#SVM
grid_svm <- expand.grid(
  sigma = 2^seq(-15, 3, 2),                
  C     = 2^seq(-5, 15, 2)                
)
mod_svm <- train(
  Type ~ ., data = train,
  method = "svmRadial",
  trControl = control, metric = "ROC",
  tuneGrid = grid_svm
)
#XGB
grid_xgb <- expand.grid(
  nrounds          = seq(200, 1000, 200),
  max_depth        = c(3,4,5,6),
  eta              = c(0.01,0.05, 0.1),
  gamma            = c(1,3,5),
  colsample_bytree = c(0.5,0.6,0.6, 0.8),
  min_child_weight = c(1,2,3),
  subsample        = c(0.5,1),
  rate_drop = c(0.5),
  skip_drop = c(0.05)
)
mod_xgb <- train(
  Type ~ ., data = train,
  method = "xgbDART",
  trControl = control, metric = "ROC",tuneGrid = grid_xgb
)
#GLM
mod_glm=train(Type ~., data = train, method = "glm", family="binomial", trControl=control,metric = "ROC")
#GBM
grid_gbm <- expand.grid(
  n.trees          = seq(100, 800, 200),
  interaction.depth= c(1, 2, 3),
  shrinkage        = c(0.01, 0.05, 0.1),
  n.minobsinnode   = c(5, 10, 20)
)
mod_gbm=train(Type ~., data = train, method = "gbm", trControl=control,metric = "ROC",tuneGrid = grid_gbm)
######KNN
grid_knn <- expand.grid(
  k = seq(3, 41, 2)  
)
mod_knn=train(Type ~., data = train, method = "knn", trControl=control,metric = "ROC",tuneGrid = grid_knn)

######NNET
grid_nnet <- expand.grid(
  size  = c(1, 3, 5, 7),     
  decay = 10^seq(-4, -1, 1)   
)
mod_nnet=train(Type ~., data = train, method = "nnet", trControl=control,metric = "ROC",tuneGrid = grid_nnet)
# dt
grid_dt <- expand.grid(
  cp = seq(0.001, 0.1, 0.005)  
)
mod_dt=train(Type ~., data = train, method = "rpart", trControl=control,metric = "ROC",tuneGrid = grid_dt)
save(control,mod_dt,,mod_glm,mod_knn,mod_lasso,mod_nnet,mod_rf,mod_svm,file = "7model_tuneParameter.Rdata")

#save(control,mod_dt,mod_gbm,mod_glm,mod_knn,mod_lasso,mod_nnet,mod_rf,mod_svm,mod_xgb,file = "9model_tuneParameter.Rdata")


p_fun=function(object, newdata){
  predict(object, newdata=newdata, type="prob")[,2]
}
yTest=as.numeric(ifelse(test$Type == "LDH", 1, 0))
yTrain = as.numeric(ifelse(train$Type == "LDH",1,0))
yall=as.numeric(ifelse(all$Type == "LDH",1,0))

dt <- data.frame(row.names = c("RF","SVM","GLM","KNN","NNET","LASSO","DT"))

# ROC(all)
pred1=predict(mod_rf, newdata=all, type="prob")
pred2=predict(mod_svm, newdata=all, type="prob")
#pred3=predict(mod_xgb, newdata=all, type="prob")
pred4=predict(mod_glm, newdata=all, type="prob")
#pred5=predict(mod_gbm, newdata=all, type="prob")
pred6=predict(mod_knn, newdata=all, type="prob")
pred7=predict(mod_nnet, newdata=all, type="prob")
pred8=predict(mod_lasso, newdata=all, type="prob")
pred9=predict(mod_dt, newdata=all, type="prob")

roc1=roc(yall, as.numeric(pred1[,2]))
roc2=roc(yall, as.numeric(pred2[,2]))
#roc3=roc(yall, as.numeric(pred3[,2]))
roc4=roc(yall, as.numeric(pred4[,2]))
#roc5=roc(yall, as.numeric(pred5[,2]))
roc6=roc(yall, as.numeric(pred6[,2]))
roc7=roc(yall, as.numeric(pred7[,2]))
roc8=roc(yall, as.numeric(pred8[,2]))
roc9=roc(yall, as.numeric(pred9[,2]))
dt$All <- c(roc1$auc,roc2$auc,roc4$auc,roc6$auc,roc7$auc,roc8$auc,roc9$auc)


pdf(file="ALLROC.pdf", width=5, height=5)
plot(roc1, print.auc=F, legacy.axes=T, main="", col="chocolate")
plot(roc2, print.auc=F, legacy.axes=T, main="", col="aquamarine3", add=T)
plot(roc4, print.auc=F, legacy.axes=T, main="", col="bisque3", add=T)
plot(roc6, print.auc=F, legacy.axes=T, main="", col="burlywood", add=T)
plot(roc7, print.auc=F, legacy.axes=T, main="", col="darkgoldenrod3", add=T)
plot(roc8, print.auc=F, legacy.axes=T, main="", col="darkolivegreen3", add=T)
plot(roc9, print.auc=F, legacy.axes=T, main="", col="dodgerblue3", add=T)

legend('bottomright',
       c(paste0('RF: ',sprintf("%.03f",roc1$auc)),
         paste0('SVM: ',sprintf("%.03f",roc2$auc)),
         paste0('GLM: ',sprintf("%.03f",roc4$auc)),
         paste0('KNN: ',sprintf("%.03f",roc6$auc)),
         paste0('NNET: ',sprintf("%.03f",roc7$auc)),
         paste0('LASSO: ',sprintf("%.03f",roc8$auc)),
         paste0('LASSO: ',sprintf("%.03f",roc9$auc))),
       
       
       col=c("chocolate","aquamarine3","bisque3",
             "burlywood","darkgoldenrod3","darkolivegreen3",
             "dodgerblue3"), lwd=2, bty = 'n')
dev.off()










#dt$All <- c(roc1$auc,roc2$auc,roc3$auc,roc4$auc,roc5$auc,roc6$auc,roc7$auc,roc8$auc,roc9$auc)

# ROC(train)
pred1=predict(mod_rf, newdata=train, type="prob")
pred2=predict(mod_svm, newdata=train, type="prob")
#pred3=predict(mod_xgb, newdata=train, type="prob")
pred4=predict(mod_glm, newdata=train, type="prob")
#pred5=predict(mod_gbm, newdata=train, type="prob")
pred6=predict(mod_knn, newdata=train, type="prob")
pred7=predict(mod_nnet, newdata=train, type="prob")
pred8=predict(mod_lasso, newdata=train, type="prob")
pred9=predict(mod_dt, newdata=train, type="prob")

roc1=roc(yTrain, as.numeric(pred1[,2]))
roc2=roc(yTrain, as.numeric(pred2[,2]))
#roc3=roc(yTrain, as.numeric(pred3[,2]))
roc4=roc(yTrain, as.numeric(pred4[,2]))
#roc5=roc(yTrain, as.numeric(pred5[,2]))
roc6=roc(yTrain, as.numeric(pred6[,2]))
roc7=roc(yTrain, as.numeric(pred7[,2]))
roc8=roc(yTrain, as.numeric(pred8[,2]))
roc9=roc(yTrain, as.numeric(pred9[,2]))

pdf(file="TRAINROC.pdf", width=5, height=5)
plot(roc1, print.auc=F, legacy.axes=T, main="", col="chocolate")
plot(roc2, print.auc=F, legacy.axes=T, main="", col="aquamarine3", add=T)
plot(roc4, print.auc=F, legacy.axes=T, main="", col="bisque3", add=T)
plot(roc6, print.auc=F, legacy.axes=T, main="", col="burlywood", add=T)
plot(roc7, print.auc=F, legacy.axes=T, main="", col="darkgoldenrod3", add=T)
plot(roc8, print.auc=F, legacy.axes=T, main="", col="darkolivegreen3", add=T)
plot(roc9, print.auc=F, legacy.axes=T, main="", col="dodgerblue3", add=T)

legend('bottomright',
       c(paste0('RF: ',sprintf("%.03f",roc1$auc)),
         paste0('SVM: ',sprintf("%.03f",roc2$auc)),
         paste0('GLM: ',sprintf("%.03f",roc4$auc)),
         paste0('KNN: ',sprintf("%.03f",roc6$auc)),
         paste0('NNET: ',sprintf("%.03f",roc7$auc)),
         paste0('LASSO: ',sprintf("%.03f",roc8$auc)),
         paste0('LASSO: ',sprintf("%.03f",roc9$auc))),
       
       
       col=c("chocolate","aquamarine3","bisque3",
             "burlywood","darkgoldenrod3","darkolivegreen3",
             "dodgerblue3"), lwd=2, bty = 'n')
dev.off()


dt$Train <- c(roc1$auc,roc2$auc,roc4$auc,roc6$auc,roc7$auc,roc8$auc,roc9$auc)
#dt$Train <- c(roc1$auc,roc2$auc,roc3$auc,roc4$auc,roc5$auc,roc6$auc,roc7$auc,roc8$auc,roc9$auc)

# ROC(test)
pred1=predict(mod_rf, newdata=test, type="prob")
pred2=predict(mod_svm, newdata=test, type="prob")
#pred3=predict(mod_xgb, newdata=test, type="prob")
pred4=predict(mod_glm, newdata=test, type="prob")
#pred5=predict(mod_gbm, newdata=test, type="prob")
pred6=predict(mod_knn, newdata=test, type="prob")
pred7=predict(mod_nnet, newdata=test, type="prob")
pred8=predict(mod_lasso, newdata=test, type="prob")
pred9=predict(mod_dt, newdata=test, type="prob")

roc1=roc(yTest, as.numeric(pred1[,2]))
roc2=roc(yTest, as.numeric(pred2[,2]))
#roc3=roc(yTest, as.numeric(pred3[,2]))
roc4=roc(yTest, as.numeric(pred4[,2]))
#roc5=roc(yTest, as.numeric(pred5[,2]))
roc6=roc(yTest, as.numeric(pred6[,2]))
roc7=roc(yTest, as.numeric(pred7[,2]))
roc8=roc(yTest, as.numeric(pred8[,2]))
roc9=roc(yTest, as.numeric(pred9[,2]))
pdf(file="TESTROC.pdf", width=5, height=5)
plot(roc1, print.auc=F, legacy.axes=T, main="", col="chocolate")
plot(roc2, print.auc=F, legacy.axes=T, main="", col="aquamarine3", add=T)
plot(roc4, print.auc=F, legacy.axes=T, main="", col="bisque3", add=T)
plot(roc6, print.auc=F, legacy.axes=T, main="", col="burlywood", add=T)
plot(roc7, print.auc=F, legacy.axes=T, main="", col="darkgoldenrod3", add=T)
plot(roc8, print.auc=F, legacy.axes=T, main="", col="darkolivegreen3", add=T)
plot(roc9, print.auc=F, legacy.axes=T, main="", col="dodgerblue3", add=T)

legend('bottomright',
       c(paste0('RF: ',sprintf("%.03f",roc1$auc)),
         paste0('SVM: ',sprintf("%.03f",roc2$auc)),
         paste0('GLM: ',sprintf("%.03f",roc4$auc)),
         paste0('KNN: ',sprintf("%.03f",roc6$auc)),
         paste0('NNET: ',sprintf("%.03f",roc7$auc)),
         paste0('LASSO: ',sprintf("%.03f",roc8$auc)),
         paste0('LASSO: ',sprintf("%.03f",roc9$auc))),
       
       
       col=c("chocolate","aquamarine3","bisque3",
             "burlywood","darkgoldenrod3","darkolivegreen3",
             "dodgerblue3"), lwd=2, bty = 'n')
dev.off()


dt$Test <- c(roc1$auc,roc2$auc,roc4$auc,roc6$auc,roc7$auc,roc8$auc,roc9$auc)
#dt$Test <- c(roc1$auc,roc2$auc,roc3$auc,roc4$auc,roc5$auc,roc6$auc,roc7$auc,roc8$auc,roc9$auc)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
dt$Mean_AUC <- rowMeans(dt)
dt <- dt[order(dt$Mean_AUC, decreasing = T),]
write.csv(dt,"model_AUC_all_train_test.csv",quote = F,row.names = T)
col_ha <- HeatmapAnnotation(
  which = "col",  # 标注在热图的列侧
  Cohort = c("All", "Train", "Test", "Mean_AUC"),  # 标注的分组内容
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),  # 标注名称的字体设置
  annotation_name_side = "left",  # 标注名称显示在左侧
  # 修正颜色映射：将 #456 改为 #445566（或其他6位合法颜色）
  col = list(
    Cohort = c(
      "All" = "#F39B7FFF",    # 8位格式（含透明度，合法）
      "Train" = "#8DA0CB",    # 6位格式（合法）
      "Test" = "#66C2A5",     # 6位格式（合法）
      "Mean_AUC" = "#445566"  # 修正：3位 #456 → 6位 #445566（合法）
    )
  ),
  show_legend = F,  # 不显示该标注的图例（若后续需显示，可改为 T）
  # 图例参数（虽 show_legend=F，但保留正确格式，便于后续调整）
  annotation_legend_param = list(
    Cohort = list(
      title = "Cohort",
      title_position = "topleft",
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_rot = 0,
      legend_height = unit(1, "cm"),
      legend_width = unit(5, "mm"),
      labels_gp = gpar(fontsize = 8, fontface = "bold")
    )
  )
)

pdf("AUC_Heatmap1.pdf", width = 7, height = 5)
heatmap_data <- as.matrix(dt[, c("All", "Train", "Test", "Mean_AUC")])

# 2. 修正：colorRamp2 断点与颜色匹配（以3断点3颜色为例）
color_map <- colorRamp2(
  breaks = c(0.7, 0.85, 1),  # 3个断点
  colors = c("#DDF1FC", "#FFF7AC", "#ECB477")  # 3个颜色（无#efa，避免非法颜色）
)

# 3. 重新构建热图（代入修正后的数据和颜色映射）
heatmap <- Heatmap(
  heatmap_data,  # 修正：矩阵格式的数据
  name = "AUC", 
  heatmap_legend_param = list(
    title = "AUC",
    title_position = "topleft", 
    labels_rot = 0,
    legend_height = unit(3, "cm"),
    legend_width = unit(5, "mm"),
    labels_gp = gpar(fontsize = 8, fontface = "bold")
  ),
  border = TRUE,
  column_split = factor(
    c("All", "Train", "Test", "Mean_AUC"),
    levels = c("All", "Train", "Test", "Mean_AUC")
  ),
  column_gap = unit(1.5, "mm"),
  show_column_names = FALSE,
  show_row_names = TRUE,
  col = color_map,  # 修正：合法的颜色映射
  column_title = "", 
  column_title_side = "top",
  row_title_side = "left",
  column_title_gp = gpar(fontsize = 12, fontface = "bold", col = "black"),
  row_title_gp = gpar(fontsize = 15, fontface = "bold", col = "black"),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  row_title_rot = 90,
  column_order = c("All", "Train", "Test", "Mean_AUC"),
  show_row_dend = FALSE,
  top_annotation = col_ha  # 需确保col_ha已修正（如Mean_AUC颜色为#445566）
)
draw(heatmap)
dev.off()

imp_tbl <- varImp(mod_rf, scale = TRUE)$importance
imp_tbl$Feature <- rownames(imp_tbl)
imp_tbl <- imp_tbl[order(-imp_tbl$Overall), c("Feature","Overall")]
write.csv(imp_tbl, "RF_importance_train.csv", row.names = FALSE)










library(kernelshap)
train_data = train
pred_fun <- function(object, newdata) {
  probs <- predict(object, newdata = newdata, type = "prob")
  as.numeric(probs[,"LDH"])
}

# 从 train_data 中抽取 100 行作为背景数据（bg_X），仅保留特征列（与解释数据一致）
bg_X <- train_data[
  sample(nrow(train_data), size = min(100, nrow(train_data))),  # 抽100行，若train_data不足100行则抽全部
  -1  # 排除标签列，仅保留特征列（与 X = train_data[, -1] 一致）
]

# 验证 bg_X 样本量（确保≥10行）
cat("背景数据（bg_X）样本量：", nrow(bg_X), "行\n")  # 输出应≥10

shap_values <- kernelshap(
  object = mod_lasso,          # 你的随机森林模型（mod_rf）
  X = train_data[, -1],     # 待解释的特征数据（解释数据，可少量）
  bg_X = bg_X,              # 手动指定的背景数据（解决“X太小”错误的核心）
  pred_fun = pred_fun       # 适配 mod_rf 的预测函数（如之前定义的针对ranger/randomForest的函数）
)


#shap_values <- kernelshap(mod_rf,train_data[, -1],pred_fun = pred_fun)
library(shapviz)
shap_vis <- shapviz(shap_values, train_data[, -1])
feature_importance <- colMeans(abs(shap_values$S))
sorted_features <- names(sort(feature_importance, decreasing = TRUE))
library(canvasXpress)
visualization_theme <- theme_minimal() + theme(plot.title = element_text(face ="bold", size =14),  
                                               axis.title = element_text(size =12))








sv_dependence(shap_vis, sorted_features) + visualization_theme 
plots <- lapply(sorted_features, function(feat) {
  sv_dependence(shap_vis, feat) + visualization_theme
})
library(cowplot)
p_waterfall <- sv_waterfall(
  shap_vis, 
  row_id = 5# 第4个样本的瀑布图
) + 
  labs(x = "SHAP value (Prediction)") +  # x轴标签
  theme_bw() +  # 建议加白色背景（与你之前的风格统一，可选）
  theme(
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.x = element_text(size = 12, color = "black", face = "bold")
  )
p_waterfall
# 2. 保存瀑布图（用 ggsave 指定宽6、高2.6）
ggsave(
  filename = "2_SHAP_waterfall_trainRow4.pdf",
  plot = p_waterfall,
  width = 5,
  height = 3.5,
  dpi = 300,
  device = "pdf"
)

# （可选）预览瀑布图
print(p_waterfall)


visualization_theme <- theme_minimal() + theme(plot.title = element_text(face ="bold", size =14),  
                                               axis.title = element_text(size =12))
sv_importance(shap_vis, kind="bar", show_numbers=TRUE) +
  labs(x = "Mean |importance|", y = NULL)+theme_bw()+theme(axis.text.y = element_text(color = "black",size = 12),axis.title.x = element_text(face = "bold",size = 14))
ggsave("2_SHAP_Feature_Importance_Barplot.pdf", width=5.6, height=3)
sv_importance(shap_vis, kind="bee", show_numbers=TRUE) +  
  visualization_theme +  labs( x = "SHAP value (impact on model output)",y = NULL)+
  theme_bw()+theme(axis.text.y = element_text(color = "black",size = 11),axis.title.x = element_text(size = 12,color = "black",face = "bold"))
ggsave("2_SHAP_BeeSwarm_Plot.pdf", width=6, height=2.6)

sv_dependence(shap_vis, sorted_features) + visualization_theme 
plots <- lapply(sorted_features, function(feat) {
  sv_dependence(shap_vis, feat) + visualization_theme
})
library(cowplot)
plot_grid(plotlist = plots, ncol = 3)
ggsave("2_SHAP_feature_dependence.pdf",width = 9,height = 4)

sv_waterfall(shap_vis, row_id=4) + 
  labs(x = "SHAP value (Prediction)")+
  theme(axis.text.y = element_text(color = "black",size = 11),axis.title.x = element_text(size = 12,color = "black",face = "bold"))
ggsave("2_SHAP_waterfall_trainRow4.pdf",width = 6,height = 2.6)


library(ggplot2)
library(shapviz)

# 1. 提取 SHAP 重要性数据（与棒棒糖图逻辑一致）
importance_data <- data.frame(
  Feature = colnames(shap_vis$S),  # 特征名
  Importance = colMeans(abs(shap_vis$S))  # 重要性（SHAP绝对值均值）
) %>%
  arrange(desc(Importance))  # 按重要性降序排序（y轴从高到低）

# 2. 绘制棒棒糖图 + 0.0045 垂直虚线 + 数值标注
p_lollipop_with_line <- ggplot(importance_data, aes(
  x = Importance, 
  y = reorder(Feature, Importance)  # 按重要性重排y轴特征
)) +
  # 1. 画“棒棒糖杆”（从x=0到x=重要性值）
  geom_segment(
    aes(xend = 0, yend = Feature),
    color = "#2E86AB",  # 杆的颜色（蓝色）
    linewidth = 1.2     # 杆的粗细
  ) +
  # 2. 画“棒棒糖头”（顶端圆点）
  geom_point(
    size = 4,
    color = "#A23B72",  # 点的边框色（紫红）
    fill = "#F18F01",   # 点的填充色（橙色）
    shape = 21          # 带填充的圆形
  ) +
  # 3. （核心）添加 X=0.0045 垂直虚线
  geom_vline(
    xintercept = 0.0045,  # 虚线位置：X轴0.0045
    color = "#666666",    # 虚线颜色（深灰，不抢焦点）
    linetype = "dashed",  # 线型：虚线
    linewidth = 0.8       # 线宽（适配图尺寸）
  ) +
  # 4. （核心）标注“0.0045”数值（在虚线旁）
  annotate(
    "text",                # 注释类型：文本
    x = 0.0045,            # 文本X位置=虚线位置
    y = nrow(importance_data) + 0.8,  # 文本Y位置：最上方特征的上方（避免遮挡）
    label = "0.0046",      # 标注内容
    color = "#666666",     # 文本颜色=虚线颜色（统一风格）
    size = 3.8,            # 文本大小（适配图尺寸）
    angle = 90,            # 文本旋转90度（垂直排列，节省横向空间）
    fontface = "bold"      # 文本加粗（更醒目）
  ) +
  # 5. （可选）添加重要性数值标签（在棒棒糖头右侧）
  geom_text(
    aes(label = round(Importance, 4)),  # 保留4位小数，可调整
    hjust = -0.5,  # 标签在点的右侧（避免重叠）
    size = 3.5,
    color = "black"
  ) +
  # 6. 轴标签与主题（与之前风格统一）
  labs(
    x = "Mean |SHAP value|",  # X轴标题（重要性指标）
    y = NULL                   # 隐藏Y轴标题（特征名已清晰）
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(color = "black", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()  # 隐藏次要网格线，更简洁
  ) +
  # 7. 调整X轴范围（避免数值标签/0.0045标注超出图外）
  xlim(0, max(importance_data$Importance) * 1.3)  # 右边界留30%空白

# 3. 保存图片（与之前尺寸一致）
ggsave(
  "2_SHAP_Feature_Importance_Lollipop_with_line.pdf",
  plot = p_lollipop_with_line,
  width = 5.6,
  height = 3
)

sv_dependence(shap_vis, sorted_features) + visualization_theme 
plots <- lapply(sorted_features, function(feat) {
  sv_dependence(shap_vis, feat) + visualization_theme
})
library(cowplot)
plot_grid(plotlist = plots, ncol = 3)
ggsave("2_SHAP_feature_dependence.pdf",width = 12,height =6)
####################



library(pcutils)
library(ggstatsplot)
library(PMCMRplus) 
library(ggplot2) 
library(cowplot) 

set.seed(123)

#包含基因名和分组这两列
data <- merged_data
data<-as.data.frame(data)
{
  # 请根据您的实际数据结构进行调整
  plist = list()
  for (i in 2:10) {
    plist[[i]] = group_box(tab = data[c("CYB5R2","MMD","CFB","SERPINA1")], group = "group", metadata = data, 
                           mode = i,p_value1 = "t.test") +#t.test，wilcox.test，anova，kruskal.test等方法进行比较
      ggtitle(paste0("mode", i)) +
      theme_classic() + 
      theme(legend.position = "none")
  }
  
  # 绘制箱式图
  plot_grid(plotlist = plist, ncol = 3)
}
#############
# 加载pROC包，用于绘制ROC曲线
library(pROC)    

# 设置工作路径，文件读写都将在这个路径下进行

# 读入CSV文件，第一列作为行名
rt=data
rt=rt[,c("group","CYB5R2","MMD","CFB","SERPINA1")]
# 选择目标变量y，通常是样本的分类标签
y=colnames(rt)[1]

# 定义绘图时使用的颜色，
#如果列数超过4，则使用彩虹颜色
bioCol=c("red","blue","green","yellow")
if(ncol(rt)>4){
  bioCol=rainbow(ncol(rt))}

#绘制
pdf("ROC(多基因).pdf",width=5,height=5)

# 绘制第一个基因的ROC曲线
roc1=roc(rt[,y], as.vector(rt[,2]))
aucText=c( paste0(colnames(rt)[2],", AUC=",sprintf("%0.3f",auc(roc1))) )
plot(roc1, col=bioCol[1])

# 绘制其他基因的ROC曲线，循环处理剩余基因
for(i in 3:ncol(rt)){
  roc1=roc(rt[,y], as.vector(rt[,i]))
  lines(roc1, col=bioCol[i-1])
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%0.3f",auc(roc1))) )
}

# 添加图例，显示各基因的AUC值
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
# 关闭PDF，完成绘图
dev.off()




# 单独展示第某一个模型
plot(plist[[2]])








#######
data<-read.csv("validation.csv")
library(pcutils)
library(ggstatsplot)
library(PMCMRplus) 
library(ggplot2) 
library(cowplot) 

set.seed(123)

#包含基因名和分组这两列

{
  # 请根据您的实际数据结构进行调整
  plist = list()
  for (i in 2:5) {
    plist[[i]] = group_box(tab = data[c("CYB5R2","MMD","CFB","SERPINA1")], group = "group", metadata = data, 
                           mode = i,p_value1 = "t.test") +#t.test，wilcox.test，anova，kruskal.test等方法进行比较
      ggtitle(paste0("mode", i)) +
      theme_classic() + 
      theme(legend.position = "none")
  }
  
  # 绘制箱式图
  plot_grid(plotlist = plist, ncol = 3)
}

plot(plist[[2]])


#############
# 加载pROC包，用于绘制ROC曲线
library(pROC)    

# 设置工作路径，文件读写都将在这个路径下进行

# 读入CSV文件，第一列作为行名
rt=data
rt=rt[,c("group","CYB5R2","MMD","CFB","SERPINA1")]
# 选择目标变量y，通常是样本的分类标签
y=colnames(rt)[1]

# 定义绘图时使用的颜色，
#如果列数超过4，则使用彩虹颜色
bioCol=c("red","blue","green","yellow")
if(ncol(rt)>4){
  bioCol=rainbow(ncol(rt))}

#绘制
pdf("ROC(多基因).pdf",width=5,height=5)

# 绘制第一个基因的ROC曲线
roc1=roc(rt[,y], as.vector(rt[,2]))
aucText=c( paste0(colnames(rt)[2],", AUC=",sprintf("%0.3f",auc(roc1))) )
plot(roc1, col=bioCol[1])

# 绘制其他基因的ROC曲线，循环处理剩余基因
for(i in 3:ncol(rt)){
  roc1=roc(rt[,y], as.vector(rt[,i]))
  lines(roc1, col=bioCol[i-1])
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%0.3f",auc(roc1))) )
}

# 添加图例，显示各基因的AUC值
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
# 关闭PDF，完成绘图
dev.off()




# 单独展示第某一个模型
plot(plist[[2]])

################
######
library(tidyverse)
library(GEOquery)
library(tidyverse)
library(GEOquery)
library(limma) 
library(affy)
library(stringr)
gset <- getGEO('GSE56081', destdir=".",
               AnnotGPL = F,     
               getGPL = F)  
exp<-exprs(gset[[1]])
#exp<-log2(exp+1)
cli<-pData(gset[[1]])	
GPL<-read.csv("36311041.csv",sep = ",",header = T)
GPL<-GPL[,c(1,2)]
#GPL<-fData(gset[[1]])
gpl<-GPL[,c(1,2)]        
gpl$GeneSymbols<-data.frame(sapply(gpl$GeneSymbols,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]
exp<-as.data.frame(exp)
exp$ProbeName<-rownames(exp)
exp_symbol<-merge(exp,gpl,by="ProbeName")
exp_symbol<-na.omit(exp_symbol)
table(duplicated(exp_symbol$Gene.Symbol))
exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$GeneSymbols)
write.table(exp_unique,"exp56081.csv",sep=",")
write.table(cli,"group56081.csv",sep=",")

#########
library(tidyverse)

# 1. 用基础包索引筛选列（替代 select()：保留 group 列和 4 个特征列）
# 先确认 group 列和特征列的列名是否在 data 中（避免列名拼写错误）
colnames(data)  # 运行此代码，检查是否有 "group", "CYB5R2", "MMD", "CFB", "SERPINA1"

# 筛选列：直接指定列名，用 [, c(列名)] 索引
data_subset <- data[, c("group", "CYB5R2", "MMD", "CFB", "SERPINA1")]
data_subset<-read.csv("rt.csv")
data_subset<-data_subset[,-2]
# 2. 转为长格式（用 pivot_longer()，此时数据已为数据框，可正常运行）
data_long <- data_subset %>%
  pivot_longer(
    cols = -group,  # 除了 group 列，其他列转为长格式
    names_to = "indexes",
    values_to = "value"
  )

# 验证结果
head(data_long)

# 2. 初始化图列表（避免空索引）

library(ggplot2)
library(ggpubr)  
library(cowplot)

target_features <- "HRTDGS1" 
plist <- list()  

for (i in 1:length(target_features)) {
  data_current <- data_long %>%
    filter(indexes == target_features[i])
  
  # 提前计算关键位置参数（确保均值文字在顶端，与P值上下呼应）
  max_y <- max(data_current$value, na.rm = TRUE)  # 数据最大值
  p_y <- max_y * 1.25  # P值位置（顶端稍下）
  mean_text_y <- max_y * 1.15  # 均值文字位置（P值正下方，顶端区域）
  
  plist[[i]] <- ggplot(data_current, aes(x = group, y = value, fill = group)) +
    # 1. 箱线图主体（底层，不遮挡上层元素）
    geom_boxplot(alpha = 0.7, width = 0.6, color = "black") +
    # 2. 散点（底层，抖动避免重叠）
    geom_jitter(width = 0.1, alpha = 0.6, size = 1.5, color = "black") +
    
    # ---------------------- 核心调整：均值数值放顶端 ----------------------
  stat_summary(
    aes(
      group = group, 
      label = round(..y.., 2),  # 均值保留2位小数
      y = mean_text_y  # 强制将均值文字放在顶端固定y位置（mean_text_y）
    ),
    fun = mean,          # 计算均值（仅用于获取label数值）
    geom = "text",       # 文字类型
    color = "darkblue",  # 文字颜色（醒目）
    size = 4,            # 文字放大，顶端更清晰
    vjust = 0.5,         # 垂直居中（因y已固定，无需偏移）
    show.legend = FALSE
  ) +
    
    # ---------------------- 可选：均值点仍保留在箱线图上方（不冲突） ----------------------
  stat_summary(
    aes(group = group),
    fun = mean,
    geom = "point",
    color = "darkred",
    size = 4,
    show.legend = FALSE
  ) +
    
    # 3. P值标注（放在均值文字正上方，顶端最外侧）
    stat_compare_means(
      method = "t.test",
      label = "p",
      label.y = p_y,  # P值在均值文字上方（p_y > mean_text_y）
      size = 4,       # P值字体放大，与均值文字匹配
      color = "darkred"
    ) +
    
    # 4. 图表标题和轴标签
    ggtitle(target_features[i]) +
    labs(x = "Group", y = "Expression Level") +
    
    # 5. 主题设置（确保y轴范围足够容纳顶端标注）
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 11, face = "bold"),
      # 可选：调整y轴边距，避免顶端标注被截断
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
    ) +
    
    # 6. 箱线图填充色
    scale_fill_brewer(palette = "Set2") +
    
    # ---------------------- 关键：扩大y轴上限，确保顶端标注不被截断 ----------------------
  ylim(min(data_current$value, na.rm = TRUE) * 0.9, p_y * 1.1)
}

# 组合多图
final_plot <- plot_grid(plotlist = plist, ncol = 1, labels = "AUTO")

# 查看最终图表
print(final_plot)
核心调整逻辑（确保均值在顶端且

# 7. 保存高分辨率图表（适合论文/汇报使用）
ggsave(
  filename = "Target_Genes_Boxplot_With_P.png",  # 保存文件名
  plot = final_plot,
  width = 8,  # 图宽（英寸）
  height = 6, # 图高（英寸）
  dpi = 300   # 分辨率（300dpi为印刷标准，避免模糊）
)

# 提示保存完成
cat("图表已保存为：Target_Genes_Boxplot_With_P.png\n")
















###ciber
library(CIBERSORT)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(pheatmap)
library(tibble)
library(tidyr)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(reshape2)
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
#表达矩阵文件
exp <- read.csv(file="expmerge.csv",row.names = 1,sep = ",")
exp<-2^exp
res <- cibersort(sig_matrix, exp, perm = 50, QN = F)

write.table(res, "immune.csv", 
            sep = ",", row.names = T, col.names = T, quote = F)








a <- read.csv("immune.csv", row.names = 1, header = T, as.is = F)
head(a)
a$X<-rownames(a)
fen<-read.csv("fen.csv")
a<-merge(fen,a,"X")
a
mydata1<-melt(a,
              id.vars=c("X","fen"),
              variable.name="immunecell",
              value.name="tpm")
ylabname <- paste("immunecell", "expression")
colnames(mydata1) <- c("Sample", "Groups", "immunecell","tpm")
mydata1<-mydata1[(1:748),]
#write.csv(mydata1,"mydata.csv")
ggplot(mydata1,aes(immunecell,tpm,fill = Groups)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "细胞类型", y = "细胞比例") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 14,face = "italic",colour = 'black'),
        axis.text.y = element_text(face = "italic",size = 14,colour = 'black'))+
  scale_fill_nejm()+
  stat_compare_means(aes(group = Groups),label = "p.format",size=3,method = "kruskal.test")

ggsave("immunecellbox.pdf", width = 18, height = 8,dpi = 300)








#####
CTSD
CYP1B1
METRNL
gene <- "CTSD"
exp <- read.csv(file="expmerge.csv",row.names = 1,sep = ",")
exp<-2^exp
exp<-t(exp)
exp=as.matrix(exp)
y <- as.numeric(exp[,gene])

### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation <- data.frame()
immu_data<-read.csv("immune.csv",row.names = 1)
### 2.批量把数据导出到容器
for(i in 1:length(colnames(immu_data))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(immu_data[,i]),y,method="spearman")
  ## 3.填充
  correlation[i,1] = colnames(immu_data)[i]
  correlation[i,2] = dd$estimate
  correlation[i,3] = dd$p.value
}
correlation<-na.omit(correlation)
### 修改名称
colnames(correlation) <- c("Cell","cor","pvalue")
View(correlation)
correlation<-correlation[(1:22),]
write.csv(correlation,"CSTD.csv")
#bangbangtangtu
p.col = c('gold','pink','orange','LimeGreen','darkgreen')
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                                                             ifelse(x>0.2,p.col[4], p.col[5])
  )))
  return(color)
} 
p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
                                                           ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)
} 
data<-read.csv("CSTD.csv")
#data<-data[-2,]
points.color = fcolor(x=data$pvalue,p.col=p.col)
data$points.color = points.color
points.cex = fcex(x=data$cor)
data$points.cex = points.cex
data=data[order(data$cor),] 
xlim = ceiling(max(abs(data$cor))*10)/10
pdf(file="CSTD棒棒糖图1.pdf", width=9, height=7)
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(data)+0.5),xlab="Correlation Coefficient of CSTD",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(data),col="white",lty=1,lwd=2) 
segments(x0=data$cor,y0=1:nrow(data),x1=0,y1=1:nrow(data),lwd=4)
#绘制图形的圆圈
points(x=data$cor,y = 1:nrow(data),col = data$points.color,pch=16,cex=data$points.cex)
#展示免疫细胞的名称
text(par('usr')[1],1:nrow(data),data$Cell,adj=1,xpd=T,cex=1.5) 
pvalue.text=ifelse(data$pvalue<0.001,'<0.001',sprintf("%.03f",data$pvalue))
redcutoff_cor=0
redcutoff_pvalue=0.05
text(par('usr')[2],1:nrow(data),pvalue.text,adj=0,xpd=T,col=ifelse(abs(data$cor)>redcutoff_cor & data$pvalue<redcutoff_pvalue,"red","black"),cex=1.5)
axis(1,tick=F)
par(mar=c(0,4,3,4))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")

#绘制圆圈颜色的图例
par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
dev.off()








gene <- "CYP1B1"
exp <- read.csv(file="expmerge.csv",row.names = 1,sep = ",")
exp<-2^exp
exp<-t(exp)
exp=as.matrix(exp)
y <- as.numeric(exp[,gene])

### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation <- data.frame()
immu_data<-read.csv("immune.csv",row.names = 1)
### 2.批量把数据导出到容器
for(i in 1:length(colnames(immu_data))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(immu_data[,i]),y,method="spearman")
  ## 3.填充
  correlation[i,1] = colnames(immu_data)[i]
  correlation[i,2] = dd$estimate
  correlation[i,3] = dd$p.value
}
correlation<-na.omit(correlation)
### 修改名称
colnames(correlation) <- c("Cell","cor","pvalue")
View(correlation)
correlation<-correlation[(1:22),]
write.csv(correlation,"CYP1B1.csv")
#bangbangtangtu
p.col = c('gold','pink','orange','LimeGreen','darkgreen')
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                                                             ifelse(x>0.2,p.col[4], p.col[5])
  )))
  return(color)
} 
p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
                                                           ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)
} 
data<-read.csv("CYP1B1.csv")
#data<-data[-2,]
points.color = fcolor(x=data$pvalue,p.col=p.col)
data$points.color = points.color
points.cex = fcex(x=data$cor)
data$points.cex = points.cex
data=data[order(data$cor),] 
xlim = ceiling(max(abs(data$cor))*10)/10
pdf(file="CYP1B1棒棒糖图1.pdf", width=9, height=7)
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(data)+0.5),xlab="Correlation Coefficient of CYP1B1",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(data),col="white",lty=1,lwd=2) 
segments(x0=data$cor,y0=1:nrow(data),x1=0,y1=1:nrow(data),lwd=4)
#绘制图形的圆圈
points(x=data$cor,y = 1:nrow(data),col = data$points.color,pch=16,cex=data$points.cex)
#展示免疫细胞的名称
text(par('usr')[1],1:nrow(data),data$Cell,adj=1,xpd=T,cex=1.5) 
pvalue.text=ifelse(data$pvalue<0.001,'<0.001',sprintf("%.03f",data$pvalue))
redcutoff_cor=0
redcutoff_pvalue=0.05
text(par('usr')[2],1:nrow(data),pvalue.text,adj=0,xpd=T,col=ifelse(abs(data$cor)>redcutoff_cor & data$pvalue<redcutoff_pvalue,"red","black"),cex=1.5)
axis(1,tick=F)
par(mar=c(0,4,3,4))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")

#绘制圆圈颜色的图例
par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
dev.off()








gene <- "METRNL"
exp <- read.csv(file="expmerge.csv",row.names = 1,sep = ",")
exp<-2^exp
exp<-t(exp)
exp=as.matrix(exp)
y <- as.numeric(exp[,gene])

### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation <- data.frame()
immu_data<-read.csv("immune.csv",row.names = 1)
### 2.批量把数据导出到容器
for(i in 1:length(colnames(immu_data))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(immu_data[,i]),y,method="spearman")
  ## 3.填充
  correlation[i,1] = colnames(immu_data)[i]
  correlation[i,2] = dd$estimate
  correlation[i,3] = dd$p.value
}
correlation<-na.omit(correlation)
### 修改名称
colnames(correlation) <- c("Cell","cor","pvalue")
View(correlation)
correlation<-correlation[(1:22),]
write.csv(correlation,"METRNL.csv")
#bangbangtangtu
p.col = c('gold','pink','orange','LimeGreen','darkgreen')
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                                                             ifelse(x>0.2,p.col[4], p.col[5])
  )))
  return(color)
} 
p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
                                                           ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)
} 
data<-read.csv("METRNL.csv")
#data<-data[-2,]
points.color = fcolor(x=data$pvalue,p.col=p.col)
data$points.color = points.color
points.cex = fcex(x=data$cor)
data$points.cex = points.cex
data=data[order(data$cor),] 
xlim = ceiling(max(abs(data$cor))*10)/10
pdf(file="METRNL棒棒糖图1.pdf", width=9, height=7)
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(data)+0.5),xlab="Correlation Coefficient of METRNL",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(data),col="white",lty=1,lwd=2) 
segments(x0=data$cor,y0=1:nrow(data),x1=0,y1=1:nrow(data),lwd=4)
#绘制图形的圆圈
points(x=data$cor,y = 1:nrow(data),col = data$points.color,pch=16,cex=data$points.cex)
#展示免疫细胞的名称
text(par('usr')[1],1:nrow(data),data$Cell,adj=1,xpd=T,cex=1.5) 
pvalue.text=ifelse(data$pvalue<0.001,'<0.001',sprintf("%.03f",data$pvalue))
redcutoff_cor=0
redcutoff_pvalue=0.05
text(par('usr')[2],1:nrow(data),pvalue.text,adj=0,xpd=T,col=ifelse(abs(data$cor)>redcutoff_cor & data$pvalue<redcutoff_pvalue,"red","black"),cex=1.5)
axis(1,tick=F)
par(mar=c(0,4,3,4))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")

#绘制圆圈颜色的图例
par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
dev.off()







##########
rm(list=ls())
#以后seurat包可根据安装路径指定版本
#.libPaths(c("D:/R语言/Rlibrary/win-library/4.2","~/seurat_v5","C:/Program Files/R/R-4.3.1/library"))
library(Seurat)
library(stringr)
library(dplyr)
library(future)
library(future.apply)
library(dplyr)
library(msigdbr)
library(clusterProfiler)
library(devtools)
library(harmony)
library(clustree)
library(ggplot2)
library(reshape2)
library(data.table)

options(future.globals.maxSize = 60000 * 1024^2)
getwd()

#define the color
library(ggsci)
cors <- pal_igv()(12) #定义颜色

getwd()
#setwd("E:/科研-合作/单细胞/单细胞联合bulk教学/codes")

#### 读取单细胞数据 ####
## 循环读取GSE184198的数据
filename <- paste('GSENP/',list.files('GSENP/'),sep = '')
sceList <- lapply(filename, function(x){
  obj <- CreateSeuratObject(counts = Read10X(x),
                            min.cells = 3, min.features = 50,
                            project = str_split(x,'/')[[1]][2])
})
names(sceList) <- list.files('GSENP/')
scRNA <- merge(sceList[[1]],sceList[-1],add.cell.ids = names(sceList),project='LDH')
dim(scRNA) # 23319个gene， 22217个cell
scRNA@meta.data[1:5,] #简单看一下


#### 质量控制 ####
##计算质控指标
#计算细胞中线粒体核糖体基因比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
library(ggplot2)
col.num <- length(levels(scRNA@active.ident))
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(scRNA))
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes)
##设置质控标准，比较随意
print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=500", "maxGene=4000", "pctMT=20"))
minGene=200
maxGene=5000
pctMT=10
pctHB=0
FeatureScatter(scRNA, "nCount_RNA", "percent.mt", group.by = "orig.ident")
FeatureScatter(scRNA, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident")

#查看质控指标
#设置绘图元素
##数据质控
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA< maxGene & percent.mt < pctMT& percent.HB < 0.00001 )
col.num <- length(levels(scRNA@active.ident))
VlnPlot(scRNA,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        cols =rainbow(col.num), 
        pt.size = 0, 
        ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

# 标准化
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA@meta.data$group <- ifelse(substr(rownames(scRNA@meta.data),1,2)=='LD','LDH','NC')

#### 降维聚类 ####
library(Seurat)
library(tidyverse)
library(patchwork)

scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
scRNA<- RunPCA(scRNA, features = VariableFeatures(scRNA))
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")
#######
pbmc=scRNA
ElbowPlot(pbmc, ndims = 20)


#设置PC
pcs = 1:15

#harmony
#pbmc <- RunHarmony(pbmc, group.by.vars="orig.ident", assay.use="RNA", max.iter.harmony = 20)

table(pbmc@meta.data$orig.ident)

#选取合适的分辨率
#从0.1-2的resolution结果均运行一遍
seq = seq(0.1,2,by=0.1)
pbmc <- FindNeighbors(pbmc,  dims = pcs) 
for (res in seq){
  pbmc = FindClusters(pbmc, resolution = res)
}

#画图
p1 = clustree(pbmc,prefix = "RNA_snn_res.")+coord_flip()
p = p1+plot_layout(widths = c(3,1))
ggsave("2.png", p1, width = 30, height = 14)

#降维聚类
#pbmc <- FindNeighbors(pbmc, reduction = "harmony",  dims = pcs) %>% FindClusters(resolution = 0.6)
#pbmc <- RunUMAP(pbmc, reduction = "harmony",  dims = pcs) %>% RunTSNE(dims = pcs, reduction = "harmony")
#降维聚类
pbmc=scRNA
pbmc <- FindNeighbors(pbmc, reduction = "pca",  dims = pcs) %>% FindClusters(resolution = 0.3)
pbmc <- RunUMAP(pbmc, reduction = "pca",  dims = pcs) #%>% RunTSNE(dims = pcs, reduction = "pca")

library(scRNAtoolVis)
library(ggpubr)
library(ggsci)
#pbmc <- readRDS("umap.rds")

P=DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters") 
pdf(file=paste0('clustersumap.pdf'),width = 8,height = 6)
clusterCornerAxes(object = pbmc, reduction = 'umap',
                  clusterCol = "seurat_clusters",
                  addCircle = TRUE, cicAlpha = 0.1, nbin = 150,  
                  cellLabel = TRUE, cellLabelSize = 5,cicDelta = 0.01) +  
  scale_color_igv() + scale_fill_igv()
dev.off()
table(pbmc$seurat_clusters)

#保存
saveRDS(yjsl,"umap.rds")
scRNA<-pbmc
Idents(scRNA)=scRNA$seurat_clusters
library(SingleR)








T cells: CD3D, CD3E, CD8A, PTPRC, CD4, CD2
B cells: CD19, CD79A, MS4A1 (CD20)
Plasma cells: IGHG1, MZB1, SDC1, CD79A
Monocytes / Macrophages: CD68, CD163, CD14
NK Cells: FGFBP2, FCGR3A, CX3CR1
Fibroblasts: FGF7, MME, GSN, LUM, DCN
Endothelial cells: PECAM1, VWF
Epithelial / Tumor: EPCAM, KRT19, PROM1, ALDH1A1, CD24
Mast cells: CPA3, CST3, KIT, TPSAB1, TPSB2, MS4A2

Idents(sce.all)=sce.$seurat_clusters


gene_marker <- list("MGP","GPX3","FN1",
                    "CD163","CD206",
                    "COL1A1","COL3A1","FAP",
                    "CD31","CD34")
features = gene_marker
length(features)
#all(features %in% rownames(sce.all))
#绘图

pdf(file=paste0('scRNA_cluster_vlnplot.pdf'),width = 6,height = 6)
VlnPlot(scRNA, features = features, #6*6
        stack = TRUE, 
        sort = TRUE, 
        cols = cors,
        split.by =  "seurat_clusters" , #每种cluster 一个颜色
        flip = TRUE,
        split.plot = TRUE) +
  theme(legend.position = "none")
dev.off()








Idents(yjsl)=yjsl@meta.data$seurat_clusters
genes <- list("chondrocyte"= c("MGP","GPX3","FN1"),
              "M1 Macrophage"="CD68",
              "M2 Macrophage"=c("CD163","CD206"),
              "Fibroblasts" = c("COL1A1","COL3A1","FAP"),
              "Endothelial cell"=c("CD31","CD34")
)
yjsl=pbmc
p1 <- DotPlot(yjsl, features = genes ,
              assay='RNA' ) + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
#另一种颜色气泡
#pbmc=scRNA
p1<-DotPlot(yjsl, features = genes,cols = "RdYlBu") +
  RotatedAxis()

DotPlot(yjsl, 
        features = genes,
        group.by = "seurat_clusters") + coord_flip()

#不要求非常准确
celltype <- c("Chondrocyte",              #cluster0
              "Chondrocyte",      #cluster1
              "Fibroblast",      #以下按顺序操作
              "Fibroblast",
              "Chondrocyte",
              "Chondrocyte",
              "Fibroblast",
              "Chondrocyte",
              "Chondrocyte",
              "Chondrocyte",
              "Fibroblast",
              "Chondrocyte")

scRNA$celltype <- recode(scRNA@meta.data$seurat_clusters,
                          "0" = "Chondrocyte",
                          "1" = "Chondrocyte",
                          "2" = "Fibroblast",
                          "3" = "Fibroblast",
                          "4" = "Chondrocyte",
                          "5" = "Chondrocyte",
                          "6" = "Fibroblast",
                          "7" = "Chondrocyte",
                          "8" = "Chondrocyte",
                          "9" = "Chondrocyte",
                          "10" = "Fibroblast",
                          "11" = "Chondrocyte")


## 重命名Idents

# 命名celltype
yjsl@meta.data$celltype <- Idents(yjsl)

## 保存文件
saveRDS(yjsl,'pbmc.rds')


scRNA=yjsl
### 聚类可视化
## umap图
DimPlot(scRNA,reduction = "umap",label = T,group.by="group",cols = cors)
DimPlot(scRNA,reduction = "umap",label = T,group.by="celltype",split.by="orig.ident")


#scRNA=pbmc
####注释后可视化####
#美化教程：https://mp.weixin.qq.com/s/1C64jsl08oTVoxbNQ49ioA
p=clusterCornerAxes(object = scRNA, reduction = "umap", clusterCol = "SingleR", cicDelta = 0.5) 
p=clusterCornerAxes(object = scRNA, reduction = "umap", clusterCol = "seurat_clusters", cicDelta = 0.5) 
##绘制nature级别umap图   #7*5
#library(igvWidgets)
library(scRNAtoolVis)
pdf(file=paste0('scRNA_celltype1.pdf'),width = 7,height = 5)
clusterCornerAxes(object = yjsl, 
                  reduction = 'umap',
                  clusterCol = "celltype",
                  addCircle = TRUE, cicAlpha = 0.1, nbin = 200,
                  cellLabel = TRUE, cellLabelSize = 3.5, cicDelta = 0.02) +
  scale_color_brewer(palette = "Set3") +  # 点/标签颜色
  scale_fill_brewer(palette = "Set3") 
dev.off()

##绘制各细胞marker小提琴图
#确定选用展示的markers（根据经验选择）


#绘图
pdf(file=paste0('scRNA_celltype_vlnplot.pdf'),width = 6,height =8)
VlnPlot(scRNA, features = features, #6*6
        stack = TRUE, 
        sort = TRUE, 
        cols = cors,
        split.by =  "seurat_clusters" , #每种cluster 一个颜色
        flip = TRUE,
        split.plot = TRUE) +
  theme(legend.position = "none")
dev.off()



#画图、
#pbmc<-readRDS("pbmc_SINGLERTUMOR.umap.rds")
library(scRNAtoolVis)
library(ggpubr)
library(ggsci)
pdf(file=paste0('celltypeumap.pdf'),width = 8,height = 6)
clusterCornerAxes(object = pbmc, reduction = 'umap',
                  clusterCol = "SingleR", #分类依据
                  #noSplit = F, groupFacet = 'group', #分组时使用
                  addCircle = TRUE, cicAlpha = 0.2, nbin = 200,  #画圈
                  cellLabel = T, cellLabelSize = 3,cicDelta = 0.1) +  #标签
  scale_color_igv() + scale_fill_igv()
dev.off()
############
scRNA<-readRDS("pbmc.rds")


df <- as.data.frame(GetAssayData(
  object = scRNA, 
  assay = "RNA",  # 通常是RNA assay
  slot = "data"   # slot可选：data(归一化后)、counts(原始计数)、scale.data(标准化后)
))


key_genes <- c("CSTD","CYP1B1","METRNL")
df_filtered <- df[rownames(df) %in% key_genes, ]  # 只保留关键基因的行
library(tidyverse)
# 或者：只保留特定细胞亚群（需要先从col_info筛选细胞名）
col_info <- scRNA@meta.data %>% 
  # 提取需要的列：细胞名、亚群、分组（列名根据你的Seurat对象调整）
  rownames_to_column("sample_full") %>%  # 把行名（细胞ID）转为sample_full列
  select(sample_full, celltype, group)    # 只保留这三列，确保列名和后续代码匹配

# 再运行你的重塑代码
dt <- df_filtered %>%
  rownames_to_column("genes") %>%
  pivot_longer(-genes, names_to = "sample_full", values_to = "expressions") %>%
  left_join(col_info, by = "sample_full") %>%
  transmute(celltype, genes, expressions, group)



write.csv(dt,"singlecellexpression.csv")

head(dt)
tail(dt)

#dt <- fread("5-pseudobulks/dt_cell.txt")
gene=c("CYP1B1","METRNL")
gene <- c("CYP1B1","METRNL")
dt2 <- dt %>% 
  filter(genes %in% gene) 
#dt2 <- dt[genes %in% gene, ]
dt2$genes <- factor(dt2$genes, levels = unique(dt2$genes))
dt <- dt2
library(ggplot2)
library(gghalves)
p <- ggplot() +
  geom_half_violin(
    data = dt[dt$group == "NC", ],
    aes(x = genes, y = expressions, fill = group),
    color = "black",
    linewidth = 0.4,
    draw_quantiles = c(0.5),
    scale = "width"
  ) +
  facet_grid(rows = vars(celltype), scales = "free_y")
p

mytheme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 60, hjust = 1),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12, angle = 0),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

fill_col <- c(NC = "#80B1D3", LDH = "#FDB462")
point_col <- c(NC = "#BF5B17", LDH = "#F0027F")

p <- ggplot() +
  geom_half_violin(
    data = dt %>% filter(group == "NC"),
    aes(x = genes, y = expressions, fill = group),
    color = "black",
    linewidth = 0.4,
    draw_quantiles = c(0.5),
    scale = "width",
    side = "l"
  ) +
  facet_grid(rows = vars(celltype), scales = "free_y")

p1 <- p +
  geom_half_violin(
    data = dt %>% filter(group == "LDH"),
    aes(x = genes, y = expressions, fill = group),
    color = "black",
    linewidth = 0.4,
    draw_quantiles = c(0.5),
    scale = "width",
    side = "r"
  ) +
  facet_grid(rows = vars(celltype), scales = "free_y")

p2 <- p1 +
  mytheme +
  scale_fill_manual(values = fill_col, breaks = c("NC", "LDH")) +
  scale_y_continuous(breaks = seq(0, 11, by = 5)) +
  labs(y = "Log Normalized Expression")

valid_panels <- dt %>%
  group_by(celltype, genes) %>%
  summarise(n_group = n_distinct(group), .groups = "drop") %>%
  filter(n_group == 2)

dt_valid <- dt %>%
  inner_join(valid_panels %>% select(celltype, genes), by = c("celltype", "genes"))
library(ggpubr)



library(rstatix)
sample_count <- dt_valid %>%
  group_by(celltype, genes, group) %>%
  summarise(n = n(), .groups = "drop_last") %>%  # 按3级分组统计样本数
  filter(n >= 2) %>%  # 只保留样本量≥2的group
  summarise(valid_group = n() >= 2, .groups = "drop")  # 要求该基因-亚群下至少有2个有效group

# 步骤2：筛选出所有有效分组的数据集
dt_valid_filtered <- dt_valid %>%
  inner_join(sample_count %>% filter(valid_group), by = c("celltype", "genes"))

# 步骤3：再做t检验 + 添加显著性
pvals <- dt_valid_filtered %>%
  group_by(celltype, genes) %>%
  t_test(expressions ~ group) %>%
  add_significance(p.col = "p") %>%
  ungroup()

pvals %>% select(celltype, genes, p, p.signif)

ypos <- dt_valid %>%
  group_by(celltype, genes) %>%
  summarise(y.position = max(expressions, na.rm = TRUE) * 1.05, .groups = "drop")

p_anno <- pvals %>%
  left_join(ypos, by = c("celltype", "genes")) %>%
  mutate(
    xmin = genes,
    xmax = genes,
    label = p.signif
  )

p_final <- p2 +
  ggpubr::stat_pvalue_manual(
    p_anno,
    label = "label",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    tip.length = 0,
    bracket.size = 0
  )

pdf("pea1.pdf", width = 6, height = 5)
p_final
dev.off()
install.packages("plot1cell")

devtools::install_github("TheHumphreysLab/plot1cell")

BiocManager::install("EnsDb.Hsapiens.v86")
library(plot1cell)
library(RColorBrewer)

sce <- scRNA
Idents(sce) <- sce$celltype
colnames(sce@meta.data)
sce <- subset(sce, subset = group != "HIS")
table(sce$group)

library(SingleCellExperiment)
library(circlize)
library(dplyr)


prepare_circlize_data <- function(sce, scale = 0.8) {
  # 1. 提取Seurat对象的元数据（细胞注释：亚群、分组等）
  # 注意：根据你的Seurat对象元数据列名修改 cluster/group 对应的列名
  meta_data <- sce@meta.data %>%
    rownames_to_column("cell") %>%  # 把行名（细胞ID）转为列
    select(cell, cluster = seurat_clusters, group = group)  # 替换为你的实际列名
  
  # 2. 提取归一化后的表达矩阵（取高变基因，减少数据量）
  # 若没有高变基因，可改为 VariableFeatures(sce) 或直接指定目标基因
  if (length(VariableFeatures(sce)) > 0) {
    expr_mat <- GetAssayData(
      object = sce,
      assay = "RNA",  # 根据你的assay修改
      slot = "data"   # data=归一化后，counts=原始计数
    )[VariableFeatures(sce)[1:500], ]  # 取前500个高变基因，避免内存溢出
  } else {
    # 无高变基因时取前200个基因
    expr_mat <- GetAssayData(sce, assay = "RNA", slot = "data")[1:200, ]
  }
  
  # 3. 表达量缩放（按scale参数调整）
  expr_mat_scaled <- t(scale(t(as.matrix(expr_mat))) * scale)
  
  # 4. 整理为circlize可用的长格式数据框
  circ_data <- as.data.frame(expr_mat_scaled) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta_data, by = "cell")  # 合并细胞注释信息
  
  return(circ_data)
}

# 现在调用函数（sce是Seurat对象）
circ_data <- prepare_circlize_data(sce, scale = 0.8)


#circ_data <- prepare_circlize_data(sce, scale = 0.8)
set.seed(1234)
cluster_colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(sce@meta.data$celltype)))
group_colors <- rand_color(length(names(table(sce$group))))
rep_colors <- rand_color(length(names(table(sce$orig.ident))))
library(circlize)
library(tidyverse)
library(RColorBrewer)

  
pdf("1circlize.pdf", height = 10, width = 10)
plot_circlize(
  circ_data,
  do.label = TRUE,
  pt.size = 0.15,
  col.use = cluster_colors,
  bg.color = FALSE,
  kde2d.n = 200,
  repel = FALSE,
  label.cex = 1
)

add_track(circ_data, group = "group", colors = group_colors, track_num = 2)
add_track(circ_data, group = "orig.ident", colors = rep_colors, track_num = 3)
dev.off()

cell_counts <- as.data.frame(table(sce@meta.data$celltype, sce$group))

p=ggplot(
  data = cell_counts,
  aes(x = forcats::fct_rev(Var2), y = Freq, fill = Var1)
) +
  geom_bar(position = "fill", stat = "identity") +
  coord_flip() +
  labs(x = NULL, y = "Cell Lineage Frequency") +
  scale_x_discrete(labels = c("HC", "Pso")) +
  scale_fill_manual(
    name = "Cell lineage",
    values = cluster_colors
  ) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(size = 16, color = "black"))
p
setwd("../")
dir.create("subT")
setwd("subT")

my_sub <- "T"
seu.obj <- sce.all.int
sub.cells <- subset(seu.obj, idents = my_sub)
f <- "obj.Rdata"

if (!file.exists(f)) {
  sub.cells <- sub.cells %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData(features = rownames(.)) %>%
    RunPCA(features = VariableFeatures(.)) %>%
    FindNeighbors(dims = 1:15) %>%
    FindClusters(resolution = 0.5) %>%
    RunUMAP(dims = 1:15)
  save(sub.cells, file = f)
}

load(f)
DimPlot(sub.cells, reduction = "umap", label = TRUE) + NoLegend()

library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

CD4_markers_list <- list(
  Tc = c("CD3D", "CD3E"),
  CD4 = c("CD4"),
  Treg = c("TNFRSF4", "BATF", "TNFRSF18", "FOXP3", "IL2RA", "IKZF2"),
  naive = c("CCR7", "SELL", "CD5"),
  Tfh = c("CXCR5", "BCL6", "ICA1", "TOX", "TOX2", "IL6ST"),
  ILC = c("TNFRSF25", "KRT81", "LST1", "AREG", "LTB", "CD69")
)

CD8_markers_list1 <- list(
  CD8 = c("CD8A", "CD8B"),
  TN_TCM = c("CCR7", "SELL", "TCF7", "LEF1"),
  TEM = c("GZMK"),
  TEFF = c("TBX21", "FCGR3A", "FGFBP2"),
  TRM = c("XCL1", "XCL2", "ITGAE", "CD69"),
  IEL_T = c("TMIGD2"),
  yT1c = c("GNLY", "PTGDS", "GZMB", "TRDC"),
  yT2c = c("TMN1", "HMGB2", "TYMS"),
  MAIT_T = c("SLC4A10")
)

CD8_markers_list2 <- list(
  CD8T = c("CD8A", "CD8B"),
  MAIT = c("ZBTB16", "NCR3", "RORA"),
  ExhaustedCD8T = c("LAG3", "TIGIT", "PDCD1", "HAVCR2", "CTLA4"),
  EffMemoryCD8 = c("EOMES", "ITM2C"),
  Resting_NK = c("XCL1", "XCL2", "KLRC1"),
  Cytotoxic_NK = c("CX3CR1", "FGFBP2", "FCGR3A", "KLRD1"),
  Pre_exhausted = c("IFNG", "PRF1", "GNLY", "GZMA", "NKG7", "GZMK")
)

cd4_and_cd8T_markers_list <- list(
  naive = c("CCR7", "SELL", "TCF7", "IL7R", "CD27", "CD28", "LEF1", "S1PR1"),
  CD8Trm = c("XCL1", "XCL2", "MYADM"),
  NKTc = c("GNLY", "GZMA"),
  Tfh = c("CXCR5", "BCL6", "ICA1", "TOX", "TOX2", "IL6ST"),
  th17 = c("IL17A", "KLRB1", "CCL20", "ANKRD28", "IL23R", "RORC", "FURIN", "CCR6", "CAPG", "IL22"),
  CD8Tem = c("CXCR4", "GZMH", "CD44", "GZMK"),
  Treg = c("FOXP3", "IL2RA", "TNFRSF18", "IKZF2"),
  naive2 = c("CCR7", "SELL", "TCF7", "IL7R", "CD27", "CD28"),
  CD8Trm2 = c("XCL1", "XCL2", "MYADM"),
  MAIT = c("KLRB1", "ZBTB16", "NCR3", "RORA"),
  yT1c = c("GNLY", "PTGDS", "GZMB", "TRDC"),
  yT2c = c("TMN1", "HMGB2", "TYMS"),
  yt = c("TRGV9", "TRDV2")
)

markers_list <- c(
  "CD4_markers_list",
  "CD8_markers_list1",
  "CD8_markers_list2",
  "cd4_and_cd8T_markers_list"
)

Idents(sce.all.int) <- "celltype"

.clean_feature_list <- function(feat_list, gene_universe) {
  feat_list <- lapply(feat_list, stringr::str_to_upper)
  allv <- unlist(feat_list, use.names = FALSE)
  dup <- names(table(allv))[table(allv) > 1]
  feat_list <- lapply(feat_list, function(v) setdiff(v, dup))
  feat_list <- lapply(feat_list, function(v) intersect(v, gene_universe))
  feat_list[lengths(feat_list) > 0]
}

invisible(lapply(markers_list, function(x) {
  genes_to_check <- get(x, inherits = TRUE)
  genes_to_check <- .clean_feature_list(genes_to_check, rownames(sub.cells))
  
  if (length(genes_to_check) == 0L) {
    message("skip ", x, ": no valid genes")
    return(NULL)
  }
  
  p <- DotPlot(sub.cells, features = genes_to_check, assay = "RNA") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0("Marker check: ", x))
  
  w <- length(unique(unlist(genes_to_check))) / 5 + 6
  ggsave(filename = paste0("check_for_", x, ".pdf"), plot = p, width = w, height = 6)
}))

library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(DESeq2)
library(ggrepel)
library(ggsci)
library(RColorBrewer)
library(qs)

sce <- sce.all.int
sce <- subset(sce, subset = group != "HIS")
colnames(sce.all@meta.data)

dir.create("5-pseudobulks")
setwd("5-pseudobulks")

sce.all <- sce
av <- AggregateExpression(
  sce.all,
  group.by = c("orig.ident", "celltype"),
  assays = "RNA",
  return.seurat = FALSE
)
av <- as.data.frame(av[[1]])
head(av)[1:3, 1:3]

cg <- names(tail(sort(apply(log(av + 1), 1, sd)), 1000))
df <- cor(as.matrix(log(av[cg, ] + 1)))
colnames(df)

ac <- as.data.frame(str_split(colnames(df), "_", simplify = TRUE))
rownames(ac) <- colnames(df)
ac$V1

ac$group <- ifelse(
  grepl("GSM6840143|144|145|146|147|148|149|150|151|152", ac$V1),
  "Control",
  "Psoriasis"
)

table(ac$group)
head(ac)

pheatmap::pheatmap(
  df,
  show_colnames = FALSE,
  show_rownames = FALSE,
  annotation_col = ac,
  filename = "cor_celltype-vs-orig.ident.pdf"
)

save(av, file = "av_for_pseudobulks.Rdata")

av <- AggregateExpression(
  sce.all,
  group.by = c("orig.ident", "celltype"),
  assays = "RNA"
)
av <- as.data.frame(av[[1]])
df <- log(av + 1)
head(ac)
celltp <- unique(ac$V2)
celltp

library(FactoMineR)
library(factoextra)
library(ggstatsplot)
library(pheatmap)

pca_list <- lapply(celltp, function(x) {
  exp <- df[, rownames(ac[ac$V2 == x, ])]
  cg <- names(tail(sort(apply(exp, 1, sd)), 1000))
  exp <- exp[cg, ]
  dat.pca <- PCA(as.data.frame(t(exp)), graph = FALSE)
  group_list <- ac[ac$V2 == x, "group"]
  this_title <- paste0(x, "_PCA")
  
  p2 <- fviz_pca_ind(
    dat.pca,
    geom.ind = "point",
    col.ind = group_list,
    palette = "Dark2",
    addEllipses = TRUE,
    legend.title = "Groups"
  ) +
    ggtitle(this_title) +
    theme_ggstatsplot() +
    theme(plot.title = element_text(size = 12, hjust = 0.5))
  
  p2
})

wrap_plots(pca_list, byrow = TRUE, nrow = 2)
ggsave("all_pca.pdf", width = 12, height = 6)
ggsave("all_pca.png", dpi = 300, width = 12, height = 6)

av <- AggregateExpression(
  sce.all,
  group.by = c("orig.ident", "group"),
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)
av <- as.data.frame(av[[1]])
head(av)[1:3, 1:3]

library(tibble)
library(DESeq2)

counts_res <- av
head(counts_res)

colData <- data.frame(samples = colnames(counts_res))
colData <- colData %>%
  mutate(
    condition = ifelse(
      grepl("Psoriasis", samples),
      "Psoriasis",
      "Control"
    )
  ) %>%
  column_to_rownames(var = "samples")

dds <- DESeqDataSetFromMatrix(
  countData = counts_res,
  colData = colData,
  design = ~condition
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, name = "condition_Psoriasis_vs_Control")
res

gene <- c("OASL","SCO2","IFI6","NMI",
          "SAMD9","MX1","UBE2L6","OAS3")

res[gene, ]
res["SAMD9", ]
res["UBE2L6", ]

resOrdered <- res[order(res$padj), ]
head(resOrdered)

DEG <- as.data.frame(resOrdered)
DEG_deseq2 <- na.omit(DEG)

DEG_deseq2 <- DEG_deseq2 %>%
  mutate(
    Type = if_else(
      padj > 0.05, "stable",
      if_else(
        abs(log2FoldChange) < 1, "stable",
        if_else(log2FoldChange >= 1, "up", "down")
      )
    )
  ) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  rownames_to_column("Symbol")

ggplot(DEG_deseq2, aes(log2FoldChange, -log10(padj))) +
  geom_point(
    size = 3.5,
    alpha = 0.8,
    aes(color = Type),
    show.legend = TRUE
  ) +
  scale_color_manual(values = c("#00468B", "gray", "#E64B35")) +
  ylim(0, 15) +
  xlim(-10, 10) +
  labs(x = "Log2(fold change)", y = "-log10(padj)") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2,
    color = "black",
    lwd = 0.8
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = 2,
    color = "black",
    lwd = 0.8
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("volcano.pdf", width = 9, height = 7)

library(tidyverse)
library(ggstatsplot)

gene <- c("OASL", "RTP4", "SCO2", "IFI6", "LGALS3BP", "NMI", "IRF7",
          "SAMD9", "IFIT3", "MX1", "UBE2L6", "ISG15", "TAP1", "OAS3")

bdata <- as.data.frame(t(counts_res[gene, ]))
td <- bdata
class(td)

td$group <- "ss"
td[rownames(td) %like% "Psoriasis", ]$group <- "Psoriasis"
td[rownames(td) %like% "Control", ]$group <- "Control"
table(td$group)

datasets <- list(td)
dataset_names <- c("GSE220116")

for (i in seq_along(datasets)) {
  df <- datasets[[i]]
  name <- dataset_names[i]
  
  df_long <- df %>%
    pivot_longer(
      cols = -group,
      names_to = "gene",
      values_to = "Value"
    )
  
  genes <- unique(df_long$gene)
  
  for (g in genes) {
    subdata <- df_long %>%
      filter(gene == g & group %in% c("Psoriasis", "Control"))
    
    p <- ggbetweenstats(
      data = subdata,
      x = UC,
      y = Value,
      title = paste0("Expression of ", g, " in ", name),
      messages = FALSE
    )
    
    ggsave(
      filename = paste0(name, "_", g, ".pdf"),
      plot = p,
      width = 6,
      height = 5
    )
  }
}


gene <- c(
  "OASL","SCO2","IFI6","NMI",
  "SAMD9","MX1","UBE2L6","OAS3"
)

dt <- readRDS("dt.rds")
dt[1:8, ]

dt$genes <- factor(dt$genes, levels = unique(dt$genes))

library(gghalves)
library(ggplot2)
library(cols4all)
library(dplyr)
library(rstatix)
library(ggpubr)

## 分组：HC vs PSO（原来是 HC vs SLE）
dt <- dt %>%
  mutate(
    Group = factor(Group, levels = c("HC","PSO"))
  )

mytheme <- theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.text.y     = element_text(size = 12),
    axis.text.x     = element_text(size = 12, angle = 60, hjust = 1),
    axis.title.y    = element_text(size = 12),
    axis.title.x    = element_blank(),
    strip.background= element_blank(),
    strip.text.y    = element_text(size = 12, angle = 0),
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 12)
  )

## 颜色命名也改成 PSO
fill_col  <- c(HC = "#80B1D3", PSO = "#FDB462")
point_col <- c(HC = "#BF5B17", PSO = "#F0027F")

## 左：HC
p <- ggplot() +
  geom_half_violin(
    data   = dt %>% filter(Group == "HC"),
    aes(x  = genes, y = expressions, fill = Group),
    color  = "black",
    linewidth = 0.4,
    draw_quantiles = c(0.5),
    scale  = "width",
    side   = "l"
  ) +
  facet_grid(rows = vars(cluster), scales = "free_y")

## 右：PSO
p1 <- p +
  geom_half_violin(
    data   = dt %>% filter(Group == "PSO"),
    aes(x  = genes, y = expressions, fill = Group),
    color  = "black",
    linewidth = 0.4,
    draw_quantiles = c(0.5),
    scale  = "width",
    side   = "r"
  ) +
  facet_grid(rows = vars(cluster), scales = "free_y")

p2 <- p1 +
  mytheme +
  scale_fill_manual(values = fill_col, breaks = c("HC","PSO")) +
  scale_y_continuous(breaks = seq(0, 11, by = 5)) +
  labs(y = "Log Normalized Expression")

## 散点：HC
p3 <- p2 +
  geom_half_point(
    data = dt %>% filter(Group == "HC"),
    aes(x = genes, y = expressions, color = Group),
    range_scale = 0.35, size = 0.1, show.legend = FALSE, side = "l"
  )

## 散点：PSO
p4 <- p3 +
  geom_half_point(
    data = dt %>% filter(Group == "PSO"),
    aes(x = genes, y = expressions, color = Group),
    range_scale = 0.35, size = 0.1, show.legend = FALSE, side = "r"
  ) +
  scale_color_manual(values = scales::alpha(point_col, 0.3),
                     breaks = c("HC","PSO"))

## 显著性检验部分不用改，只是对 Group 这列做检验
valid_panels <- dt %>%
  group_by(cluster, genes) %>%
  summarise(n_group = n_distinct(Group), .groups = "drop") %>%
  filter(n_group == 2)

dt_valid <- dt %>%
  inner_join(valid_panels %>% select(cluster, genes),
             by = c("cluster","genes"))

pvals <- dt_valid %>%
  group_by(cluster, genes) %>%
  t_test(expressions ~ Group) %>%    # 或 wilcox_test
  adjust_pvalue(method = "BH") %>%
  add_significance(p.col = "p") %>%
  add_significance(p.col = "p.adj") %>%
  ungroup()

ypos <- dt_valid %>%
  group_by(cluster, genes) %>%
  summarise(y.position = max(expressions, na.rm = TRUE) * 1.05,
            .groups = "drop")

p_anno <- pvals %>%
  left_join(ypos, by = c("cluster","genes")) %>%
  mutate(
    xmin  = genes,
    xmax  = genes,
    label = p.adj.signif
  )

p_final <- p4 +
  ggpubr::stat_pvalue_manual(
    p_anno,
    label       = "label",
    xmin        = "xmin",
    xmax        = "xmax",
    y.position  = "y.position",
    tip.length  = 0,
    bracket_size = 0
  )

print(p_final)




##############
library(ggstatsplot)
library(ggplot2)
dt<-read.csv("CYPB1Bexpression.csv",sep = ",",check.names = T,row.names = 1)
p1 <- ggbetweenstats(
  data  = dt,
  x= group,
  y= expressions,
  plot.type = "box",#图标类型，不加参数则默认为boxviolin，其它可选为：box、violin
  title = ""
)
p2 <- ggbetweenstats(
  data  = dt,
  x     = group,
  y     = expressions,
  plot.type = "boxviolin",
  title = "expression level of CYP1B1",
  results.subtitle = TRUE,#决定是否将统计检验的结果显示为副标题（默认TRUE）;如果设置为FALSE,则仅返回绘图
  subtitle = NULL,#副标题,默认显示统计检验结果,自定义则results.subtitle=FALSE
  outlier.tagging = TRUE,#是否标记离群异常值，默认FALSE
  outlier.shape = 19,#异常值形状,可设置为NA将其隐藏（不是删除，因此不会影响统计检验结果）
  outlier.color = "pink",#异常值颜色
  outlier.label.args = list(size = 4),#异常值标签大小
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6),
                    alpha = 0.5, size = 3, stroke = 0),#传递给geom_point的参数设置
  violin.args = list(width = 0.4, alpha = 0.2),#传递给geom_violin的参数设置
  ggtheme = theme_bw(),#主题修改，可直接调用ggplot2的主题，默认主题为ggstatsplot::theme_ggstatsplot()
  package = "ggsci",#提取调色板所需的包
  palette = "uniform_startrek"#选择提取包中的调色板
)
p2

dt<-read.csv("mtnernl expression.csv",sep = ",",check.names = T,row.names = 1)
p1 <- ggbetweenstats(
  data  = dt,
  x= group,
  y= expressions,
  plot.type = "box",#图标类型，不加参数则默认为boxviolin，其它可选为：box、violin
  title = ""
)
p2 <- ggbetweenstats(
  data  = dt,
  x     = group,
  y     = expressions,
  plot.type = "boxviolin",
  title = "expression level of METRNL",
  results.subtitle = TRUE,#决定是否将统计检验的结果显示为副标题（默认TRUE）;如果设置为FALSE,则仅返回绘图
  subtitle = NULL,#副标题,默认显示统计检验结果,自定义则results.subtitle=FALSE
  outlier.tagging = TRUE,#是否标记离群异常值，默认FALSE
  outlier.shape = 19,#异常值形状,可设置为NA将其隐藏（不是删除，因此不会影响统计检验结果）
  outlier.color = "pink",#异常值颜色
  outlier.label.args = list(size = 4),#异常值标签大小
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6),
                    alpha = 0.5, size = 3, stroke = 0),#传递给geom_point的参数设置
  violin.args = list(width = 0.4, alpha = 0.2),#传递给geom_violin的参数设置
  ggtheme = theme_bw(),#主题修改，可直接调用ggplot2的主题，默认主题为ggstatsplot::theme_ggstatsplot()
  package = "ggsci",#提取调色板所需的包
  palette = "uniform_startrek"#选择提取包中的调色板
)
p2

