# RNA-Seq-associated-research
general analysis for all RNA-Seq

# 目录

- [从GEO数据库中根据GEO编号获取micro Array矩阵及注释数据](#micro-Array)
- [基因集打分策略](#基因集打分策略)
  - [使用VISION](#VISION)
  - [使用GSVA](#micro-Array)
  
## micro-Array

### 将GEO的micro array文件转换为matrix文件

~~如果网络不好可以直接下载seris里的matrix文件到本地，然后自动加载。注意matrix文件要放在当前工作目录下。~~
原来并不需要提前下载这个文件，有GEO号就可以自动下载到`destdir`设定的路径

```
library(GEOquery)
gset <- getGEO("GSE80697", GSEMatrix =TRUE, getGPL=FALSE,destdir = ".",AnnotGPL = F)
```

第一步提取样本分组信息，`gset`中包含了meta信息，使用`pData`可以读取。

```
meta <- pData(gset[[1]])
```

第二步提取基因表达矩阵，`gset`中包含了matrix信息，使用`Biobase::exprs()`可以转换为矩阵格式，并利用`limma::normalizeBetweenArrays()`进行矩阵矫正。

```
gset2 <- gset[[1]]
ex <- exprs(gset2)
ex=normalizeBetweenArrays(ex) 
boxplot(ex,outline=FALSE, notch=T,col=factor(meta$`virus:ch1`), las=2) ## 查看表达水平
```

第三步将探针名称转换为基因名，首先`gset`中可以提取该芯片是哪个平台的，根据平台号，寻找对应的**symbol-ID**对应信息所属的包，一般在`biomanager`上面可以下载到。

具体查询平台号与对应包的关系，可以关注[jmzeng1314/idmap1](https://github.com/jmzeng1314/idmap1)，有详细的对应信息以及他开发的小工具。

加载了对应包之后，找到末尾是由**SYMOBL**结尾的文件，使用`BiocGenerics::toTable()`可以将其转换为**symbol-ID**的data frame。

```
index = gset[[1]]@annotation
library(HsAgilentDesign026652.db) # plat GPL13497
ids <- toTable(HsAgilentDesign026652SYMBOL)
```

这三步搞完以后，将探针文件和矩阵文件合并，根据基因名去重，生成最终矩阵文件

```
ex <- as.data.frame(ex) 

exp <- ex %>% mutate(probe_id=rownames(ex)) %>% inner_join(ids,by="probe_id") %>% 
  dplyr::select(probe_id, symbol,everything()) ### everything() 加括号
exp <- exp[!duplicated(exp$symbol),] 
rownames(exp) <- exp$symbol 
exp <- exp[,-(1:2)]
```

有一种情况，可能因为微阵列的数据比较久远，对应平台的注释文件在`biomanager`下载不到，但是GEO数据库本身其实提供了微阵列数据的注释信息。
同样根据`index`得到对应平台号，然后使用`getGEO`获得对应注释信息，同时将注释信息和矩阵合并。要注意的是这种文件中还是存在很多**id**对应不到**symbol**的，要把他们去除。

```
index = gset[[1]]@annotation
gpl <- getGEO(index)
platform <- Table(gpl)
platform$ID <- as.character(platform$ID)

ex <- as.data.frame(ex) 

exp <- ex %>% mutate(ID=rownames(ex)) %>% inner_join(platform,by="ID") %>% 
  dplyr::select(ID, GENE_SYMBOL,everything()) ### everything() 加括号
exp <- exp[!duplicated(exp$GENE_SYMBOL),] 
exp <- exp[-which(exp$SPOT_ID == "empty"),]
rownames(exp) <- exp$GENE_SYMBOL 
exp <- exp[,grep("^GSM",colnames(exp))]

```

## 基因集打分策略

无论在bulkRNA-seq还是scRNA-seq中，基于基因的表达进行特定功能基因集的打分都是必不可少的，以下总结了几种打分的方法。

### VISION 

[VISION原文](https://www.nature.com/articles/s41467-019-12235-0)以及[VISION github 页面](https://github.com/YosefLab/VISION)详细介绍了相应的功能及原理。

简单来说，VISION通过预先的降维（例如PCA）构建细胞在多维空间的相似性图，然后基于Geary-C value评估在多维空间的相似性图上是否存在特征一致性，Geary-C value越大代表对应特征集存在一组特定基因在细胞的多维空间的相似性图具有一致性表达，提示该特定功能可能有特定激活的趋势。
这个方法理论上可以帮助挖掘新的细胞类型，且可应用到空间组学的研究。
在确定了对应特征集在多维空间具有特异性活性后，基于细胞注释信息，使用**AUCell**评估不同细胞类型中该特征集的AUC score。

```
count <- Assay.obj@assays$RNA@counts
n.umi <- colSums(count)
scaled_counts <- t(t(count) / n.umi) * median(n.umi)

meta <- Assay.obj@meta.data

table(mock.ctrl.mac.obj$deve_anno)/dim(mock.ctrl.mac.obj@assays$RNA)[2]
## 鉴于macrophage中最小的群体所占比例为0.5%，所以默认的预制0.1%是可取的

vis <- Vision(scaled_counts,
              signatures = c("KEGG_metabolism_pathway.gmt"), # 需要记忆集gmt文件
              meta = meta[,c("deve_anno","State")])
options(mc.cores = 96)

# 一步法分析
vis <- analyze(vis,hotspot =F) # hotspot是空间数据专用参数

vis <- addProjection(vis, "UMAP", Assay.obj@reductions$harmony_umap@cell.embeddings)

saveRDS(vis,file = "VISION_result.rds")
```
由于VISION本身基于web开发了交互式的可视化工具，但是没有便捷的绘图函数，所以需要自己编写画图的函数。

1. 从结果中提取数据,包括Geary-C value，AUC value以及AUC对应的FDR：
```
sigScores <- getSignatureScores(vis)

full_C_value <- vis@LocalAutocorrelation$Signatures

full_C_value <- full_C_value[order(full_C_value[,1],decreasing = T),]

full_AUC <- do.call(cbind, lapply(vis@ClusterComparisons$Signatures$deve_anno, `[[`, 1))
rownames(full_AUC) <- rownames(vis@ClusterComparisons$Signatures$deve_anno[[1]])

full_AUC_FDR <- do.call(cbind, lapply(vis@ClusterComparisons$Signatures$deve_anno, `[[`, 3))
rownames(full_AUC_FDR) <- rownames(vis@ClusterComparisons$Signatures$deve_anno[[1]])

full_AUC <- full_AUC[rownames(full_C_value),]
full_AUC_FDR <- full_AUC_FDR[rownames(full_C_value),]
```

2. 热图的绘制：
```
library(ComplexHeatmap)
library(circlize)

# 标注特定的行名，可自定义
label_row <- c("Arachidonic acid metabolism","Linoleic acid metabolism","Pyrimidine metabolism",
               "Glycosphingolipid biosynthesis - lacto and neolacto series","Tryptophan metabolism",
               "Inositol phosphate metabolism","Arginine biosynthesis","Nicotinate and nicotinamide metabolism","Sphingolipid metabolism","Oxidative phosphorylation",
               "Cholesterol metabolism","Steroid biosynthesis","Glycosphingolipid biosynthesis - globo and isoglobo series")
label_row_index <- match(label_row,rownames(full_AUC),)

p <- Heatmap(full_AUC,cluster_rows = T,name = "AUC",show_row_dend = F,show_column_dend = F,
             col = colorRamp2(c(0,0.5,1),c("#3690C0","white","#fa5252")), # 非常漂亮的热图color
             # 左注释展示Geary-C value以及特定行名
             left_annotation = rowAnnotation( label = anno_mark(at= label_row_index,labels = label_row,side = "left"),
                                              C_value = anno_simple(full_C_value[,1],
                                                                   #pch = as.numeric(as.character(round(full_C_value[,1],2))),
                                                                   col = colorRamp2(c(0, 0.5, 1), c("#3690C0","white","#fa5252")),
                                                                   border = T, width = unit(1.5, "cm"),height = unit(2, "cm"))
                                            ),
             row_names_side = "left",border = T,row_names_max_width = unit(4, "cm"),show_row_names = F,
             # 热图单元格展示FDR
             cell_fun = function(j, i, x, y, width, height, fill) {
               if(full_AUC_FDR[i, j] < 0.05)
                 grid.text("*", x, y, gp = gpar(fontsize = 12))
             })

p
```
<img width="705" alt="image" src="https://github.com/user-attachments/assets/753f05fb-1d10-46b6-b6c6-887afda61b90" />

[pathway_AUC_C_value_Heatmap2.pdf](https://github.com/user-attachments/files/19001200/pathway_AUC_C_value_Heatmap2.pdf)






