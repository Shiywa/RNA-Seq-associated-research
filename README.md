# RNA-Seq-associated-research
general analysis for all RNA-Seq

# 目录

- [从GEO数据库中根据GEO编号获取micro Array矩阵及注释数据](#micro-Array)
- [基因集打分策略](#基因集打分策略)
  - [使用VISION](#VISION)
  - [使用GSVA](#GSVA&ssGSEA)
  
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

简单来说，**VISION**通过预先的降维（例如PCA）构建细胞在多维空间的相似性图，然后基于**Geary-C value**评估在多维空间的相似性图上是否存在特征一致性，Geary-C value越大代表对应特征集存在一组特定基因在细胞的多维空间的相似性图具有一致性表达，提示该特定功能可能有特定激活的趋势。
这个方法理论上可以帮助挖掘新的细胞类型，且可应用到空间组学的研究。
在确定了对应特征集在多维空间具有特异性活性后，基于细胞注释信息，使用**AUCell**评估不同细胞类型中该特征集的AUC score。

```
library(VISION)
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
由于**VISION**本身基于web开发了交互式的可视化工具，但是没有便捷的绘图函数，所以需要自己编写画图的函数。

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
<img width="667" alt="image" src="https://github.com/user-attachments/assets/b7acdce7-a43c-4441-a4a9-883310e4df3f" />

### GSVA&ssGSEA

两种方法都被包装在`GSVA`这个包里面了，这两种方法的思想是非常相似的，都是给予单个样本里面的基因表达rank来进行富集分析，但是两种方法所用到的累积函数是不一样的，此外GSVA方法依赖其他样本的表达情况
ssGSEA则更加独立。

新的`GSVA`更新了对于这两种方法的使用方式：

1. GSVA
```
# 读取基因集
gene_sets <- getGmt("KEGG_metabolism_pathway.gmt")

# 提取基因集名称和对应基因转换为list
gene_sets_list <- lapply(gene_sets, geneIds)

# 设置 list 的名称为基因集名称
names(gene_sets_list) <- names(gene_sets)

# 设置并行参数，使用所有可用 CPU 核心
bpparam <- MulticoreParam(workers = parallel::detectCores() - 20)  

# 构建input，设置所用参数，使用gsva方法分析
gsvapar <- gsvaParam(mtx, gene_sets_list,
                     kcdf = "Gaussian",
                     maxDiff = T, 
                     absRanking=FALSE,
                     tau = 1,
                     sparse = T)

# 运行 GSVA
gsva_results <- gsva(
  gsvapar,BPPARAM = bpparam
)

```
2. ssGSEA

```
ssgseapar <- ssgseaParam(mtx, gene_sets_list)

ssgsea_results <- gsva(
  ssgseapar,BPPARAM = bpparam
)
```

对两种方法进行比较可以发现，**GSVA**更易受到低表达基因的影响，此外**GSVA**由于依赖所有样本的数据，对于批次效应更加敏感，所以**ssGSEA**更加适合单细胞多样本数据，当然两种方法可以整合起来一起分析。

### 通路核心基因挖掘

确定通路活性以后往往希望挖掘是否存在critical feature对于通路有重要作用，目前比较基本的思路就是：
1. 跟通路活性评分高相关的feature；
2. 基底表达量比较高的feature。（更易于解释）

热图方便展示，以下是我写的一个展示潜在重要基因的热图
```
# 首先是SCP中的PrepareDB函数获得小鼠的KEGG数据
library(SCP)
db_list <- PrepareDB(species = "Mus_musculus", db = c("KEGG","Reactome"))
db_KEGG <- db_list$Mus_musculus$KEGG

# 然后获得通路全部基因
ge1 <- db_KEGG$TERM2GENE$symbol[which(db_KEGG$TERM2GENE$Term == db_KEGG$TERM2NAME$Term[which(db_KEGG$TERM2NAME$Name == i)])]

# 计算各个通路平均表达量以绘制热图主体
library(Seurat)
avg_exp_matrix <- AverageExpression(mock.ctrl.mac.obj, return.seurat = FALSE,assays = "RNA",features = ge1,group.by = "deve_anno")
avg_exp_matrix <- avg_exp_matrix$RNA
avg_exp_matrix <- avg_exp_matrix[rowSums(avg_exp_matrix) != 0,]
avg_exp_matrix2 <- t(scale(t(avg_exp_matrix)))

# 计算每个基因与三种富集方法得到的通路活性的相关性，这里用了spearman，不限于线性相关  
cor_results1 <- apply(vis@exprData[rownames(avg_exp_matrix),], 1, function(row) cor(row, VISION_score[,i], use = "pairwise.complete.obs", method = "spearman"))
cor_results2 <- apply(vis@exprData[rownames(avg_exp_matrix),], 1, function(row) cor(row, ssGSEA_score[,i], use = "pairwise.complete.obs", method = "spearman"))
cor_results3 <- apply(vis@exprData[rownames(avg_exp_matrix),], 1, function(row) cor(row, GSVA_score[,i], use = "pairwise.complete.obs", method = "spearman"))

# 计算高相关基因，我这里是经验数值，cor > 0.3且cor位于top 25%  
high_cor_genes1 <- names(which(cor_results1 > quantile(cor_results1, 0.75) & cor_results1 > 0.3))  # 取前25%高AUC细胞
high_cor_genes2 <- names(which(cor_results2 > quantile(cor_results2, 0.75) & cor_results2 > 0.3))  # 取前25%高AUC细胞
high_cor_genes3 <- names(which(cor_results3 > quantile(cor_results3, 0.75) & cor_results3 > 0.3))  # 取前25%高AUC细胞
high_cor_genes <- unique(c(high_cor_genes1,high_cor_genes2,high_cor_genes3))

# 构建热图左注释，包括label、平均表达量以及相关性三部分  
left_anno = rowAnnotation(label = anno_mark(at= match(high_cor_genes,rownames(avg_exp_matrix)),labels = high_cor_genes,side = "left"),
                            
                            avg1 = anno_barplot(avg_exp_matrix[,c("Resting Conv PM")],gp = gpar(fill = new.color3["Resting Conv PM"]),
                                                axis_param = list(direction = "reverse"),width = unit(1, "cm")),
                            avg2 = anno_barplot(avg_exp_matrix[,c("Resting Mono-like SPM")],gp = gpar(fill = new.color3["Resting Mono-like SPM"]),
                                                axis_param = list(direction = "reverse"),width = unit(1, "cm")),
                            avg3 = anno_barplot(avg_exp_matrix[,c("Resting GLPM")],gp = gpar(fill = new.color3["Resting GLPM"]),
                                                axis_param = list(direction = "reverse"),width = unit(1, "cm")),
                            VISION = anno_barplot(cor_results1, baseline = 0,axis_param = list(direction = "reverse"),width = unit(2, "cm"),
                                               gp = gpar(fill = "#1B9E77")),
                            ssGSEA = anno_barplot(cor_results2, baseline = 0,axis_param = list(direction = "reverse"),width = unit(2, "cm"),
                                               gp = gpar(fill = "#D95F02")),
                            GSVA = anno_barplot(cor_results3, baseline = 0,axis_param = list(direction = "reverse"),width = unit(2, "cm"),
                                               gp = gpar(fill = "#7570B3")))
  
  # 热图绘制
p <- Heatmap(avg_exp_matrix2,cluster_rows = T,cluster_columns = T,name = "zscore",show_row_dend = F,show_column_dend = F,
               col = colorRamp2(c(-2,0,2),c("#3690C0","white","#fa5252")),top_annotation = top_anno,
               left_annotation = left_anno,
               row_names_side = "left",border = T,row_names_max_width = unit(4, "cm"),show_row_names = F,
               width = unit(dim(avg_exp_matrix)[2]*0.8, "cm"), height = unit(dim(avg_exp_matrix)[1]*0.3, "cm"))

lgd_title = Legend(pch = " ", type = "points", labels = i)

# 绘制的时候控制好热图的边界距离，同时添加一个pathway label
pdf("plot.pdf",height = 15,width = 15)
draw(p, heatmap_legend_side = "right", padding = unit(c(4, 5, 2, 2), "mm"),annotation_legend_list = list(lgd_title))
dev.off()
```
<img width="925" alt="image" src="https://github.com/user-attachments/assets/fb21830d-0d51-481f-86f5-eeb4e366dffb" />



