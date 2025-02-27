# RNA-Seq-associated-research
general analysis for all RNA-Seq

# 目录

- [从GEO数据库中根据GEO编号获取micro Array矩阵及注释数据](#micro-Array)
- [基因集打分策略](#micro-Array)
 - [使用VISION](#micro-Array)
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








