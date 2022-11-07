# WGCNA其译为加权基因共表达网络分析。
# 该分析方法旨在寻找协同表达的基因模块(module)，并探索基因网络与关注的表型之间的关联关系，以及网络中的核心基因。
# 适用于复杂的数据模式，推荐5组(或者15个样品)以上的数据。
# 一般可应用的研究方向有：不同器官或组织类型发育调控、同一组织不同发育调控、非生物胁迫不同时间点应答、病原菌侵染后不同时间点应答。
# 第一步计算任意两个基因之间的相关系数（Person Coefficient）。
# WGCNA分析时采用相关系数加权值，即对基因相关系数取N次幂，使得网络中的基因之间的连接服从无尺度网络分布(scale-freenetworks)
# 第二步通过基因之间的相关系数构建分层聚类树，聚类树的不同分支代表不同的基因模块，不同颜色代表不同的模块。
# 基于基因的加权相关系数，将基因按照表达模式进行分类，将模式相似的基因归为一个模块。
# 这样就可以将几万个基因通过基因表达模式被分成了几十个模块，是一个提取归纳信息的过程。

# 参考文章
# https://cloud.tencent.com/developer/article/1516749

# 安装R包
# BiocManager::install("WGCNA")
# install.packages("WGCNA",type = 'binary')
# BiocManager::install("GEOquery")


# 1. 输入数据准备
# 主要是表达矩阵，如果是芯片数据，那么常规的归一化矩阵即可
# 如果是转录组数据，最好是RPKM/TPM值或者其它归一化好的表达量
# 然后就是临床信息或者其它表型，总之就是样本的属性
# 为了保证后续脚本的统一性，表达矩阵统一用datExpr标识，临床 信息统一用datTraits标识

# 56 breast cancer cell lines were profiled to identify patterns of gene expression associated with subtype and response to therapeutic compounds.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48213
# wget -c ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48213/suppl/GSE48213_RAW.tar
# tar -xf GSE48213_RAW.tar
#g zip -d *.gz
# 提示文件有损坏，只有22种细胞信息
# 首先在GSE48213_RAW目录里面生成tmp.txt文件，使用shell脚本
# awk '{print FILENAME"\t"$0}' GSM*.txt |grep -v EnsEMBL_Gene_ID >tmp.txt

# 然后把tmp.txt导入R语言里面用reshape2处理即可
a=read.table('tmp.txt',sep = '\t',stringsAsFactors = F)
library(reshape2)
# dcast将长格式数据转换成宽格式数据,并且输出结果有行标签
# formula：形如x ~ y,x为行标签，y为列标签
fpkm <- dcast(a,formula = V2~V1) 
rownames(fpkm)=fpkm[,1]
fpkm=fpkm[,-1]
# sapply一次性对一堆数据执行某个函数
# strsplit用于分割字符串
colnames(fpkm)=sapply(colnames(fpkm),function(x) strsplit(x,"_")[[1]][1])

library(GEOquery)
a=getGEO('GSE48213')
# pData: Combine data, type, comments, and metadata information to create a new pdata object, or check such an object for consistency
metadata_row=pData(a[[1]])
metadata=metadata_row[,c(2,10,12)]
#  trimws() 函数用于修剪前导和尾随空格
datTraits = data.frame(gsm=metadata[,1],
                       cellline=trimws(sapply(as.character(metadata$characteristics_ch1),function(x) strsplit(x,":")[[1]][2])),
                       subtype=trimws(sapply(as.character(metadata$characteristics_ch1.2),function(x) strsplit(x,":")[[1]][2]))
)
save(fpkm,datTraits,file = 'GSE48213-wgcna-input.RData')

load('GSE48213-wgcna-input.RData')
library(WGCNA)
RNAseq_voom <- fpkm 
## 因为WGCNA针对的是基因进行聚类，而一般我们的聚类是针对样本用hclust即可，所以这个时候需要转置。
# apply函数经常用来计算矩阵中行或列的均值、和值的函数,1代表每一行
# mad () 用于计算中位数绝对偏差，即与中位数的绝对偏差的 (lo-/hi-)中位数。
WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:5000],])
datExpr0 <- WGCNA_matrix  ## top 5000 mad genes
datExpr <- datExpr0 

## 下面主要是为了防止临床表型与样本名字对不上，下载的数据损坏，只有22种细胞信息
sampleNames = rownames(datExpr);
traitRows = match(sampleNames, datTraits$gsm)
datTraits=datTraits[traitRows, ] 
rownames(datTraits) = datTraits[, 1] 


# 2. 确定最佳beta值
# Call the network topology analysis function
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# 设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
# 最佳的beta值就是sft$powerEstimate
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
##sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
# sign()用于查找作为参数提供的数字向量的元素的符号，返回：1 代表正数，-1 代表负数
# 纵轴代表对应的网络中log(k)与log(p(k))相关系数的平方。相关系数的平方越高，说明该网络越逼近无网路尺度的分布。
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
# 纵轴代表对应的基因模块中所有基因邻接函数的均值
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# 3.一步法构建共表达矩阵
# 把输入的表达矩阵的几千个基因组归类成了几十个模块
# 计算基因间的邻接性，根据邻接性计算基因间的相似性，然后推出基因间的相异性系数，并据此得到基因间的系统聚类树。
# 然后按照混合动态剪切树的标准，设置每个基因模块最少的基因数目为30。
# 根据动态剪切法确定基因模块后，再次分析，依次计算每个模块的特征向量值，然后对模块进行聚类分析，将距离较近的模块合并为新的模块。
net = blockwiseModules(
        datExpr,
        power = sft$powerEstimate,
        maxBlockSize = 6000,
        TOMType = "unsigned", minModuleSize = 30,
        reassignThreshold = 0, mergeCutHeight = 0.25,
        numericLabels = TRUE, pamRespectsDendro = FALSE,
        saveTOMs = F, 
        verbose = 3
)
table(net$colors) 


# 4. 模块可视化
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
# plotDendroAndColors函数，它接受一个聚类的对象，以及该对象里面包含的所有个体所对应的颜色。
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
## assign all of the gene to their corresponding module 
## hclust for the genes.

# 此聚类的代码可以不运行，跟WGCNA本身关系不大。
#明确样本数和基因数
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#首先针对样本做个系统聚类树
datExpr_tree<-hclust(dist(datExpr), method = "average")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)
## 如果这个时候样本是有性状，或者临床表型的，可以加进去看看是否聚类合理
#针对前面构造的样品矩阵添加对应颜色
sample_colors <- numbers2colors(as.numeric(factor(datTraits$subtype)), 
                                colors = c("white","blue","red","green"),signed = FALSE)
## 这个给样品添加对应颜色的代码需要自行修改以适应自己的数据分析项目。
#  sample_colors <- numbers2colors( datTraits ,signed = FALSE)
## 如果样品有多种分类情况，而且 datTraits 里面都是分类信息，那么可以直接用上面代码，当然，这样给的颜色不明显，意义不大。
#构造10个样品的系统聚类树及性状热图
par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(datExpr_tree, sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")

# 5. 模块和性状的关系
## 这一步主要是针对于连续变量，如果是分类变量，需要转换成连续变量方可使用
table(datTraits$subtype)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
design=model.matrix(~0+ datTraits$subtype)
colnames(design)=levels(as.factor(datTraits$subtype))
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩 (样本vs模块)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
png("step5-Module-trait-relationships.png",width = 800,height = 1200,res = 120)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
# 因为一些历史遗留问题，这个图片缺乏X轴的标记
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# 除了上面的热图展现形状与基因模块的相关性外
# 还可以是条形图,但是只能是指定某个形状
# 或者自己循环一下批量出图。
Luminal = as.data.frame(design[,3]);
names(Luminal) = "Luminal"
y=Luminal
GS1=as.numeric(cor(y,datExpr, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance,
                          moduleColors, mean, na.rm=T)
sizeGrWindow(8,7)
par(mfrow = c(1,1))
# 如果模块太多，下面的展示就不友好
# 不过，我们可以自定义出图。
plotModuleSignificance(GeneSignificance,moduleColors)
dev.off()

# 6. 感兴趣性状的模块的具体基因分析
# 首先计算模块与基因的相关性矩阵 
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
## 算出每个模块跟基因的皮尔森相关系数矩阵
## MEs是每个模块在每个样本里面的值
## datExpr是每个基因在每个样本的表达量
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

# 再计算性状与基因的相关性矩阵 
## 只有连续型性状才能只有计算
## 这里把是否属于 Luminal 表型这个变量用0,1进行数值化。
Luminal = as.data.frame(design[,3]);
names(Luminal) = "Luminal"
geneTraitSignificance = as.data.frame(cor(datExpr, Luminal, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Luminal), sep="");
names(GSPvalue) = paste("p.GS.", names(Luminal), sep="");

#最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析 
module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Luminal",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# 7. 网络的可视化












