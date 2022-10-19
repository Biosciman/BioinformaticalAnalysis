# 加载R包
library('devtools')
library('TCGAbiolinksGUI.data')
library('AnnotationHub')
library('ExperimentHub')
library('SummarizedExperiment')
library('GenomeInfoDbData')
library('GenomicDataCommons')
library('TCGAbiolinks')
library('limma')
library('edgeR')
library('statmod')
library('GO.db')
library('org.Hs.eg.db')

### unstanded方法
# 导入csv文件
expr <- read.csv(file='STAD_Gene_unstranded.csv',
      header=T)
rownames(expr) <- expr$X
expr <- expr[, -1]

clinicaldata_definition<-read.csv(file='STAD_clinicaldata_definition.csv',
                                  header=T)
rownames(clinicaldata_definition) <- clinicaldata_definition$X
clinicaldata_definition <- clinicaldata_definition[, -1]
rownames(clinicaldata_definition) <- substr(rownames(clinicaldata_definition), 
                         start = 1,stop = 15)
rownames(clinicaldata_definition) <- gsub("-",".",rownames(clinicaldata_definition)) # 字符串替换

# 分析差异基因——数据准备
group <- clinicaldata_definition$definition # 选取临床信息中的definition（正常组织or癌组织）
group <- gsub(' ','_',group) # 字符串替换，将空格替换为_，防止代码bug
dgelist <- DGEList(counts = expr, group = group) # 构建DGEList对象
dgelist_norm <- calcNormFactors(dgelist) # calcNormFactors()函数对数据标准化，以消除由于样品制备或建库测序过程中带来的影响，method = "TMM"
cutoff <- 1 # 设置阈值，这个阈值对于大多数数据集都能很好地分隔表达的基因与不表达的基因，CPM值为1意味着如果一个基因在测序深度最低的样品中有至少20个计数
drop <- which(apply(cpm(dgelist_norm), 1, max) < cutoff) # 过滤低表达，根据CPM值(counts per million reads)进行过滤
dgelist_norm_keep <- dgelist_norm[-drop,] # 保留高表达

# 分析差异基因——limma分析数据
mm <- model.matrix(~0 + group) # 分组，是为1，不是为0。自变量为是否为癌症，因变量为基因表达。
par(mfrow = c(1, 3)) # 实现一页多图，1行，3列
# voom 方法估计对数计数的均值-方差关系，为每个观测值生成一个精确权重，并将这些权重输入到 limma 经验贝叶斯分析管道
# 假设log-CPM值符合正态分布，并使用由voom函数计算得到的精确权重来调整均值与方差的关系，从而对log-CPM值进行线性建模
# 如果横坐标接近0的位置出现迅速上升，说明low counts数比较多
# 得到的曲线不平滑表示数据需要再进行过滤，比如提高CPM的cutoff
# 每个黑点表示一个基因，红线为对这些点的拟合
# https://doi.org/10.1186/gb-2014-15-2-r29
y <- voom(dgelist_norm_keep, mm, plot = T) 
# lmFit使用每个基因的加权最小二乘拟合线性模型
# lmFit为每个基因的表达值拟合一个模型。然后，通过利用所有基因的信息来进行经验贝叶斯调整，这样可以获得更精确的基因水平的变异程度估计
fit <- lmFit(y, mm)

# makeContrasts实际就是定义比较分组信息
# coef 是一个通用函数，它从建模函数返回的对象中提取模型系数。coefficients是它的别名。
contr <- makeContrasts(groupPrimary_solid_Tumor - groupSolid_Tissue_Normal, levels = colnames(coef(fit)))
# 比较每个基因
tmp <- contrasts.fit(fit, contr)
# 再利用经验贝叶斯（Empirical Bayes），将比其他基因大得多或小得多的标准误差缩小到平均标准误差
# https://doi.org/10.2202/1544-6115.1027
tmp <- eBayes(tmp)

# 使用plotSA 绘制了log2残差标准差与log-CPM均值的关系。平均log2残差标准差由水平蓝线标出
plotSA(tmp, main="Final model: Mean-variance trend")

# 基因表达上调或下调总览
summary(decideTests(tmp))
# topTable 列出差异显著基因
# logFC: log2 fold change of Primary_solid_Tumor/Solid_Tissue_Normal
# AveExpr: Average expression across all samples, in log2 CPM
# t: logFC divided by its standard error
# P.Value: Raw p-value (based on t) from test that logFC differs from 0
# adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
# B: log-odds that gene is DE (arguably less useful than the other columns)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

# 输出文件
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.csv(top.table, file = "Primary_solid_TumorvsSolid_Tissue_Normal_fpkm_unstrand.csv", row.names = F, quote = F)



### FPKM or TPM方法
# 导入csv文件
expr <- read.csv(file='STAD_Gene_unstranded.csv',
                 header=T)
rownames(expr) <- expr$X
expr <- expr[, -1]

clinicaldata_definition<-read.csv(file='STAD_clinicaldata_definition.csv',
                                  header=T)
rownames(clinicaldata_definition) <- clinicaldata_definition$X
clinicaldata_definition <- clinicaldata_definition[, -1]
rownames(clinicaldata_definition) <- substr(rownames(clinicaldata_definition), 
                                            start = 1,stop = 15)
rownames(clinicaldata_definition) <- gsub("-",".",rownames(clinicaldata_definition)) # 字符串替换

# 分析差异基因——数据准备
group <- clinicaldata_definition$definition # 选取临床信息中的definition（正常组织or癌组织）
group <- gsub(' ','_',group) # 字符串替换，将空格替换为_，防止代码bug
dgelist <- DGEList(counts = expr, group = group) # 构建DGEList对象

logCPM <- cpm(dgelist, log=TRUE, prior.count=3)
mm <- model.matrix(~0 + group) 
fit <- lmFit(logCPM, mm)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))


# makeContrasts实际就是定义比较分组信息
# coef 是一个通用函数，它从建模函数返回的对象中提取模型系数。coefficients是它的别名。
contr <- makeContrasts(groupPrimary_solid_Tumor - groupSolid_Tissue_Normal, levels = colnames(coef(fit)))
# 比较每个基因
tmp <- contrasts.fit(fit, contr)
# 再利用经验贝叶斯（Empirical Bayes），将比其他基因大得多或小得多的标准误差缩小到平均标准误差
# https://doi.org/10.2202/1544-6115.1027
tmp <- eBayes(tmp)

# 使用plotSA 绘制了log2残差标准差与log-CPM均值的关系。平均log2残差标准差由水平蓝线标出
plotSA(tmp, main="Final model: Mean-variance trend")

# 基因表达上调或下调总览
summary(decideTests(tmp))
# topTable 列出差异显著基因
# logFC: log2 fold change of Primary_solid_Tumor/Solid_Tissue_Normal
# AveExpr: Average expression across all samples, in log2 CPM
# t: logFC divided by its standard error
# P.Value: Raw p-value (based on t) from test that logFC differs from 0
# adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
# B: log-odds that gene is DE (arguably less useful than the other columns)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

# 输出文件
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.csv(top.table, file = "Primary_solid_TumorvsSolid_Tissue_Normal_fpkm_unstrand_limma-trend.csv", row.names = F, quote = F)


## 可视化
# 热图
library(pheatmap)
library(RColorBrewer)
color <- colorRampPalette(c('green','white','red'))(100) # 定义热图颜色
cg = rownames(top.table[which(top.table$adj.P.Val<0.05),]) # 提取差异表达基因
dim(top.table)
mat = expr[cg,]
n=t(scale(t(mat))) # 使用scale方法来对数据进行中心化和标准化，t() 函数行列转换
n[n>1]=1 # 高表达定义为1
n[n<1]=-1 # 低表达低于为-1
ac=data.frame(group=group)
rownames(ac) = colnames(mat)
ht <- pheatmap(n,show_rownames=F,show_colnames=F,
               cluster_rows=F,cluster_cols=T,
               annotation_col=ac, color=color
               )

# 火山图
dev.off()
library(EnhancedVolcano)
library(airway)
# pCutoff和FCcutoff可设置阈值
v = EnhancedVolcano(top.table,
                    lab =rownames(top.table),
                    x = 'logFC',
                    y = 'P.Value',
                    xlim = c(-5, 5),
                    title = 'Differential expression',
                    col=c('black','blue','green','red1'),
                    colAlpha = 1,

)

plot(v)
