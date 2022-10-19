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
library('stringr')

# 读取基因注释文件
Ginfo<-read.table("overlapTable27.txt",
                  header=T,sep="\t",check.names=F,row.names = 1)

# 开始下载
library(TCGAbiolinks)
projects <- c("TCGA-STAD")#选择STAD数据集
# RNA Transcriptome
query <- GDCquery(
  project = projects,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",#选择STAR-Count数据格式
  #barcode = ids
)

# 下载至当前目录
GDCdownload(query,
            method = "api", 
            directory = getwd(),
            files.per.chunk  =  10)

# 数据转换
rna <- GDCprepare(query, save = FALSE, 
                  summarizedExperiment = T, directory = getwd(), 
                  remove.files.prepared = F)
#expMatrix<-rna@assays@data@listData$tpm_unstrand
#expMatrix<-rna@assays@data@listData$fpkm_unstrand
expMatrix<-rna@assays@data@listData$unstranded
#上面选择了tpm_unstrand数据格式，可以更换为其他格式
# 用正则提取simple的ENSEMBLE ID
rownames(expMatrix)<-substr(str_extract(rna@rowRanges$gene_id, pattern = ".*\\."), 
                            1, 
                            nchar(str_extract(rna@rowRanges$gene_id, pattern = ".*\\."))-1)
colnames(expMatrix)<-rna@colData@listData$barcode
colnames(expMatrix) <- substr(colnames(expMatrix), 
                              start = 1,stop = 15)


# 基因注释
rownames(Ginfo)<-Ginfo$simple # 选择simple的ENSEMBLE ID
comgene <- intersect(rownames(expMatrix),rownames(Ginfo)) #取可注释基因
expr <- as.data.frame(expMatrix)[comgene,]
Ginfo <- Ginfo[comgene,]
expr$gene <- Ginfo$genename
expr <- expr[!duplicated(expr$gene),] #去重复
Ginfo <- Ginfo[rownames(expr),]
rownames(expr) <- expr$gene
expr <- expr[,-ncol(expr)] # ? 减多少表示删除第多少行

# 获取临床数据
clinicaldata<-as.data.frame(rna@colData)
clinicaldata_definition<-as.data.frame(clinicaldata)[,c('barcode','patient','sample','definition')]
clinicaldata_definition2<-as.data.frame(clinicaldata)
# 导出csv文件
write.csv(expr, file='STAD_Gene_unstranded.csv',quote=F,row.names=T)
write.csv(clinicaldata_definition, file='STAD_clinicaldata_definition.csv',row.names=T)
write.table(clinicaldata_definition2, file='STAD_clinicaldata_definition2.txt',sep =";",row.names=T) # csv的分隔符为‘,’，表格内容有此符号会出错
