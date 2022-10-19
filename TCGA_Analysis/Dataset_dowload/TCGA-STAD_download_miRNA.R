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

# 开始下载
library(TCGAbiolinks)
projects <- c("TCGA-STAD") # 选择STAD数据集
# miRNA Transcriptome
query <- GDCquery(
  project = projects,
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling",
)

GDCdownload(query,
            method = "api", 
            directory = getwd(),
            files.per.chunk  =  10)

# 数据转换
miRna <- GDCprepare(query, save = FALSE, 
                  summarizedExperiment = T, directory = getwd(), 
                  remove.files.prepared = F)


list <- seq(3,length(colnames(miRna)),3)

list <- c(1, list)
list  
  
miRNA_rpm <- miRna[,list]
rownames(miRNA_rpm) <- miRNA_rpm$miRNA_ID  

miRNA_rpm <- miRNA_rpm[,-1]
colnames(miRNA_rpm) <- sub('reads_per_million_miRNA_mapped_',"",colnames(miRNA_rpm))
colnames(miRNA_rpm) <- substr(colnames(miRNA_rpm), 
                                            start = 1,stop = 12)
colnames(miRNA_rpm) <- gsub("-",".",colnames(miRNA_rpm)) # 字符串替换
# write.csv(miRNA_rpm,file='TCGA-STAD_miRNA_expression.csv')

candidates <- c('hsa-mir-20a','hsa-mir-92a-1','hsa-mir-92a-2')
write.csv(miRNA_rpm[candidates,],file = 'miRNA_candidates.csv')
