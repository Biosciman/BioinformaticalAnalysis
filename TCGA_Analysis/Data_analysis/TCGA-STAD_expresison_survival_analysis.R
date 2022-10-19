library(stringr)
library(survival)
library(survminer)
# 导入csv文件
expr <- read.csv(file='STAD_Gene_unstranded.csv',
                 header=T)
rownames(expr) <- expr$X
expr <- expr[, -1]

clinicaldata_definition<-read.table(file='STAD_clinicaldata_definition2.txt',
                                    sep=';',
                                  header=T) # 写入临床数据的分隔符为；
rownames(clinicaldata_definition) <- clinicaldata_definition$X
clinicaldata_definition <- clinicaldata_definition[, -1]
rownames(clinicaldata_definition) <- substr(rownames(clinicaldata_definition), 
                                            start = 1,stop = 15)
rownames(clinicaldata_definition) <- gsub("-",".",rownames(clinicaldata_definition)) # 字符串替换

# 截取肿瘤组织表达数据
exprTumor <- expr[, which(clinicaldata_definition$definition=='Primary solid Tumor')]
summary(clinicaldata_definition$definition=='Primary solid Tumor') # 肿瘤样品数总览
dim(expr)
dim(exprTumor)

# 整理临床数据
meta = clinicaldata_definition[,c( 'sample_submitter_id', 
                    'vital_status', 
                    'days_to_death', 
                    'days_to_last_follow_up', 
                    'race', 
                    'days_to_birth', 
                    'gender' , 
                    'paper_TNM.Stage' ,
                    'definition')]

colnames(meta)=c('ID','event','death','last_followup','race','age','gender','stage', 'definition') # 简化meta的列名
meta[meta==""]=NA # 空着的值改为NA
meta[meta=="X"]=NA # stage中X值改为NA
meta <- meta[which(meta$definition=='Primary solid Tumor'),]
dim(meta)
rownames(meta) <- substr(meta$ID, 
                                            start = 1,stop = 15)
rownames(meta) <- gsub("-",".",rownames(meta)) # 字符串替换

# 整理生存分析的输入数据
# 由随访时间和死亡时间计算生存时间(月)
table(meta$event) 
meta$time = ifelse(meta$event=="Alive", meta$last_followup, meta$death)
meta$time = as.numeric(meta$time)/30

# 根据生死定义event，活着是0，死的是1
meta$event=ifelse(meta$event=='Alive', 0, 1)
table(meta$event)

# 年龄和年龄分组
meta$age=ceiling(abs(as.numeric(meta$age))/365) # ceiling将返回最接近输入值但大于输入值的值。
meta$age_group=ifelse(meta$age>median(meta$age,na.rm = T),'older','younger')
table(meta$age_group)

# stage
table(meta$stage)
a = str_extract_all(meta$stage,"I|V") # 提取分期信息
b = sapply(a,paste,collapse = "") # 将a中的信息无缝连接
meta$stage = b

# 去掉生存信息不全或者生存时间小于0.1(月)的病人，样本纳排标准不唯一，且差别很大
k1 = meta$time>=0.1
table(k1)
k2 = !(is.na(meta$time)|is.na(meta$event)) # is.na()用于检测确实值是否缺失，缺失返回True，未缺失返回False
table(k2)

exprTumor = exprTumor[,k1&k2] # 取k1和k2均为true的column
meta = meta[k1&k2,] # 取k1和k2均为true的row

# 绘制生存曲线
# 格式：
# sfit <- survfit(Surv(time,event)~group,data= meta)
# time指总生存期
# event是终点事件
# group是分组，是临床信息中的一列
# data是临床信息表格

#年龄
sfit <- survfit(Surv(time, event)~age_group, data=meta)
ggsurvplot(sfit, conf.int=F, pval=T)

# 曲线美化
ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)

# 单个基因的生存曲线
dev.off() # dev.off()删除当前图片信息
g = rownames(exprTumor)[1]
meta$gene = ifelse(as.numeric(exprTumor[g,])> median(as.numeric(exprTumor[g,])),'high','low')
sfit1=survfit(Surv(time, event)~gene, data=meta)
ggsurvplot(sfit1,pval =TRUE, data = meta)

# 以基因 TNMD 为例
dev.off()
g = 'AC138207.1'
meta$gene = ifelse(as.numeric(exprTumor[g,])> median(as.numeric(exprTumor[g,])),'high','low')
sfit1=survfit(Surv(time, event)~gene, data=meta)
ggsurvplot(sfit1,pval =TRUE, data = meta)
png(filename = "survival_probability.png",width =800, height =800, res=1200)
dev.off()

# 多个基因的生存曲线
gs=rownames(exprTumor)[1:4]
splots <- lapply(gs, function(g){
  meta$gene=ifelse(as.numeric(exprTumor[g,])> median(as.numeric(exprTumor[g,])),'high','low')
  sfit1=survfit(Surv(time, event)~gene, data=meta)
  ggsurvplot(sfit1,pval =TRUE, data = meta)
}) 
arrange_ggsurvplots(splots, print = TRUE,  
                    ncol = 2, nrow = 2, risk.table.height = 0.4)

# p value
log_rank_p <- lapply(gs, function(gene){
  meta$group=ifelse(as.numeric(exprTumor[gene,])> median(as.numeric(exprTumor[gene,])),'high','low')
  data.survdiff=survdiff(Surv(time, event)~group,data=meta)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  a <- data.frame(HR)
  b <- data.frame(p.val)
  l <- c(a, b)
  return(l)
})

HR <- unlist(lapply(log_rank_p, '[[', 1))
PVALUE <- unlist(lapply(log_rank_p, '[[', 2))

# 批量生存分析
# 删除表达量全是0的数据，否则logrank会报错
dim(exprTumor)
drop <- which(apply(exprTumor, 1, max) == 0)
exprTumor <- exprTumor[-drop,]
dim(exprTumor)
# logrank批量生存分析
# Cox 回归的重要统计指标： 风险比（ hazard ratio）
# 当 HR>1 时，说明研究对象是一个危险因素。
# 当 HR<1 时，说明研究对象是一个保护因素。
# 当 HR=1 时，说明研究对象对生存时间不起作用。
log_rank_p <- lapply(rownames(exprTumor), function(gene){
  meta$group=ifelse(as.numeric(exprTumor[gene,])> median(as.numeric(exprTumor[gene,])),'high','low')
  data.survdiff=survdiff(Surv(time, event)~group,data=meta)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  a <- data.frame(HR)
  b <- data.frame(p.val)
  l <- c(a, b)
  return(l)
  })

HR <- unlist(lapply(log_rank_p, '[[', 1)) # 提取每个列表元素的第1个元素
PVALUE <- unlist(lapply(log_rank_p, '[[', 2)) # 提取每个列表元素的第2个元素

exprTumor_KM_pvalue <- as.data.frame(rownames(exprTumor))
exprTumor_KM_pvalue$HR <- HR
exprTumor_KM_pvalue$PVALUE <- PVALUE
colnames(exprTumor_KM_pvalue) <- c('Gene', 'HR', 'pvalue')
dim(exprTumor_KM_pvalue)
object.size(exprTumor_KM_pvalue)

# 导出数据
write.csv(exprTumor_KM_pvalue,file='KM_HR_PVALUE.csv',row.names = F,quote=F)

