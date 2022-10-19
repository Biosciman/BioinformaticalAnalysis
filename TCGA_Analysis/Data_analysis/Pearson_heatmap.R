library("ggplot2")
library("dplyr")
library("reshape")
library("psych")
library("aplot")

filename = 'Final_gene_expression.csv'
table1 <- read.csv(filename,header=T,row.names = 1)
table2 <- table1

pp <- corr.test(table1,table2,method="pearson",adjust = "fdr")

#抽提出相关性系数与P值
cor <- pp$r 
pvalue <- pp$p
#根据P值用循环将P值做成*
display <- pvalue
l1 <- nrow(display);l2 <- ncol(display)
for(i in 1:l1){
  for(k in 1:l2){
    a <- as.numeric(display[i,k])
    if(a <= 0.001){
      a <- "***"
    }
    if( 0.001 < a && a <= 0.01){
      a <- "**"
    }
    if(0.01 < a && a < 0.05){
      a <- "*"
    }
    if(a >= 0.05){
      a <- ""
    }
    display[i,k] <- a
  }
}
#对数据进行格式转换适用于ggplot2绘图
heatmap <- melt(cor)
heatmap[,4] <- melt(pvalue)[,3]
heatmap[,5] <- melt(display)[,3]
names(heatmap) <- c("sample","gene","cor","pvalue","display")
#导出数据
write.table (heatmap,file ="heatmap.xls", sep ="\t", row.names = F)  

#ggplot2绘制热图
pp <- ggplot(heatmap,aes(gene,sample,fill=cor)) + 
  geom_tile()+
  theme_minimal()+scale_fill_viridis_c()+
  geom_text(aes(label=display),size=5,color="white")+
  scale_y_discrete(position="right")+xlab(NULL) + ylab(NULL)+
  theme(axis.text.x = element_text(angle = 90,hjust=0.8,vjust=0.6))+
  theme(axis.text.x=element_text(family = "Times",
                                 face = "plain",colour = "black",size=10))+
  theme(axis.text.y=element_text(family = "Times",
                                 face = "plain",colour = "black",size=10))+
  theme(legend.text=element_text(face="plain",family = "Times",
                                 colour = "black",size = 10))+
  labs(fill = "expr")
pp

