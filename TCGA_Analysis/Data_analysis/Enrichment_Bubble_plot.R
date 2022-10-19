library(GOplot)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

circ2 <- read.csv(file='GO_enrichment_up.csv',
                  header=T)


## 第一种方案
dev.off
ggplot(circ2,aes(x=1,y=-log10(adj_pval)))+
  geom_point(aes(size=count,color=category),alpha=0.6)+
  scale_size(range=c(1,12))+
  scale_color_brewer(palette = "Accent")+
  theme_bw()+
  theme(
    legend.position = c("none")
  )+
  geom_text_repel(
    data = circ2[-log10(circ2$adj_pval)>3,],
    aes(label = ID),
    size = 3,
    segment.color = "black", show.legend = FALSE )+
  facet_grid(.~category)+
  geom_text_repel(
    data = circ2[which(circ2$category=='GO:MF'),],
    aes(label = ID),
    size = 3,
    segment.color = "black", show.legend = FALSE )+
  facet_grid(.~category)


## 第二种方案
dev.off
ggplot(circ2,aes(x=random,y=-log10(adj_pval)))+
  geom_point(aes(size=count,color=category),alpha=0.6)+
  scale_size(range=c(1,36))+
  scale_color_brewer(palette = "Accent")+
  theme_bw()+
  theme(
    #legend.position = c("none")
  )+
  geom_text_repel(
    data = circ2[-log10(circ2$adj_pval)>3 & circ2$count>90,],
    aes(label = ID),
    size = 3,
    segment.color = "black", show.legend = FALSE )+
  xlab('GC-upregulated HUGs')+
  theme(panel.grid =element_blank())+ # 删除网格线
  theme(axis.text.x = element_blank())   ## 删去x刻度标签



circ1 <- read.csv(file='/Users/zhangliming/Documents/Bioinformatics/GSE129219/REAC_enrichment_up.csv',
                  header=T)
dev.off
ggplot(circ1,aes(x=random,y=-log10(adj_pval)))+
  geom_point(aes(size=count, color=1),alpha=0.6)+
  scale_size(range=c(1,36))+
  scale_color_gradient(low="red")+
  theme_bw()+
  theme(
    #legend.position = c('none')
  )+
  geom_text_repel(
    data = circ1[-log10(circ1$adj_pval)>1 & circ1$count>35,],
    aes(label = ID),
    size = 3,
    segment.color = "black", show.legend = FALSE )+
  xlab('GC-upregulated HUGs')+
  theme(panel.grid =element_blank())+ # 删除网格线
  theme(axis.text.x = element_blank())   ## 删去x刻度标签
  
  
