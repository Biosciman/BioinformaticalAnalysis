# 参考文章
# https://cloud.tencent.com/developer/article/1516749

# 安装R包
# BiocManager::install("WGCNA")
# install.packages("WGCNA",type = 'binary')

# 加载R包
library("WGCNA")

# 输入数据准备
# 主要是表达矩阵，如果是芯片数据，那么常规的归一化矩阵即可
# 如果是转录组数据，最好是RPKM/TPM值或者其它归一化好的表达量
# 然后就是临床信息或者其它表型，总之就是样本的属性
# 为了保证后续脚本的统一性，表达矩阵统一用datExpr标识，临床 信息统一用datTraits标识
filename1 = 'STAD_Gene_fpkm_unstrand.csv'
filename2 = 'STAD_clinicaldata_definition2.txt'