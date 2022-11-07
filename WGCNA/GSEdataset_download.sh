#!/bin/bash
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48213
# 下载，-c断点续传
wget -c ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48213/suppl/GSE48213_RAW.tar
# 解压 -x从备份文件中还原文件（extract）,-f指定备份文件
tar -xf GSE48213_RAW.tar
# gzip -d 对压缩文件进行解压缩
gzip -d *.gz
