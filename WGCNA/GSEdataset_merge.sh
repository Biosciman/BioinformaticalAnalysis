#!/bin/bash
# 在GSE48213_RAW目录里面生成tmp.txt文件
# AWK 是一种处理文本文件的语言
# $0 输出文本中的第0项
# |运算符，把左边的命令的执行结果作为输入传递给右边命令
# grep -v 显示不包含匹配文本的所有行。
awk '{print FILENAME"\t"$0}' GSM*.txt |grep -v EnsEMBL_Gene_ID >tmp.txt
