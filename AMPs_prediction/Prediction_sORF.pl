#!/usr/bin/perl -w
use strict; # 使用strict后，声明变量需要加上my
my $in = "fasta_sequence_names.txt"; # my定义变量in，in为序列名称文件
my ($a, $cmd); # my定义变量a和cmd
open I, "<$in"; # I为文件名，"<$in"中<表示仅读，前提文件in必须已经存在
while(defined($a=<I>)){
    # defined()返回bool,从文件中读取一行数据存储到变量$a中并把Perl文件指针向后移动一行。
    $a =~ s/\n//igm;
    # 正则表达式，~表示匹配，s为替换模式，i为忽略模式中的大小写，g为全局匹配，m为多行模式了；此处删除换行\n
    my $out_name = $a.".orf.fa"; # .为字符串的连接
    $cmd = "getorf -sequence /fasta_sequence_path/".$a." -find 0 -table 11 -minsize 15 -maxsize 150 -outseq /output_sequnces_path/".$out_name; # EMBOSS软件命令
    print "$cmd\n";
    system("$cmd");
    # 调用系统函数，即将$cmd中字符串命令在计算机终端中运行,双引号里面$变量为真变量
}
close I;
