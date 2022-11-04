#!/usr/bin/perl -w
use strict;
=pod
Transformation of c_AMPs sequences into fragments of a certain length, here we setting the length were contain from full length to half length

将AMP序列进行kmer分割，这边进行50%序列长度至100%序列长度的kmer

usage perl kmer.pl Sequence.fa > mers.fa
=cut

my $in = $ARGV[0]; # Sequence.fa from "Result integration" section, and remove 2 head line.
my $b;
my %h;
open I,"<$in";
while(defined($b=<I>)){
    $b =~ s/\n//igm;
    $b =~ s/\s//igm;
    $b =~ s/\*//igm;
    if($b =~ />/){}else{
        my @pr = &kmer($b); # &为调用子函数
        my $mers = join(',',@pr);
#        print "$b:$mers\n";
    }
}
close I;

# sub定义一个子函数，kmer是名称
sub kmer{
    (my $sequ) = @_; # @_为子函数的私有变量
    my @a = split //, $sequ;
    my $le = @a;
    my $len = int($le/2)-1;
    my $ii = 0;
    my @se;
    while($ii < $le){
       my $leng_set = $len;
       while($leng_set <= $le){
           if($ii+$leng_set < $le){
               # join()函数使用VAR的值将LIST的元素组合为单个字符串
               my $temp = join('',@a[$ii..$leng_set+$ii]); # .. 为范围运算符,先计算$leng_set+$ii
               # push()函数用于将值列表推送到数组的末尾
               push @se, $temp;
           }
           $leng_set++; # ++ 自增运算符 值加1
       }
       $ii++;
    }
return @se;
}
