#!/usr/bin/perl -w
use strict;

=pod
usage perl k_mer_dict.pl mers.fa > key_values_01

Set mers.fa to dictionary data format and every single mers from c_AMPs were setted as dictionary keys, the original c_AMPs were dictionary values.

输出kmers：全长序列的 键值对
=cut

my $in = $ARGV[0];
my $b;
my %h;
my $b1;
open I,"<$in";
while(defined($b=<I>)){
    $b =~ s/\n//igm;
    $b =~ s/(?<se>.*:)(?<k>.*)//igm;
    $b1 = $+{se}; # 定义se，即AMP全长序列
    my $pp = $+{k}; # 定义k，即AMP的kmers
    my @pr = split /,/, $pp;
    $b1 =~ s/://igm;
    dic1(@pr,$b1);
}
close I;

sub dic1{
    (my @ks) = @_; # 传参,@pr和$b1一起传入
    my $va = $b1; # va为AMP全长
    foreach my $kk (@ks){
        $h{$kk}{$va} = 1; #多维哈希，第一个键为$kk，第二个键为$va，值为1
    }
}


# keys函数返回hash所有键的列表
foreach my $ke (keys %h){ # 获得第一维键，即kmers
    my $h2 = $h{$ke}; # 将“AMP全长：1”此键值对赋值给h2
    my $k2 = join(',',keys %$h2); #获得第二维键，即AMP全长序列
    print "$ke: ", $k2,"\n";
}
