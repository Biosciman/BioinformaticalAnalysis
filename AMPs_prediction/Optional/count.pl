#!/usr/bin/perl -w
# 选择小于50个氨基酸的序列
# perl count.pl all_seq.fa > less_50_data.fa
use strict;
my $in = $ARGV[0];
my $a;
my $temp = "";
open I,"<$in";
while(defined($a=<I>)){
  $a =~ s/\n//gm;
  $a =~ s/\s/_/gm; #\s表示空格
  if($a =~ />/){
    $temp = $a;
  }else{
      # @定义数组，并用split函数将每个氨基酸分割进数组seq
    my @len = split //, $a;
    my $le = @len;
      # &&如果两个操作数都为true，则条件为true
    if(4 < $le && $le < 51){
      if($temp){
        print"$temp\n$a\n";
      }else{
        print "$a\n";
      }
    }
  }
}
close I; 
