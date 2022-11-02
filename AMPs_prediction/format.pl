#!/usr/bin/perl -w
use strict;

=pod
Usage: perl format.pl AMP.te.fa P > amp.txt
Here, P means positive true lable, N means negative false lable, and
none means do not need using this function, using for predicted
file re-format shuold silence this function
=cut

my $in = $ARGV[0];
# %定义hash，即键值对
my %aacode =
   (
   A => "1", C => "2", D => "3", E => "4",
   F => "5", G => "6", H => "7", I => "8",
   K => "9", L => "10", M => "11", N => "12",
   P => "13", Q => "14", R => "15", S => "16",
   T => "17", V => "18", W => "19", Y => "20",
   );
my ($a,$b);
open I,"<$in";
while(defined($a=<I>)){
    #chomp函数，负责删除标量型变量中的最后一个字符，或者数组中每个字的最后一个字符，并保证只有该行末字符是换行符时才进行删除操作。它会返回删除后的字符数目。
   chomp($a);
    #uc函数将其转换为大写后返回传递给它的字符串。
   $a = uc $a;
    # 配平￥a中“>”符号,"\"为转义
   if($a =~ /\>/){
       $b = "0";
       }else{
       if($a =~ /B/){}
       elsif($a =~ /J/){}
       elsif($a =~ /O/){}
       elsif($a =~ /U/){}
       elsif($a =~ /X/){}
       elsif($a =~ /Z/){}
           # 排除BJOUXZ这些非氨基酸字母
       else{
           # @定义数组，并用split函数将每个氨基酸分割进数组seq
           my @seq = split //,$a;
           my $len = @seq; #获得数组元素个数，即氨基酸序列长度
           my $i = 0;
           while($i < $len){ #perl索引从0开始
               $b = $b.",".$aacode{$seq[$i]}; #将氨基酸替换为数字，并以“,”连接
               $i = $i + 1;
           }
           my @bb = split /,/, $b; #将数字化的氨基酸分割进数组bb
           my $lb = @bb; #获得数组元素个数，即氨基酸序列长度
           my $cha = 300 - $lb; #计算300与氨基酸个数的差值
           my $j = 0;
           while($j < $cha){
               $b = "0".",".$b; # 填充至300个氨基酸
               $j = $j + 1;
           }
           #eq为字符串比较运算符，为等于的意思；|为位算符，针对二进制数据的位运算，表示两个位都为0时，结果才为0
           if($ARGV[1] eq "P" | $ARGV[1] eq "p"){
               $b = $b.",1";
               print "$b\n";
           }elsif($ARGV[1] eq "N" | $ARGV[1] eq "n"){
               $b = $b.",0";
               print "$b\n";
           }elsif($ARGV[1] eq "none"){
               print "$b\n";
           }
       }
   }
}
close I;
