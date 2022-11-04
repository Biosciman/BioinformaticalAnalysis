#!/usr/bin/perl -w
use strict;
# Selected the c_AMPs that appears in the proteome
# usage perl run.pl key_values_01 proteom_less_50.fa key_values_01.out

my $in = $ARGV[0];
my $b;
my %h;
open I,"<$in";
while(defined($b=<I>)){
    $b =~ s/\n//igm;
    $b =~ s/(?<se>.*:)(?<k>.*)//igm;
    my $b1 = $+{se};
    my $pp = $+{k};
    $b1 =~ s/\s//igm;
    $b1 =~ s/://igm;
    $pp =~ s/\s//igm;
    $h{$b1} = $pp;
}
close I;

my $in1 = $ARGV[1];
my $a;
open II,"<$in1";
while(defined($a=<II>)){
    $a =~ s/\n//igm;
    $a =~ s/\s//igm;
    $a =~ s/\*//igm;
    if($a =~ />/){}else{
        if($h{$a}){
            print"$a:$h{$a}\n";
        }
    }
}
close II;
