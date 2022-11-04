#!/usr/bin/perl -w
use strict;
# usage perl selec_seq_frmt.pl selected_c_AMPs_out.txt > selected_c_AMPs.fa

my $a;
my %h;
my $i = 1;
open I,"<$ARGV[0]";
while(defined($a=<I>)){
    chomp($a);
    my @aa = split(/:/, $a);
    my $pp = $aa[1];
    $pp =~ s/,/\n/igm;
    $h{$pp} = ">c_APM".$i;
    $i++;
}
close I;

foreach my $k (keys %h){
    print "$h{$k}\n$k\n";
}
