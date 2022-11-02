#!/usr/bin/perl -w
use strict;
my $cmd = "getorf -sequence GCF_000005845.2_ASM584v2_genomic.fna -find 0 -table 11 -minsize 15 -maxsize 150 -outseq ASM584v2.orf.fa";
system("$cmd");
