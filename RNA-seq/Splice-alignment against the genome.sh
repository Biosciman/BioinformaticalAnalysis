#!/bin/bash
tophat -p 3 --bowtie1 --no-novel-juncs --read-mismatches 2 -G corrected.gtf --transcriptome-index transcriptome_data_Aaegl5 --output-dir PBS_HD_1 index/AaegL5 PBS_HD_1.fq
tophat -p 3 --bowtie1 --no-novel-juncs --read-mismatches 2 -G corrected.gtf --transcriptome-index DENV2_transcriptome_Aaegl5 --output-dir DENV2_HD_1 index/AaegL5 DENV2_HD_1.fq
