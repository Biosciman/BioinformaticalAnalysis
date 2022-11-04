#!/bin/bash
cuffdiff -o diff_out -p 2 -b GCF_002204515.2_AaegL5.0_genomic.fna -u corrected.gtf -L PBS,DENV2 PBS_HD_1/accepted_hits.bam DENV2_HD_1/accepted_hits.bam