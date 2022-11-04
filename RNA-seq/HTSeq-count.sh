#!/bin/bash
samtools view -h -o accepted_hits.sam accepted_hits.bam
htseq-count --stranded no -i gene_id accepted_hits.sam corrected.gtf > accepted_hits.gene.counts