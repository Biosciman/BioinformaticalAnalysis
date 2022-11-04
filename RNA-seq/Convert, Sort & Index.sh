#!/bin/bash
samtools view -b -S RNAseq.sam -o RNAseq.bam
samtools sort RNAseq.bam -o RNAseq.sorted.bam
samtools index RNAseq.bam
samtools tview RNAseq.bam genome