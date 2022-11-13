#!/bin/bash
blastp -query Alb_up_protein.fa -db Aeg_up_protein.db.fa -out Alb_up_protein_blast.tsv -outfmt 6 -culling_limit 1
blastp -query Alb_down_protein.fa -db Aeg_down_protein.db.fa -out Alb_down_protein_blast.tsv -outfmt 6 -culling_limit 1
