#!/bin/bash
makeblastdb -in Aeg_up_protein.fa -dbtype prot -out Aeg_up_protein.db.fa
makeblastdb -in Aeg_down_protein.fa -dbtype prot -out Aeg_down_protein.db.fa
