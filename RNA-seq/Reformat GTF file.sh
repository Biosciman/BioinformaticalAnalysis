#!/bin/bash
perl -ne 's/(gene_name\s\"[\d|\w|-]+)(;)([\d|\w|-]+\")/$1.$3/; print' Aedes_aegypti_lvpagwg.AaegL5.52.gtf > corrected.gtf
