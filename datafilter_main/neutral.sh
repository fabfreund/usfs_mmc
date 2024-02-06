#!/bin/bash
##transtorm a SNP file with the allele in 1st and 3rd position into a neutral SNP file
#Select for A-T or G-C alleles
## to use for Aptenodytes (after running the first filter script): ./neutral.sh SNPfile_Aptenodytes.txt
file=$1
awk '{if (($1=="A" && $3=="T") || ($1=="T" && $3=="A") || ($1=="C" && $3=="G") || ($1=="G" && $3=="C")) print}' ${file} > neutral_${file}
