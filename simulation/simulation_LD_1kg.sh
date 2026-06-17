#!/bin/bash

Dir=~/Documents/GWAS/data/1000GENOME/phase3

for pid in EUR AFR EAS
do
  plink --file $Dir/$pid/chr1 --extract data/SNP_1kg.txt --recode \
  --out data/1kg_${pid}
  
  ls data/LD |while read block
  do
    if [ ! -e data/LD/$block/LD_1kg_${pid}.ld.gz ]
    then
      plink --file data/1kg_${pid} --extract data/LD/$block/SNP_1kg.txt \
      -r square gz --out data/LD/$block/LD_1kg_${pid}
    fi
  done
done
