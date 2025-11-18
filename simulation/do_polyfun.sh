#!/bin/bash

for s in 1 2 3
do
  mkdir -p polyfun_output/setting_${s}

  for suffix in N1-200000
  do
    python /home/r9user9/Documents/GWAS/packages/polyfun/polyfun.py \
    --compute-h2-L2 \
    --no-partitions \
    --output-prefix polyfun_output/UKBB/setting_${s}/$suffix \
    --sumstats polyfun_input/sumstat/setting_${s}/sumstat.${suffix}.parquet \
    --ref-ld-chr polyfun_input/EUR/annotations. \
    --w-ld-chr polyfun_input/EUR/weights.
  done

  for suffix in N2-20000 N2-200000
  do
    python /home/r9user9/Documents/GWAS/packages/polyfun/polyfun.py \
    --compute-h2-L2 \
    --no-partitions \
    --output-prefix polyfun_output/UKBB/setting_${s}/$suffix \
    --sumstats polyfun_input/sumstat/setting_${s}/sumstat.${suffix}.parquet \
    --ref-ld-chr polyfun_input/EAS/annotations. \
    --w-ld-chr polyfun_input/EAS/weights.
  done
done


for s in 1 2 3
do
  mkdir -p polyfun_output/${t}_1kg/setting_${s}
  
  for suffix in N1-200000
  do
    python /home/r9user9/Documents/GWAS/packages/polyfun/polyfun.py \
    --compute-h2-L2 \
    --no-partitions \
    --output-prefix polyfun_output/1kg/setting_${s}/$suffix \
    --sumstats polyfun_input/sumstat/setting_${s}/sumstat.${suffix}.parquet \
    --ref-ld-chr polyfun_input/EUR_1kg/annotations. \
    --w-ld-chr polyfun_input/EUR_1kg/weights.
  done
  
  for suffix in N2-20000 N2-200000
  do
    python /home/r9user9/Documents/GWAS/packages/polyfun/polyfun.py \
    --compute-h2-L2 \
    --no-partitions \
    --output-prefix polyfun_output/1kg/setting_${s}/$suffix \
    --sumstats polyfun_input/sumstat/setting_${s}/sumstat.${suffix}.parquet \
    --ref-ld-chr polyfun_input/EAS_1kg/annotations. \
    --w-ld-chr polyfun_input/EAS_1kg/weights.
  done
done
