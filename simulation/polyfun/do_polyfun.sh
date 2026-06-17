#!/bin/bash

for t in UKBB 1kg
do
  for s in 1 2 3
  do
    mkdir -p output/$t/setting_${s}

    for suffix in N1-200000
    do
      python /home/r9user9/Documents/GWAS/packages/polyfun/polyfun.py \
      --compute-h2-L2 \
      --no-partitions \
      --output-prefix output/$t/setting_${s}/$suffix \
      --sumstats input/sumstat/setting_${s}/sumstat.${suffix}.parquet \
      --ref-ld-chr input/annotations/$t/EUR/annotations. \
      --w-ld-chr input/annotations/$t/EUR/weights.
    done

    for suffix in N2-20000 N2-200000
    do
      python /home/r9user9/Documents/GWAS/packages/polyfun/polyfun.py \
      --compute-h2-L2 \
      --no-partitions \
      --output-prefix output/$t/setting_${s}/$suffix \
      --sumstats input/sumstat/setting_${s}/sumstat.${suffix}.parquet \
      --ref-ld-chr input/annotations/$t/AFR/annotations. \
      --w-ld-chr input/annotations/$t/AFR/weights.
    done

    for suffix in N3-20000 N3-200000
    do
      python /home/r9user9/Documents/GWAS/packages/polyfun/polyfun.py \
      --compute-h2-L2 \
      --no-partitions \
      --output-prefix output/$t/setting_${s}/$suffix \
      --sumstats input/sumstat/setting_${s}/sumstat.${suffix}.parquet \
      --ref-ld-chr input/annotations/$t/EAS/annotations. \
      --w-ld-chr input/annotations/$t/EAS/weights.
    done
  done
done