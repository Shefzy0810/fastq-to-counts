#!/bin/bash
# Quantify each sample against the GENCODE v50 Salmon index
for s in $(tail -n +2 data/samples.tsv | cut -f1); do
  if [ -f quant/$s/quant.sf ]; then echo "$s already done"; continue; fi
  echo "Quantifying $s..."
  salmon quant -i ref/salmon_index -l A \
    -1 data/${s}_1.fastq.gz -2 data/${s}_2.fastq.gz \
    -p 4 -o quant/$s
done
