#!/bin/bash
# Download the first 2M read pairs for each airway sample
for s in $(tail -n +2 data/samples.tsv | cut -f1); do
  if [ -f data/${s}_1.fastq.gz ]; then echo "$s already done"; continue; fi
  echo "Downloading $s..."
  fastq-dump -X 2000000 --split-files --gzip -O data $s
done
