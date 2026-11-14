#!/bin/bash
# Count reads in every FASTQ file in data/
for f in data/*.fastq; do
  echo "$f: $(( $(wc -l < $f) / 4 )) reads"
done
