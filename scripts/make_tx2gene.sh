#!/bin/bash
# Build transcript -> gene table from GENCODE FASTA headers
# Header format: >ENST|ENSG|OTTHUMG|OTTHUMT|tx_name|gene_name|length|biotype|
zcat ref/gencode.v50.transcripts.fa.gz | grep '^>' | sed 's/^>//' \
  | awk -F'|' 'BEGIN{OFS="\t"; print "tx","gene_id","gene_name"} {print $1,$2,$6}' \
  > ref/tx2gene.tsv
