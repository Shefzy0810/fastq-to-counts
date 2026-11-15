#!/bin/bash
# Download GENCODE v50 human transcriptome and build the Salmon index
mkdir -p ref
wget -nc -P ref https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.transcripts.fa.gz
salmon index -t ref/gencode.v50.transcripts.fa.gz -i ref/salmon_index --gencode -p 4
