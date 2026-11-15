#!/bin/bash
# Phase 4: quality control
fastqc -t 4 --outdir qc data/*.fastq.gz
# then in the multiqc env: multiqc qc --outdir qc
