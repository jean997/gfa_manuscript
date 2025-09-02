#!/bin/bash

mkdir -p log
snakemake \
   -s Snakefile_gfa \
   --keep-going \
   --notemp \
   --jobs 96 \
   --max-jobs-per-second 5 \
   --latency-wait 30 \
   --default-resources mem_mb=5000 runtime=120 account=jvmorr0 \
   --executor slurm

