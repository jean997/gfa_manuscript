#!/bin/bash

mkdir -p log
snakemake \
   -s Snakefile_gfa \
   --keep-going \
   --notemp \
   --jobs 96 \
   --max-jobs-per-second 5 \
   --latency-wait 30 \
   --executor slurm \
   --default-resources mem_mb=5000 runtime=120 \
   --cluster-submit-arg "--account=jvmorr0" \
   --job-name "gfa-{rule}-{wildcards}-{jobid}" \
   --jobscript "log/jobscript.{rule}-{wildcards}-{jobid}.sh" \
   --cluster-output "log/log.{rule}-{wildcards}-{jobid}.out"

