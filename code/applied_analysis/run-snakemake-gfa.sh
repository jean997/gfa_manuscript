#!/bin/bash

mkdir -p log
snakemake \
   -s Snakefile_gfa \
   --keep-going \
   --jobs 96 \
   --max-jobs-per-second 5 \
   --latency-wait 30 \
   --notemp \
   --cluster-config cluster.yaml  \
   --cluster "sbatch \
              --output={cluster.log}_%j.out \
              --error={cluster.log}_%j.err \
              --job-name={cluster.name} \
              --time={cluster.time}  \
              --cpus-per-task={cluster.cpus}  \
              --mem={cluster.mem}"
