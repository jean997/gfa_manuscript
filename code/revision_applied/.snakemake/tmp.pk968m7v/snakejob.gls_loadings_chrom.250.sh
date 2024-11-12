#!/bin/sh
# properties = {"type": "single", "rule": "gls_loadings_chrom", "local": false, "input": ["../../data/gfa_intermediate_data/bc_zmat.16.RDS", "../../results/applied_analysis/gfa_results/bc_gfa_gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.final.RDS", "../../data/gfa_intermediate_data/bc_R_estimate.R_none.RDS"], "output": ["../../results/applied_analysis/gfa_results/bc_gls_loadings.gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.16.RDS"], "wildcards": {"prefix": "bc", "fs": "1", "ldstring": "r20.01_kb1000_pvalue", "Rstring": "none", "chrom": "16"}, "params": {}, "log": [], "threads": 1, "resources": {"mem_mb": 1000, "mem_mib": 954, "disk_mb": 1000, "disk_mib": 954, "tmpdir": "<TBD>"}, "jobid": 250, "cluster": {"mem": "5G", "cpus": "1", "name": "gls_loadings_chrom-Rstring=none,chrom=16,fs=1,ldstring=r20.01_kb1000_pvalue,prefix=bc", "log": "log/snake-gls_loadings_chrom-Rstring=none,chrom=16,fs=1,ldstring=r20.01_kb1000_pvalue,prefix=bc", "time": "2:00:00"}}
cd /gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_analysis/gfa_manuscript/code/revision_applied && /home/jvmorr/miniconda3/envs/snakemake/bin/python3.11 -m snakemake --snakefile '/gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_analysis/gfa_manuscript/code/revision_applied/Snakefile_gfa' --target-jobs 'gls_loadings_chrom:prefix=bc,fs=1,ldstring=r20.01_kb1000_pvalue,Rstring=none,chrom=16' --allowed-rules 'gls_loadings_chrom' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=1000' 'mem_mib=954' 'disk_mb=1000' 'disk_mib=954' --wait-for-files '/gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_analysis/gfa_manuscript/code/revision_applied/.snakemake/tmp.pk968m7v' '../../data/gfa_intermediate_data/bc_zmat.16.RDS' '../../results/applied_analysis/gfa_results/bc_gfa_gfaseed1.ldpruned_r20.01_kb1000_pvalue.R_none.final.RDS' '../../data/gfa_intermediate_data/bc_R_estimate.R_none.RDS' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'mtime' 'input' 'code' 'software-env' 'params' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --latency-wait 30 --scheduler 'ilp' --scheduler-solver-path '/home/jvmorr/miniconda3/envs/snakemake/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && touch '/gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_analysis/gfa_manuscript/code/revision_applied/.snakemake/tmp.pk968m7v/250.jobfinished' || (touch '/gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_analysis/gfa_manuscript/code/revision_applied/.snakemake/tmp.pk968m7v/250.jobfailed'; exit 1)

