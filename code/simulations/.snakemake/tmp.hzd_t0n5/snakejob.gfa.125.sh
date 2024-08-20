#!/bin/sh
# properties = {"type": "single", "rule": "gfa", "local": false, "input": ["../../results/simulation_results/simulate/simulate.A.3.0.RDS", "../../results/simulation_results/R_ests/R_ests.A.3.0.RDS"], "output": ["../../results/simulation_results/gfa/gfa_pt0.05_pval_1.A.3.0.RDS"], "wildcards": {"Rstring": "pt0.05", "ldtype": "pval", "pthresh": "1", "scenario": "A", "ndense": "3", "rep": "0"}, "params": {}, "log": [], "threads": 1, "resources": {"mem_mb": 1895, "mem_mib": 1808, "disk_mb": 1895, "disk_mib": 1808, "tmpdir": "<TBD>"}, "jobid": 125, "cluster": {"mem": "5G", "cpus": "1", "name": "gfa-Rstring=pt0.05,ldtype=pval,ndense=3,pthresh=1,rep=0,scenario=A", "log": "log/snake-gfa-Rstring=pt0.05,ldtype=pval,ndense=3,pthresh=1,rep=0,scenario=A", "time": "2:00:00"}}
cd /gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_code_ocean/git_repository/gfa_manuscript/code/simulations && /home/jvmorr/miniconda3/envs/snakemake/bin/python3.11 -m snakemake --snakefile '/gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_code_ocean/git_repository/gfa_manuscript/code/simulations/Snakefile_sims' --target-jobs 'gfa:Rstring=pt0.05,ldtype=pval,pthresh=1,scenario=A,ndense=3,rep=0' --allowed-rules 'gfa' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=1895' 'mem_mib=1808' 'disk_mb=1895' 'disk_mib=1808' --wait-for-files '/gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_code_ocean/git_repository/gfa_manuscript/code/simulations/.snakemake/tmp.hzd_t0n5' '../../results/simulation_results/simulate/simulate.A.3.0.RDS' '../../results/simulation_results/R_ests/R_ests.A.3.0.RDS' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'software-env' 'input' 'params' 'mtime' 'code' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --latency-wait 30 --scheduler 'ilp' --scheduler-solver-path '/home/jvmorr/miniconda3/envs/snakemake/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && touch '/gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_code_ocean/git_repository/gfa_manuscript/code/simulations/.snakemake/tmp.hzd_t0n5/125.jobfinished' || (touch '/gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_code_ocean/git_repository/gfa_manuscript/code/simulations/.snakemake/tmp.hzd_t0n5/125.jobfailed'; exit 1)

