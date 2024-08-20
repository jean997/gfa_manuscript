#!/bin/sh
# properties = {"type": "single", "rule": "gfa_R", "local": false, "input": ["../../results/simulation_results/simulate/simulate.A.3.0.RDS", "../../results/simulation_results/R_ests/R_ests.A.3.0.RDS"], "output": ["../../results/simulation_results/gfa/gfa_none_random_1.A.3.0.RDS"], "wildcards": {"Rstring": "none", "ldtype": "random", "pthresh": "1", "scenario": "A", "ndense": "3", "rep": "0"}, "params": {}, "log": [], "threads": 1, "resources": {"mem_mb": 1895, "mem_mib": 1808, "disk_mb": 1895, "disk_mib": 1808, "tmpdir": "<TBD>"}, "jobid": 243, "cluster": {"mem": "5G", "cpus": "1", "name": "gfa_R-Rstring=none,ldtype=random,ndense=3,pthresh=1,rep=0,scenario=A", "log": "log/snake-gfa_R-Rstring=none,ldtype=random,ndense=3,pthresh=1,rep=0,scenario=A", "time": "2:00:00"}}
cd /gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_code_ocean/git_repository/gfa_manuscript/code/simulations && /home/jvmorr/miniconda3/envs/snakemake/bin/python3.11 -m snakemake --snakefile '/gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_code_ocean/git_repository/gfa_manuscript/code/simulations/Snakefile_sims' --target-jobs 'gfa_R:Rstring=none,ldtype=random,pthresh=1,scenario=A,ndense=3,rep=0' --allowed-rules 'gfa_R' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=1895' 'mem_mib=1808' 'disk_mb=1895' 'disk_mib=1808' --wait-for-files '/gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_code_ocean/git_repository/gfa_manuscript/code/simulations/.snakemake/tmp.ko241ogt' '../../results/simulation_results/simulate/simulate.A.3.0.RDS' '../../results/simulation_results/R_ests/R_ests.A.3.0.RDS' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'input' 'params' 'software-env' 'code' 'mtime' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --latency-wait 30 --scheduler 'ilp' --scheduler-solver-path '/home/jvmorr/miniconda3/envs/snakemake/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && touch '/gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_code_ocean/git_repository/gfa_manuscript/code/simulations/.snakemake/tmp.ko241ogt/243.jobfinished' || (touch '/gpfs/accounts/jvmorr_root/jvmorr0/jvmorr/GFA_code_ocean/git_repository/gfa_manuscript/code/simulations/.snakemake/tmp.ko241ogt/243.jobfailed'; exit 1)

