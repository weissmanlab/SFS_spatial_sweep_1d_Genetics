#!/bin/bash
#SBATCH -n 20                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 1-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=5000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o L=500_rho=20000_s=5.00e-02_m=2.50e-01_r=1.00e-07_nsample=1000_nback=1000_tfix=0_0.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e L=500_rho=20000_s=5.00e-02_m=2.50e-01_r=1.00e-07_nsample=1000_nback=1000_tifx=0_0.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH  --account=desai_lab

module load python/3.7.7-fasrc01
python backwards_combined.py 500 20000 0.05000000 0.25000000 0 1000 1000 0 0 1 0.00000010   # last digit is for copy number of the foward time sim
