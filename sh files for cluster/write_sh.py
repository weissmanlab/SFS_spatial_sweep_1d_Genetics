# -*- coding: utf-8 -*-

prog_str = '''#!/bin/bash
#SBATCH -n 20                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 7-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=10000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_sample_all_{}.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_sample_all_{}.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH  --account=arielamir_lab

module load python/3.7.7-fasrc01
python SFS_from_frequency_parallel_recombination_vectorized_sample_all_multiple_forward.py {} {} {:.6f} {:.6f} {:.6f} {} {} {} {} {} # last digit is for copy number of the foward time sim
'''

import io

L = 5000
N = 200
s = 0.05
m = 0.25
r = 0
tfinal = 1000000
Un = 1
nsample = 100000
nSFS = 1000

for n in range(100):
    filename = 'L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_Un={}_nsample={}_sample_all_{}.sh'.format(
            L, N, s, m, r, tfinal, Un, nsample, n)
    print('Saving {}'.format(filename))
    with io.open(filename, 'w', newline='\n') as f:
        f.write(prog_str.format(
                L, N, s, m, r, tfinal, nsample, n,
                L, N, s, m, r, tfinal, nsample, n, 
                L, N, s, m, r, tfinal, Un, nsample, nSFS, n))