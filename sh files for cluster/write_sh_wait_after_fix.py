# -*- coding: utf-8 -*-

prog_str = '''#!/bin/bash
#SBATCH -n 20                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 7-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=8000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tfix={}_sample_all_{}.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_nsample={}_tifx={}_sample_all_{}.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH  --account=arielamir_lab

module load python/3.7.7-fasrc01
python expected_SFS_parallel_recombination_vectorized_sample_all_multiple_forward_wait_after_fixation.py {} {} {:.6f} {:.6f} {:.6f} {} {} {} {} {} {}   # last digit is for copy number of the foward time sim
'''

import io

L = 500
N = 20000
s = 0.05
m = 0.25
r = 0.2
tfinal = 1000000
Un = 1
nsample = 100000
nSFS = 1000
tfix = 100000

for n in range(100):
    filename = 'expected_L={}_N={}_s={:.6f}_m={:.6f}_r={:.6f}_tfinal={}_Un={}_nsample={}_tfix={}_sample_all_{}.sh'.format(
            L, N, s, m, r, tfinal, Un, nsample, tfix, n)
    print('Saving {}'.format(filename))
    with io.open(filename, 'w', newline='\n') as f:
        f.write(prog_str.format(
                L, N, s, m, r, tfinal, nsample, tfix, n,
                L, N, s, m, r, tfinal, nsample, tfix, n, 
                L, N, s, m, r, tfinal, Un, nsample, nSFS, tfix, n))