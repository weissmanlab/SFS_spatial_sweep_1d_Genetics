# -*- coding: utf-8 -*-

prog_str = '''#!/bin/bash
#SBATCH -n 25                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 7-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=5000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o L={}_rho={}_s={:.2e}_m={:.2e}_r={:.2e}_tfinal={}_nsample={}_tfix={}_sample_all_{}.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e L={}_rho={}_s={:.2e}_m={:.2e}_r={:.2e}_tfinal={}_nsample={}_tifx={}_sample_all_{}.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH  --account=desai_lab

module load python/3.7.7-fasrc01
python 1d_spatial_adding_recombination_test2.py {} {} {:.8f} {:.8f} {:.8f} {} {} {} {} {}   # last digit is for copy number of the foward time sim
'''

import io

L = 500
rho = 20000
s = 0.05
m = 0.25
r = 0.001
tfinal = 1000000
nsample = 100000
nSFS = 1000
tfix = 0

for n in range(100):
    filename = 'expected_L={}_rho={}_s={:.2e}_m={:.2e}_r={:.2e}_tfinal={}_nsample={}_tfix={}_sample_all_{}.sh'.format(
            L, rho, s, m, r, tfinal, nsample, tfix, n)
    print('Saving {}'.format(filename))
    with io.open(filename, 'w', newline='\n') as f:
        f.write(prog_str.format(
                L, rho, s, m, r, tfinal, nsample, tfix, n,
                L, rho, s, m, r, tfinal, nsample, tfix, n, 
                L, rho, s, m, r, tfinal, nsample, nSFS, tfix, n))