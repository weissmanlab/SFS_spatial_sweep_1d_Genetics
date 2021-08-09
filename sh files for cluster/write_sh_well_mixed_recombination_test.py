# -*- coding: utf-8 -*-

prog_str = '''#!/bin/bash
#SBATCH -n 25                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 7-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=500           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o well_mixed_N={}_Tfix={}_s={:.2f}_r={:.2e}.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e well_mixed_N={}_Tfix={}_s={:.2f}_r={:.2e}.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH  --account=desai_lab

module load python/3.7.7-fasrc01
python well_mixed_simulation_add_recombination2.py {} {:.2f} {} {} {} {:.9f}
'''

import io
import numpy as np
N = 10 ** 6
s = 0.05
r = 0.001
nsample = 10 ** 6
Nsim = 10000
Tfix = 0
rlist = 1 / 10 ** np.arange(3, 9)
for r in rlist:
    filename = 'expected_SFS_well_mixed_N={}_Tfix={}_s={:.2f}_r={:.2e}.sh'.format(
            N, Tfix, s, r)
    print('Saving {}'.format(filename))
    with io.open(filename, 'w', newline='\n') as f:
        f.write(prog_str.format(
                N, Tfix, s, r,
                N, Tfix, s, r, 
                N, s, Nsim, nsample, Tfix, r))