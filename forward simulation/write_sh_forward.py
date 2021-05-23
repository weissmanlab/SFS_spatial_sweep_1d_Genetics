# -*- coding: utf-8 -*-

prog_str = '''#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 1-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=10000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH  --account=arielamir_lab

module load gcc/7.1.0-fasrc01 gsl/2.5-fasrc01
gcc plain_sweep_original_su.c -o plain_sweep -lgsl -lm
./plain_sweep.exe {} {} {:.6f} {:.6f} {} 0
'''
#./plain_sweep.exe L N s m tfinal ranseed
import io

L = 5000
N = 200
s = 0.05
m = 0.25
tfinal = 100000

filename = 'L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.sh'.format(
        L, N, s, m, tfinal)
print('Saving {}'.format(filename))
with io.open(filename, 'w', newline='\n') as f:
    f.write(prog_str.format(
            L, N, s, m, tfinal,
            L, N, s, m, tfinal, 
            L, N, s, m, tfinal))