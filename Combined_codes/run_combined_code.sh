#!/bin/bash
#
#SBATCH -p desai #serial_requeue # Partition to submit to (comma separated)
#SBATCH -n 20 # Number of cores (number of backward simulations you want to run for one forward simulation, parameter that goes in pool)
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-90:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 50000 # Memory in MB

#SBATCH -J Forward_4 # Job name

#SBATCH -o out_file4.txt
#SBATCH -e err_file4.txt

##Read this from parameters file and change N_forw by hand
#Nforw = 1
#L = 50
#N = 50
#s = 0.05
#m = 0.25
#l0 = 25
#nbase = 10
#n_SFS = 10
#T_after_fix = 3


#gcc -std=c99 Sweep_in-2D_version-m.c -o plain_sweep -lgsl -lgslcblas -lm -g
##Usage is L N s m T_final l0 Nforw 
#./plain_sweep 500 500 0.05 0.25 10000 250 1
##Usage is L N s m l0 n_base N_SFS T_after_fix Nforw dimension
python backwards_combined.py 500 500 0.05 0.250 250 500000 10 7300 4 2
#python Plot_SFS.py




#####READ ME
##### Each submission of this script runs backward simulation code one time (each simulation averages over n backwards).
##### What you need to change is the number at the end of Job name, name of output and error file, and the second last file parameter in the input (Nforw). They must be the same 
##### This is currently set to 4 and I submit this script 4 times, setting it to 1, 2, 3, and 4 separately.  
##### This ensures that the correct file names are written for 2-D. The 1-D file name string formatting for the input forward simulation has not been done right now. 
##### This also ensures that the code writes out the correct output files with corect names. 