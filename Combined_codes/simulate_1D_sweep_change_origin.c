/*
 *  Simple simulation of a sweep in space + origin is distributed uniformly
 *
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>

// run by typing : gcc plain_sweep_change_sweep_origin.c -o plain_sweep -lgsl -lm
// ./plain_sweep.exe L N s m tfinal l0



// Global variables
const gsl_rng_type * T;
gsl_rng * R;

void next_gen(unsigned int n[], double s, double mig, int L, unsigned int N) {

    double x[L]; // frequencies
    double xmig[L]; // after migration
    int i;

    // Initialize frequencies:
    for (i = 0; i < L; i++){
        x[i] = (double)n[i] / N;
    }
    // Migration:
    for (i = 1; i < L - 1; i++){
        xmig[i] = x[i] + mig * (0.5 * (x[i-1] + x[i+1]) - x[i]);
    }
    xmig[0] = x[0] + mig * 0.5 * (x[1] - x[0]);
    xmig[L-1] = x[L-1] + mig * 0.5 * (x[L-2] - x[L-1]);

    // Sampling and selection within demes
    for (i = 0; i < L; i++) {
        n[i] = gsl_ran_binomial(R, xmig[i] + s * xmig[i] * (1 - xmig[i]), N);
    }

}

/*function roundint to round a double to an integer*/
int roundint(double x) {
    int result;

    result = floor(x);
    if (x - result >= 0.5) result++;

    return result;
} /*roundint*/



int main(int argc, char  *argv[]) {
    double mig, s; // migration rate, selection
    unsigned int N, L, t, tfinal, l0; // deme size, number of demes, time, max number of generations, origin of the sweep
    FILE *datafile, *paramfile;
    char *outfile = malloc(1000);
    char *outfile1 = malloc(1000);
    unsigned int i, j, ntot, m;
    //Initialize variables:
    if (argc != 7) {
        printf("Usage: L N s m tfinal l0\n");
        return 1;
    }

    j = 1;
    L = (unsigned int) roundint(atof(argv[j++]));
    N = (unsigned int) roundint(atof(argv[j++]));
    s = atof(argv[j++]);
    mig = atof(argv[j++]);
    tfinal = (unsigned int) roundint(atof(argv[j++]));
    l0 = (unsigned int) roundint(atof(argv[j++]));


    sprintf(outfile1, "L=%u_N=%u_s=%f_m=%f_l0=%u_%d", L, N, s, mig, l0, m);

    // Print out variables to parameter file:
    paramfile = fopen(strcat(outfile1, "_params.txt"), "w");
    fprintf(paramfile, "L = %u\nN = %u\ns = %f\nm = %f\ntfinal = %u\nl0 = %u", L, N, s, mig, tfinal, l0);
    fclose(paramfile);


    for (m = 0; m < 5; m++) {
        sprintf(outfile, "L=%u_N=%u_s=%f_m=%f_l0=%u_%d", L, N, s, mig, l0, m);
        outfile = strcat(outfile, ".txt");


        // gsl random setup:
        gsl_rng_env_setup();
        T = gsl_rng_default;
        R = gsl_rng_alloc (T);
        gsl_rng_set(R,m);


        // Initialize population:
        unsigned int n[L];
        // Leftmost demes fixed for sweeping allele
        // Fill enough demes so that it is very unlikely to die out:
        if (1 + 10 / (N * s) > L / 2){
            printf("Warning: meta-population is too small for selection to be strong.\n");
            return 1;
        }
        for (i = 0; i < L; i++){
            n[i] = 0;
        }
        

        while (n[l0] == 0){
            for (i = 0; i < L; i++){
                n[i] = 0;
            }
            n[l0] = 1;
            //Open the datafile for recording:
            datafile = fopen(outfile, "w");
            // Run until alleles fixes (or goes extinct):
            for (t = 0; t < tfinal; t++) {

                // Record the status of the population and check for fixation:
                ntot = 0;
                for (i = 0; i < L; i++){
                    fprintf(datafile, " %d", n[i]);
                    ntot += n[i];
                }
                fprintf(datafile, "\n");

                // Stop if one of the alleles is fixed or extinct:
                if ((ntot == 0) || (ntot == N * L)) {break;}

                // Evolve for one generation
                next_gen(n, s, mig, L, N);

            }

            fclose(datafile);
        }
    }




    if (t == tfinal){
        printf("Simulation finished without fixation.\n");
        return 1;
    }

    // if (n[0] == 0){
    //  printf("Allele went extinct.\n");
    //  return 2;
    // }


    return 0;

}
