/*
 *  Simple simulation of a sweep in space
 *
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>



// For misha's laptop, create make file
// run by typing : gcc Sweep_in-2D_version-m.c -o plain_sweep -lgsl -lgslcblas -lm
// ./plain_sweep L N s m tfinal l0



//Todo : when beneficial mutation is fixed in 10% of the demes, save n[i] and t. Get the right most deme with a beneficial allele, choose a random individual among them, and give a neutral mutation and track it.
// When the neutral mutation goes extinct, set n[i] and t to the saved values and start tracking again. Save n[i] and n_neutral[i]'s.
// Later, there should be one more input value (argv) that sets how far out from the center of mass (i_cm = sum(n[i] * i) / N ) to put in a neutral mutation.

// Global variables
const gsl_rng_type * T;
gsl_rng * R;

void next_gen(int L, unsigned int n[L][L], double s, double mig, unsigned int N) {

    double x[L][L]; // frequencies
    double xmig[L][L]; // after migration
    int i;
    int j;

    // Initialize frequencies:
    
    for (i = 0; i < L; i++)
    {
        for(j = 0; j < L; j++)
        {
            x[i][j] = (double)n[i][j] / N;
        }
    }

    // Migration:
    for (i = 1; i < L - 1; i++)
    {
        for(j = 1; j < L - 1; j++)
        {
            xmig[i][j] = x[i][j] + mig * (0.25 * (x[i-1][j] + x[i+1][j] + x[i][j-1] + x[i][j+1] ) - x[i][j]);
        }
    }

    //Boundary Conditions (PBC)
    for (i = 1; i < L - 1; i++)
    {
        xmig[0][i] = x[0][i] + mig * 0.25 * (x[1][i] + x[L-1][i] + x[0][i-1] + x[0][i+1]) - mig*x[0][i];
        xmig[L-1][i] = x[L-1][i] + mig * 0.25 * (x[L-2][i] + x[0][i] + x[L-1][i-1] + x[L-1][i+1]) - mig*x[L-1][i];
        xmig[i][0] = x[i][0] + mig * 0.25 * (x[i][1] + x[i][L-1] + x[i-1][0] + x[i+1][0]) - mig*x[i][0];
        xmig[i][L-1] = x[i][L-1] + mig * 0.25 * (x[i][0] + x[i][L-2] + x[i-1][L-1] + x[i+1][L-1]) - mig*x[i][L-1];

    }
    
    xmig[0][0] = x[0][0] + mig * 0.25 * (x[1][0] + x[L-1][0]+ x[0][1] + x[0][L-1]) - mig*x[0][0];  // Origin at bottom right
    xmig[L-1][0] = x[L-1][0] + mig * 0.25 * (x[L-2][0] + x[0][0] + x[L-1][1] + x[L-1][L-1]) - mig*x[L-1][0];   //Top right
    xmig[0][L-1] = x[0][L-1] + mig * 0.25 * (x[0][L-2] + x[0][0] + x[1][L-1] + x[L-1][L-1])- mig*x[0][L-1];   //Bottom Left
    xmig[L-1][L-1] = x[L-1][L-1] + mig * 0.25 * (x[L-2][L-1] + x[0][L-1] +  x[L-1][L-2] + x[L-1][0]) - mig*x[L-1][L-1];  //Top Left

    // Sampling and selection within demes
    for (i = 0; i < L; i++) 
    {

    for(j = 0; j < L; j++)
        {
            n[i][j] = gsl_ran_binomial(R, xmig[i][j] + s * xmig[i][j] * (1 - xmig[i][j]), N);  // Logistic Growth and Sampling within a deme 
        }
        
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

    //sprintf(outfile1, "L=%u_N=%u_s=%f_m=%f_tfinal=%u_l0=%u_%d", L, N, s, mig, tfinal, l0, m);
    //Print out variables to parameter file:
    //paramfile = fopen(strcat(outfile1, "_params.txt"), "w");
    //fprintf(paramfile, "L = %u\nN = %u\ns = %f\nm = %f\ntfinal = %u\nl0 = %u", L, N, s, mig, tfinal, l0);
    //fclose(paramfile);


    for (m = 0; m < 1; m++) {
        sprintf(outfile, "L=%u_N=%u_s=%0.3f_m=%0.2f_l0=%u", L, N, s, mig, l0);
        outfile = strcat(outfile, ".txt");
        
        // gsl random setup:
        gsl_rng_env_setup();
        T = gsl_rng_default;
        R = gsl_rng_alloc (T);
        gsl_rng_set(R,m);


        // Initialize population:
        unsigned int n[L][L];
        // Leftmost demes fixed for sweeping allele
        // Fill enough demes so that it is very unlikely to die out:
        if (1 + 10 / (N * s) > (L*L) / 4)
        {
            printf("Warning: meta-population is too small for selection to be strong.\n");
            return 1;
        }


        for (i = 0; i < L; i++)
        {
            for(int k = 0; k < L; k++)
            {
                n[i][k] = 0;
            }
        }

        //Note: 0 is wildtype, 1 is neutral mutation 

        while (n[l0][l0] == 0)
        {
            for (i = 0; i < L; i++)
            {
                for(int k = 0; k < L ; k++)
                { 
                    n[i][k] = 0;
                }
            }
            n[l0][l0] = 1;

            //Open the datafile for recording:
            datafile = fopen(outfile, "w");
            // Run until alleles fixes (or goes extinct):
            for (t = 0; t < tfinal; t++) {

                // Record the status of the population and check for fixation:
                ntot = 0;
                for (i = 0; i < L; i++)
                {
                    for(int k=0; k < L; k++)
                    {    
                        fprintf(datafile, " %d", n[i][k]);
                        ntot += n[i][k];
                    }
                    fprintf(datafile, "\n");
                }
                fprintf(datafile, "\n");

                // Stop if one of the alleles is fixed or extinct:
                if ((ntot == 0) || (ntot == N * L * L)) {break;}

                // Evolve for one generation
                next_gen(L, n, s, mig, N);
                    printf("Current Time = %d \n",t);


            }
            fclose(datafile);
        }
    }

    printf("Tfinal = %d",t);

    //Sanity Check
    if (t == tfinal){
        printf("Simulation finished without fixation.\n");
        return 1;
    }

    return 0;

}







