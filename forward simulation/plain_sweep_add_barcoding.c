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

// run by typing : gcc plain_sweep_add_barcoding.c -o plain_sweep_barcode -lgsl -lm
// ./plain_sweep_barcode.exe L N s m tfinal tlabel nbarcode seed

// At some time tlabel, label 10 individuals at fround (l = 2 * sum(n[i] * i) / N). nlabel is 11 * L long, i * L to (i + 1) L - 1 entries being the number of individuals with i-th label. 
// i goes from 0 to 9. 10 L to 11L - 1 entries are for the unlabeled population. 
// Note that this is a multinomial sampling! (sum of 11 cannot be bigger than N)


// Todo1 : change next_gen so that it goes through 11 'rows' separately (for loop)
// Todo2 : change main so that it saves n_label
// Todo3 : calculate the location of the wavefront when t == tlabel, 
// Todo4 : while loop should be sum of nlabel[iL] == 0

//Todo : when beneficial mutation is fixed in 10% of the demes, save n[i] and t. Get the right most deme with a beneficial allele, choose a random individual among them, and give a neutral mutation and track it.
// When the neutral mutation goes extinct, set n[i] and t to the saved values and start tracking again. Save n[i] and n_neutral[i]'s.
// Later, there should be one more input value (argv) that sets how far out from the center of mass (i_cm = sum(n[i] * i) / N ) to put in a neutral mutation.

// Global variables
const gsl_rng_type * T;
gsl_rng * R;

void next_gen(unsigned int n[], double s, double mig, int L, unsigned int N, unsigned int nbarcode) {

    double x[(nbarcode + 1) * L]; // frequencies
    double xmig[(nbarcode + 1) * L]; // after migration
    unsigned int ntotal[L]; // number of individuals with beneficial allele in each deme
    double xmigtotal[L];
    int i;
    int row;

    // Initialize frequencies:
    for (i = 0; i < nbarcode * L; i++){
        x[i] = (double)n[i] / N;
    }

    for (row = 0; row < nbarcode; row++){
        // Migration:
        for (i = 1; i < L - 1; i++){
            xmig[row * L + i] = x[row * L + i] + mig * (0.5 * (x[row * L + i - 1] + x[row * L + i + 1]) - x[row * L + i]);
        }
        xmig[row * L] = x[row * L] + mig * 0.5 * (x[row * L + 1] - x[row * L]);
        xmig[row * L + L - 1] = x[row * L + L - 1] + mig * 0.5 * (x[row * L + L - 2] - x[row * L + L - 1]);

    }

    for (i = 0; i < L; i++){
        ntotal[i] = 0;
        xmigtotal[i] = 0;
        for (row = 0; row < nbarcode; row++) {
            ntotal[i] += n[row * L + i];
            xmigtotal[i] += xmig[row * L + i];
        }
    }

    // Sampling and selection within demes
    for (i = 0; i < L; i++) {
        for (row = 0; row < nbarcode; row++) {
            n[row * L + i] = gsl_ran_binomial(R, xmig[row * L + i] + s * xmig[row * L + i] * (1 - xmigtotal[i]), N);
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
    unsigned int N, L, t, tfinal, tlabel, nbarcode, seed; // deme size, number of demes, time, max number of generations, generation when 10 individuals at front are barcoded, seed for rand number
    FILE *datafile, *paramfile;
    char *outfile = malloc(1000);
    char *outfile1 = malloc(1000);
    unsigned int i, j, ntot, row;
    //Initialize variables:
    if (argc != 9) {
        printf("Usage: L N s m tfinal tlabel nbarcode seed\n");
        return 1;
    }

    j = 1;
    L = (unsigned int) roundint(atof(argv[j++]));
    N = (unsigned int) roundint(atof(argv[j++]));
    s = atof(argv[j++]);
    mig = atof(argv[j++]);
    tfinal = (unsigned int) roundint(atof(argv[j++]));
    tlabel = (unsigned int) roundint(atof(argv[j++]));
    nbarcode = (unsigned int) roundint(atof(argv[j++]));
    seed = (unsigned int) roundint(atof(argv[j++]));


    sprintf(outfile1, "L=%u_N=%u_s=%f_m=%f_tfinal=%u_tlabel=%u_nbarcode=%u", L, N, s, mig, tfinal, tlabel, nbarcode);

    // Print out variables to parameter file:
    paramfile = fopen(strcat(outfile1, "_params.txt"), "w");
    fprintf(paramfile, "L = %u\nN = %u\ns = %f\nm = %f\ntfinal = %u\ntlabel=%u\nnbarcode=%u\n", L, N, s, mig, tfinal, tlabel, nbarcode);
    fclose(paramfile);


    sprintf(outfile, "L=%u_N=%u_s=%f_m=%f_tfinal=%u_tlabel=%u_nbarcode=%u", L, N, s, mig, tfinal, tlabel, nbarcode);
    outfile = strcat(outfile, ".txt");


    // gsl random setup:
    gsl_rng_env_setup();
    T = gsl_rng_default;
    R = gsl_rng_alloc (T);
    gsl_rng_set(R, seed);


    // Initialize population:
    unsigned int n[L];
    // labeled population (10 labels)
    unsigned int nlabel[nbarcode * L];
    // Leftmost demes fixed for sweeping allele
    // Fill enough demes so that it is very unlikely to die out:
    if (1 + 10 / (N * s) > L / 2){
        printf("Warning: meta-population is too small for selection to be strong.\n");
        return 1;
    }

    for (i = 0; i < L; i++){
        n[i] = 0;
    }
    for (i = 0; i < nbarcode * L; i++){
        nlabel[i] = 0;
    }

    unsigned int n0 = 0;
    for (i = 0; i < nbarcode; i++){
        n0 += nlabel[i * L];
    } 


    while (n0 == 0){
        for (i = 0; i < nbarcode * L; i++){
            nlabel[i] = 0;
        }
        nlabel[0] = 1;
        
        n0 = 0;
        for (i = 0; i < nbarcode; i++){
            n0 += nlabel[i * L];
        } 
        //Open the datafile for recording:
        datafile = fopen(outfile, "w");
        // Run until alleles fixes (or goes extinct):
        for (t = 0; t < tfinal; t++) {

            // Record the status of the population and check for fixation:
            ntot = 0;
            n0 = 0;
            for (i = 0; i < nbarcode * L; i++){
                fprintf(datafile, " %d", nlabel[i]);
                ntot += nlabel[i];
            }
            fprintf(datafile, "\n");
            for (i = 0; i < nbarcode; i++) {
                n0 += nlabel[i * L];
            }
            // Stop if one of the alleles is fixed or extinct:
            if ((ntot == 0) || (ntot == N * L)) {break;}
            if (t == tlabel) {
                unsigned int lfront = 0;
                for (i = 0; i < L; i++) {
                    lfront += roundint(2 * nlabel[i] * i / N);
                }
                if (nbarcode < nlabel[lfront]) {
                    for (row = 1; row < nbarcode; row++) {
                        nlabel[row * L + lfront] = 1;
                        nlabel[lfront] -= 1;
                    }

                }
                else {
                    for (row = 1; row < nlabel[lfront]; row++){
                        nlabel[row * L + lfront] = 1;
                    }
                    for (row = nlabel[lfront]; row < nbarcode; row++){
                        nlabel[row * L + lfront - 1] = 1;
                        nlabel[lfront - 1] -= 1;
                    }
                    nlabel[lfront] = 0;
                }
            }


            // Evolve for one generation
            next_gen(nlabel, s, mig, L, N, nbarcode);

        }

        fclose(datafile);
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
