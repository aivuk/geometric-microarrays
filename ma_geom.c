#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_vector.h> 
#include <gsl/gsl_randist.h> 
#include <gsl/gsl_combination.h>  
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h> 

gsl_vector *calc_metric_tensor(gsl_matrix *gs, gsl_vector *c1, gsl_vector *c2) {

    gsl_vector *eg = gsl_vector_alloc(gs->size2);
    gsl_combination *comb1, *comb2;
    int gi, i1, i2, exp1i, exp2i;       
    double eg_f, eg_c, egi, nf, nc;

    nf = 1/(c1->size*c2->size);
    nc = 1/(gsl_sf_choose(c1->size, 2) + gsl_sf_choose(c2->size, 2));

    comb1 = gsl_combination_alloc(c1->size, 2);
    comb2 = gsl_combination_alloc(c2->size, 2); 

    for (gi = 0; gi < gs->size2; ++gi) {
        eg_c = 0;
        eg_f = 0;
        egi = 0;
        // Calculate the energy part from same classes experiments 

        // Class 1 near experiments 

        gsl_combination_init_first(comb1);

        do {
            i1 = gsl_combination_get(comb1, 0);
            i2 = gsl_combination_get(comb1, 1);

            eg_c += abs(gsl_matrix_get(gs, gsl_vector_get(c1, i1), gi) 
                        - gsl_matrix_get(gs, gsl_vector_get(c1, i2), gi));

        } while (gsl_combination_next (comb1) == GSL_SUCCESS);

        // Class 2 near experiments  

        gsl_combination_init_first(comb2); 
        
        do {
            i1 = gsl_combination_get(comb2, 0);
            i2 = gsl_combination_get(comb2, 1);

            eg_c += abs(gsl_matrix_get(gs, gsl_vector_get(c2, i1), gi) 
                        - gsl_matrix_get(gs, gsl_vector_get(c2, i2), gi));

        } while (gsl_combination_next (comb2) == GSL_SUCCESS);
 
        // Calculate the energy part from different classes experiments

        for (i1 = 0; i1 < c1->size; ++ i1) {
            exp1i = gsl_vector_get(c1, i1);
            for (i2 = 0; i2 < c2->size; ++ i2) { 
                exp2i = gsl_vector_get(c2, i2); 
                eg_f += abs(gsl_matrix_get(gs, exp1i, gi) - gsl_matrix_get(gs, exp2i, gi));
            }
        }

        egi = nc*eg_c - nf*eg_f;
        gsl_vector_set(eg, gi, egi);
    }

    gsl_combination_free(comb1);
    gsl_combination_free(comb2); 

    return eg;
     
}

int main(int argc, char **argv) {

    int i, j;
    int ngenes = 23000;
    int nexp = 100;
    gsl_matrix *gexp = gsl_matrix_alloc(nexp, ngenes);
    gsl_vector *eg;
    gsl_vector *c1 = gsl_vector_alloc(50);
    gsl_vector *c2 = gsl_vector_alloc(50); 
    gsl_rng *r = gsl_rng_alloc (gsl_rng_mt19937);

    gsl_rng_env_setup();

    for (i = 0; i < nexp; ++i) {
        for (j = 0; j < ngenes; ++j) { 
            gsl_matrix_set(gexp, i, j, gsl_ran_gaussian(r, 1)); 
        }
    }

    for (i = 0; i < 50; ++i) { 
        gsl_vector_set(c1, i, i);
        gsl_vector_set(c2, i, 50 + i); 
    }

    eg = calc_metric_tensor(gexp, c1, c2);

    for (i = 0; i < ngenes; ++i) { 
        printf("%lf\n", gsl_vector_get(eg, i));
    }
 
    return 0;
}
