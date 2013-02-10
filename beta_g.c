/*#include <stdlib.h>*/
#include <math.h>

void dist_beta(double **X, double **Y, int genes, int nsamplex, int nsampley, double *beta_g) {
    int g; 
    for ( g = 0; g < genes; g++) {
        double x = 0, y = 0, z = 0;
        int i, j;

          for ( i = 0; i < nsamplex - 1; i++) {
              for ( j = i + 1; j < nsamplex; j++) {
                  x += (X[g][i] - X[g][j])*(X[g][i] - X[g][j]);
              }
          }
          x = x /( nsamplex * ( nsamplex - 1. ) );
  
          for ( i = 0; i < nsampley - 1; i++) {
              for ( j = i + 1; j < nsampley; j++) {
                  y += (Y[g][i] - Y[g][j])*(Y[g][i] - Y[g][j]);
              }
          }
          y = y /( nsampley * ( nsampley - 1. ) );
  
          for ( i = 0; i < nsampley; i++) {
              for ( j = 0; j < nsamplex; j++) {
                  z += (Y[g][i] - X[g][j])*(Y[g][i] - X[g][j]);
              }
          }
          z = z / ( nsamplex * nsampley );

        beta_g[g] = x + y - z;
    }

}
