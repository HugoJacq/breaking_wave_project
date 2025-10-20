#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "interpolate.h"



int main() {
    int N = 10;
    int Nx = 1000;
    int Ny = 1000;
    double r = 1.0;


    double xp[N];
    double yp[N];
    double fp[N];

    // Generate points in a circle
    for (int i = 0; i < N; ++i) {
        double theta = 2 * PI * i / N;
                xp[i] = r * cos(theta);
        yp[i] = r * sin(theta);
        fp[i] = sin(theta);  // some fake spectrum value
    }

    // double xi[Nx*Ny];
    // double yi[Ny*Ny];
    // double F_interp[Nx*Ny];
    //
    // for (int i = 0; i < Ny; ++i) {
    //     for (int j = 0; j < Nx; ++j) {
    //         xi[i * Nx + j] = -1.0 + 2.0 * j / (Nx - 1); // from -1 to 1
    //         yi[i * Ny + j] = -1.0 + 2.0 * i / (Ny - 1); // from -1 to 1
    //     }
    // }

    // interp 
    double Fi;
    double xi[Nx];
    double yi[Ny];
    for (int i=0; i<Nx; ++i){
      xi[i] = -1.0 +2.0*i / (Nx-1);
      yi[i] = -1.0 +2.0*i / (Nx-1);
    } 
    double F[Nx*Ny];
    for (int i=0; i<Nx;++i){
      for (int j=0; j<Ny;++j){
        F[i*Nx+j] = xi[i]*xi[i] + yi[i]*yi[i];
        //F[i*Nx+j] = 0.5*( xi[i] + yi[j] );
      } 
    }

    Fi = interp_lin(xi, yi, Nx, Ny, 0.05, 0.05, F); //, double F[])
    
    printf("true result F(0.11,0.14)=%f\n",0.05*0.05+0.05*0.05);
    //printf("true result F(0.11,0.14)=%f\n",0.5*(0.11+0.14));

    printf("interpolated Fi = %f\n",Fi);
  
  
    // Print a few results
    // for (int i = 0; i < Ny; i += 10) {
    //     for (int j = 0; j < Nx; j += 10) {
    //         printf("F[%d][%d] = %f\n", i, j, F_interp[i * Nx + j]);
    //     }
    // }

    return 0;
}
