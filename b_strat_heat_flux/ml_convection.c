/**
  Field scale wave breaking (multilayer solver)

No stratification


USAGE

  * Compilation and run (with mpirun)
  make

  * HPC source file generation
  make hpc

  * Running on HPC

  -> first compile on hpc with
  mpicc -std=c99 -O2 _*.c

  -> then run (in a slurm file)
  srun ./ml_convection


HOW TO CREATE A RESTART
  ncks -d time,N out.nc restart.nc

  with N = index of timestep in the file "out.nc"
*/


#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
//#include "aaubert/layered_perso/remap.h"
// #include "layered/perfs.h"
#include "bderembl/libs/extra.h"      // parameters from namlist
//#include "bderembl/libs/netcdf_bas.h" // read/write netcdf files
#include "bderembl/libs/netcdf_bas.h"
//#include "view.h" // Basilisk visualization
// #include "display.h"
#define g_ 9.81





/*

TO DO : 
- add dimensions so that they can be saved in the netcdf
- add ability to restart from netcdf file

*/

char my_name[40] = "ml_breaking";

/*
DEFAULT PARAMETERS

Dimensions : [Length, Time]
*/
char namlist[80] = "namelist.toml";    // file name of namlist
char file_out[20] = "out.nc";         // file name of output
char file_restart[20] = "restart.nc"; // file name of restart
// -> Initial conditions
double strat = 0.002;       // [s-1] N^2 stratification
double qt = 100;            // [W] Heat flux
double Ts = 20;             // [K] Surface temperature (arbitrary)
// -> Domain definition
int N_grid = 5 [];       // 2^N_grid : number of x and y gridpoints
double L = 200.0 [1];       // domain size
int N_layer = 2 [];         // number of layers
double h0 = 1.0 [1];        // depth of water
// -> Runtime parameters
int restart = 0;            // 1: restart, 0: no restart
double tend = 2.0;          // end time of simulation
// -> saving outputs
double dtout = 2.0;        // dt for output in netcdf
int pad = 4;                // number of 0-padding for ouput files
int nout = 1;               // number of the outfile
char fileout[100];          // name of outfile
// -> physical properties
double nu0 = 0.;            // Viscosity for vertical diffusion
double thetaH = 0.5;        // theta_h for dumping fast barotropic modes


double rho0 = 1025.     // [kg.m-3] reference density
double cp = 4.2e3       // [J.kg-1.K-1] heat capacity water
double betaT = 2e-4;    // Linear equation of state: drho = betaT*(T0-T) (Vallis 2.4)
#define drho(T) (betaT*(T0-T);
#define T0(z) (T1 + (T0 - T1)*(z + H0)/H0)
#include "layered/dr.h"

// diag
//double *dudz, *eps, *u_profile;
double *u_profile;
double dt_mean = 0.4;
double U=0.0;
int l;
static FILE * fp;



int main(int argc, char *argv[])  
{
 
  // Building a 'params' array with all parameters from the namlist
  params = array_new();
  add_param("N_grid", &N_grid, "int");
  add_param("L", &L, "double");
  add_param("N_layer", &N_layer, "int");
  add_param("h0", &h0, "double");
  add_param("tend", &tend, "double");
  add_param("nu0", &nu0, "double");
  add_param("thetaH", &thetaH, "double");
  add_param("restart", &restart, "int");
  add_param("dtout", &dtout, "double");
  
  // Search for the configuration file with a given path or read params.in
  if (argc == 2)
    strcpy(file_param, argv[1]);
  else
    strcpy(file_param, namlist);
  read_params(file_param);
  
  
  // Settings solver values from namlist values
  L0 = L;
  nu = nu0;
  N = 1 << N_grid; // 1*2^N_grid
  nl = N_layer;
  G = g_;
  theta_H = thetaH;
  CFL_H = 1; // Smaller time step

  // Boundary condition
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  
  // allocate diag
  // dudz = (double*)malloc(nl * sizeof(double));
  // eps = (double*)malloc(nl * sizeof(double));
  u_profile = (double*)calloc(nl, sizeof(double));
  l=nl-1;
  fprintf(stderr, "u%d ini = %f\n",l, U);
  fp  = fopen("u_profile.dat","w"); // reset file
  fclose(fp);
  
  fprintf (stderr, "Read in parameters!\n");
  run();
}

event init(i =  0) {
  geometric_beta (1./3., true); // Varying layer thickness
  if (restart!=1) {

    // step 1: set eta and h
    foreach() {
      zb[] = -h0;
      eta[] = 0.;
      double H = - zb[];
      foreach_layer() {
        h[] = H/nl;
      } 
    }
    
    // step 2: remap
    vertical_remapping (h, tracers);
    // step 3: set currents
    foreach() {
      foreach_layer() {
        u.x[] = 0.;
        u.y[] = 0.;
        w[] = 0.;
      }
    }
  }
  else {
    //fprintf(stderr, "restart = %d\n",restart);
    // Restarting from netcdf ...
    // "test, %s", strcat(strcat(my_name,"/"),file_restart));
    // read_nc({zb, h, u, w}, strcat(strcat(my_name,"/"),file_restart));
    fprintf(stderr, "Restarting : reading from file\n"); 
    fprintf(stderr, "->no yet coded ...\n");
    // read_nc({zb, h, u, w}, file_restart);
    //
  }


  fprintf (stderr,"Done initialization!\n");
  create_nc({zb, h, u, w, eta}, file_out);

}

event forcing(i++){
  // This event adds a heat flux forcing at surface
  foreach(){
    T[:,:,nl-1] = dt*(T[:,:,nl-1] + qt/h[:,:,nl-1])
  }

}

// This event compute layer average of u.x
event compute_horizontal_avg (i++; t<=tend+1e-10){
  foreach(reduction(+:u_profile[:nl]))
    foreach_layer(){
      u_profile[point.l] += u.x[] / (N*N) * dt / dt_mean;
    }
}

// This even writes to a file the layer average
event write_diag(t=0., t+=dt_mean){
    // main worker is writing the file
    if (pid()==0) {
      fp  = fopen("u_profile.dat","a");
      if (fp == NULL){
        fprintf(stderr, "Error opening file u_profile.dat");
        return 2;
      }
      for (int i=0; i<nl; ++i) {
        fprintf (fp, "%f %d %g\n", t, i, u_profile[i]);
      }
      fprintf(fp,"\n");
      fclose(fp);
    }
    // Reset the profile for all workers
    for (int i=0; i<nl; ++i) {
      u_profile[i] = 0.0;
    }


}

// Writing a 4D netcdf file
event output(t = 0.; t<= tend+1e-10; t+=dtout){
  write_nc();
}

// Clean for my diag of layer avg
event cleanup(t=end){
  free(u_profile);
}

/**
Results: plots
~~~gnuplot 
# Plot the heatmap
set pm3d map
set view map
set xlabel "Time (s)"
set ylabel "Layer"
set cblabel "Value"
set yrange [0:15]
set terminal pngcairo size 800,600 enhanced font 'Verdana,12'
set output 'u_profile.png'
set size 0.9, 0.9
splot "u_profile.dat" using 1:2:3 with pm3d
unset output
~~~
**/

