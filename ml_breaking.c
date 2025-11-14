/**
  Field scale wave breaking (multilayer solver)

No stratification


USAGE

  * Compilation:
  make

  * running:
    make run
  or
    ./ml_breaking config_name

  * restart:
    make restart
    then edit the .toml config file in yourfolder_restart/

  -> if config_name is not given, the default is "namlist.toml"


HOW TO CREATE A RESTART
  ncks -d time,N out.nc restart.nc

  with N = index of timestep in the file "out.nc"
*/

#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
// #include "layered/perfs.h"
#include "bderembl/libs/extra.h"      // parameters from namlist
#include "bderembl/libs/netcdf_bas.h" // read/write netcdf files
//#include "netcdf_bas.h"
#include "view.h" // Basilisk visualization

#define g_ 9.81
#include "spectrum.h" // Initial conditions generation


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
char namlist[80] = "namlist.toml";    // file name of namlist
char file_out[20] = "out.nc";         // file name of output
char file_restart[20] = "restart.nc"; // file name of restart
// -> Initial conditions
double P = 0.2 [1, -1];     // energy level (estimated so that kpHs is reasonable)
int coeff_kpL0 = 10 [];     // kpL0 = coeff_kpL0 * pi
int N_mode = 32 [];         // Number of modes in wavenumber space
int N_power = 5 [];         // directional spreading coeff
int F_shape = 0 [];         // shape of the initial spectrum
// -> Domain definition
int N_grid = 5 [];       // 2^N_grid : number of x and y gridpoints
double L = 200.0 [1];       // domain size
int N_layer = 2 [];         // number of layers
double kp = PI*10/200.0 [-1];// peak wave number
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
int RANDOM = 2;             // For random number generator
double thetaH = 0.5;        // theta_h for dumping fast barotropic modes


// diag
//double *dudz, *eps, *u_profile;
double *u_profile;
static FILE * fp;




// void my_gradzl(scalar * s, vector * v){
//   // Gradient in multilayer
//   assert (list_len(f) == vectors_len(g));
//
//   scalar dsdx, dsdy, dsdzl, dzdzl; // temp var
//
//   foreach(){
//     foreach_layer(){
//       dsdx[] = (s[1,0] - s[-1,0])/Delta;
//       dsdy[] = (s[0,1] - s[0,-1])/Delta;
//       dzdx[] = (h[point.l])
//       }
//     }
//   }













int main(int argc, char *argv[])  
{

  // Building a 'params' array with all parameters from the namlist
  params = array_new();
  add_param("P", &P, "double");
  add_param("coeff_kpL0", &coeff_kpL0, "int");
  add_param("N_mode", &N_mode, "int");
  add_param("N_power", &N_power, "int");
  add_param("F_shape", &F_shape, "int");
  add_param("N_grid", &N_grid, "int");
  add_param("L", &L, "double");
  add_param("N_layer", &N_layer, "int");
  add_param("h0", &h0, "double");
  add_param("tend", &tend, "double");
  add_param("nu0", &nu0, "double");
  add_param("RANDOM", &RANDOM, "int");
  add_param("thetaH", &thetaH, "double");
  add_param("restart", &restart, "int");
  add_param("dtout", &dtout, "double");
  kp = PI * coeff_kpL0 / L; // kpL=coeff x pi peak wavelength
  
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
  u_profile = (double*)malloc(nl * sizeof(double));
  fp  = fopen("u_profile.dat","w"); // reset file
  fclose(fp);

  fprintf (stderr, "Read in parameters!\n");
  run();
}

event init(i =  0) {
  // TO DO
  geometric_beta (1./3., true); // Varying layer thickness
  if (restart!=1) {
  // Generate linealy spaced kx, ky according to specified # of modes, and
  //  interpolated F(kx,ky)
  T_Spectrum spectrum;
  spectrum = spectrum_gen_linear(N_mode, N_power, L, P, kp);
  
  foreach() {
    zb[] = -h0;
    eta[] = wave(x, y, N_grid, spectrum);
    double H = wave(x, y, N_grid, spectrum) - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H*beta[point.l];
      z += h[]/2.;
      u.x[] = u_x(x, y, z, N_grid, spectrum);
      u.y[] = u_y(x, y, z, N_grid, spectrum);
      w[] = u_z(x, y, z, N_grid, spectrum);
      z += h[]/2.;
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
  // sprintf(fileout, "%0*d.nc", pad, 0); // add the padding
  create_nc({zb, h, u, w, eta}, file_out);
  write_nc(); 
}


// event film surface




event compute_horizontal_avg (i++; t<= tend+1e-10){
  fprintf(stderr, "dt %f\n",dt);
  foreach(reduction(+:u_profile[:nl])){
    foreach_layer(){

      // fprintf(stderr, "t %f, i %d, value %f\n", t, i,  u.x[0,0,_layer] / (N*N) * dt / 1.0);
      u_profile[point.l] += u.x[] / (N*N) * dt / 1.0;
    }
  }
  fprintf(stderr, "t %f, i %d, value %f\n", t, i,  u_profile[nl-2]);

}
event write_diag(t=0., t+=1.){
  #pragma omp critical
  {
  fp  = fopen("u_profile.dat","a");
  fprintf(fp, "%f ",t);
  for (int i=0; i<nl; ++i) {
    fprintf (fp, " %g", u_profile[i]);
  }
  fprintf(fp,"\n");
  fclose(fp);
  }
}
  // vector gradU[];
  // vector gradV[];
  // dudx = u.x[]
  //
  //
  // gradient(u.x, gradU);
  // gradient(u.y, gradV);
  //
  //
  //
  // foreach(reduction(+:dudz), reduction(+:eps)){
  //   foreach_layer(){
  //
  //   }
  // }
 

// event image (t = end) {
//   clear();
//   static FILE * fp = fopen ("image.ppm", "w");
//   output_ppm(eta, fp, min=-0.1, max=0.1);
//
// }

event snapshot (t = end)
{
  clear();
	//  view (fov = 17.3106, quat = {0.475152,0.161235,0.235565,0.832313},
	// tx = -0.0221727, ty = -0.0140227, width = 1200, height = 768);
  view( fov = 44, camera = "back", width = 1200, height = 768);
  char s[80];
  sprintf (s, "t = %.2f", t);
  draw_string (s, size = 80);
  for (double x = -1; x <= 1; x++)
    translate (x) {
      squares ("eta", linear = true, z = "eta", min = -0.15, max = 0.6, map=gray);
    }
  colorbar(map=gray, label="eta (m)", min=-1.0,max=2.0);
  save ("snap.png");
}


event output(t = 0.; t<= tend+1e-10; t+=dtout){
  fprintf(stdout, "output at t=%f, i=%d\n", t, i);
  write_nc();
}

event cleanup(t=end){
  free(u_profile);
}
//
// event end (t = tend) {
//   fprintf (stderr,"tend = %f \n", tend);
//   // create_nc({zb, h, u, w}, "end.nc");
//   write_nc();
//   // To Fix 31/10/25: This output is filled with nan. 
//
// }
