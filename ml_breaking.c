/**
# Field scale wave breaking (multilayer solver)
*/

#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
// #include "layered/perfs.h"
#include "bderembl/libs/extra.h"      // parameters from namlist
#include "bderembl/libs/netcdf_bas.h" // read/write netcdf files
#define g_ 9.81
#include "spectrum.h" // Initial conditions generation


/*

TO DO : 
- add dimensions so that they can be saved in the netcdf
- add ability to restart from netcdf file

*/



/*
DEFAULT PARAMETERS

Dimensions : [Length, Time]
*/
char namlist[80] = "namlist.toml"; // file name of namlist
char file_out[20] = "out.nc";
// -> Initial conditions
double P = 0.2 [1, -1];     // energy level (estimated so that kpHs is reasonable)
int coeff_kpL0 = 10 [];     // kpL0 = coeff_kpL0 * pi
int N_mode = 32 [];         // Number of modes in wavenumber space
int N_power = 5 [];         // directional spreading coeff
int F_shape = 0 [];         // shape of the initial spectrum
// -> Domain definition
int N_grid = 1024 [];       // number of x and y gridpoints
double L = 200.0 [1];       // domain size
int N_layer = 2 [];         // number of layers
double kp = PI*10/200.0 [-1];// peak wave number
double h0 = 1.0 [1];        // depth of water
// -> Runtime parameters
int restart = 0;            // 1: restart, 0: no restart
double tend = 2.0;          // end time of simulation
// -> saving outputs
int pad = 4;                // number of 0-padding for ouput files
int nout = 1;               // number of the outfile
char fileout[100];          // name of outfile
// -> physical properties
double nu0 = 0.;            // Viscosity for vertical diffusion
int RANDOM = 2;             // For random number generator
double thetaH = 0.5;        // theta_h for dumping fast barotropic modes






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
  add_param("L", &N_grid, "double");
  add_param("N_layer", &N_layer, "int");
  add_param("h0", &h0, "double");
  add_param("tend", &tend, "double");
  add_param("nu0", &nu0, "double");
  add_param("RANDOM", &RANDOM, "int");
  add_param("thetaH", &thetaH, "double");
  add_param("restart", &restart, "int");
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
  N = 1 << 7; // idk what this is used for 
  nl = N_layer;
  G = g_;
  theta_H = thetaH;
  CFL_H = 1; // Smaller time step

  // Boundary condition
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);

  fprintf (stderr, "Read in parameters!\n");
  run();
}

event init(i =  0)
{
  // TO DO
  // *> add an option to be able to restart the simu


  // Generate linealy spaced kx, ky according to specified # of modes, and
  //  interpolated F(kx,ky)
  T_Spectrum spectrum;
  spectrum = spectrum_gen_linear(N_mode, N_power, L, P, kp);
  geometric_beta (1./3., true); // Varying layer thickness
  // printf("kx[0] %f, ky[0] %f\n", spectrum.kx[0], spectrum.ky[0]); // OK
  // printf("f_kxky[0] %f\n", spectrum.F_kxky[0]); // OK but weird output -220458615281.473083
  // printf("eta %f, ux %f, uy %f, uz %f\n", ETA, UX, UY, UZ);
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
  fprintf (stderr,"Done initialization!\n");
  // sprintf(fileout, "%0*d.nc", pad, 0); // add the padding
  create_nc({zb, h, u, w}, file_out);
  write_nc(); 
  // Here add save to netcdf 
}

event end (t = tend) {
  fprintf (stderr,"tend = %f \n", tend);
  // create_nc({zb, h, u, w}, "end.nc");
  write_nc();
  // To Fix 31/10/25: This output is filled with nan. 

}
