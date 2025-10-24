/**
# Field scale wave breaking (multilayer solver)
*/

// #include "grid/multigrid.h"
// #include "layered/hydro.h"
// #include "layered/nh.h"
// #include "layered/remap.h"
// #include "layered/perfs.h"

#include <bderembl/libs/extra.h>
#include <bderembl/libs/netcdf_bas.h>
#define g_ 9.81
#include "spectrum.h"

/*
DEFAULT PARAMETERS
*/
char namlist[80] = "namlist.toml"; // file name of namlist 
// -> Initial conditions
double P = 0.2;             // energy level (estimated so that kpHs is reasonable)
int coeff_kpL0 = 10;        // kpL0 = coeff_kpL0 * pi
int N_mode = 32;            // Number of modes in wavenumber space
int N_power = 5;            // directional spreading coeff
int F_shape = 0;    // shape of the initial spectrum
// -> Domain definition
int N_grid = 1024;         // number of x and y gridpoints
double L = 200.0;           // domain size
int N_layer = 2;           // number of layers

int main(int argc,char* argv[]) {
  params = array_new();
  add_param ("P", &P, "double");
  add_param ("coeff_kpL0", &coeff_kpL0, "int");
  add_param ("N_mode", &N_mode, "int");
  add_param ("N_power", &N_power, "int");
  add_param ("F_shape", &F_shape, "int");
  add_param ("N_grid", &N_grid, "int");
  add_param ("L", &N_grid, "double");
  add_param ("N_layer", &N_layer, "int");
  double kp = PI * coeff_kpL0 / L; // kpL=coeff x pi peak wavelength

  // Search for the configuration file with a given path or read params.in
  if (argc == 2)
    strcpy(file_param,argv[1]); 
  else 
    strcpy(file_param, namlist);

  read_params(file_param);
  // create_outdir(); // Create a directory 'outdir_000X'
  // backup_config(file_param); 

  

  T_Spectrum spectrum;
  // Generate linealy spaced kx, ky according to specified # of modes, and 
  //  interpolated F(kx,ky)
  spectrum = spectrum_gen_linear(N_mode, N_power, L, P, kp);


}
