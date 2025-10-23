/**
# Field scale wave breaking (multilayer solver)
*/

// #include "grid/multigrid.h"
// #include "layered/hydro.h"
// #include "layered/nh.h"
// #include "layered/remap.h"
// #include "layered/perfs.h"
// #include "input.h" // Used for older input method (matrices of quantities)

// Updated 2024/12/20: the newest Basilisk output_matrix is already compatible with MPI
// #include "output_mpi.h" // Antoon's function for MPI compatible matrix output
#include <bderembl/libs/extra.h>
//#include <spectrum.h>



double nu = 0.;
double dh[2] = {0., 1.};


int main(int argc,char* argv[]) {
  params = array_new();
  add_param ("N", &N, "int");
  add_param ("nu", &nu, "double");
  add_param ("dh", &dh[0], "array");

  // Search for the configuration file with a given path or read params.in
  if (argc == 2)
    strcpy(file_param,argv[1]); // default: params.in
  // printf("nu before read %f\n", nu);
  read_params(file_param);
  // printf("nu after read %f\n", nu);

  create_outdir(); // Create a directory 'outdir_000X'
  backup_config(file_param); 
}
