#include "build_tree.h"
#include "kernel.h"
#include "timer.h"
#include "string.h"
#include <algorithm>
#include <iostream>
//#include "traverse_eager.h"
#include "traverse_lazy.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  P = atoi(argv[1]);                                            // Order of expansions
  ncrit = 64;                                                   // Number of bodies per leaf cell
  theta = 0.4;                                                  // Multipole acceptance criterion
  int err = 0, N = 1 << atoi(argv[2]);
  printf("--- %-16s ------------\n", "FMM Profiling");          // Start profiling
  start("Initialize bodies");                                   // Start timer
  Bodies bodies(N);
  FILE * file = fopen("../../sandbox/direct/xi.dat","r");
  for (int i=0; i<N; i++) {
    float dummy;
    err = fscanf(file,"%f ", &dummy);
    err = fscanf(file,"%f ", &dummy);
    err = fscanf(file,"%f ", &dummy);
    err = fscanf(file,"%f ", &bodies[i].X[0]);
    err = fscanf(file,"%f ", &bodies[i].X[1]);
    err = fscanf(file,"%f ", &bodies[i].X[2]);
    err = fscanf(file,"%f ", &bodies[i].q);
    bodies[i].p = 0;
    for (int d=0; d<3; d++) bodies[i].F[d] = 0;
    bodies[i].i = i;
  }
  Bodies bodies2 = bodies;
  fclose(file);
  stop("Initialize bodies");                                    // Stop timer

  //! Build tree
  start("Build tree");                                          // Start timer
  Cells cells = buildTree(bodies);                              // Build tree
  stop("Build tree");                                           // Stop timer

  //! FMM evaluation
  start("P2M & M2M");                                           // Start timer
  initKernel();                                                 // Initialize kernel
  upwardPass(cells);                                            // Upward pass for P2M, M2M
  stop("P2M & M2M");                                            // Stop timer
  start("M2L & P2P");                                           // Start timer
  horizontalPass(cells, cells);                                 // Horizontal pass for M2L, P2P
  stop("M2L & P2P");                                            // Stop timer
  start("L2L & L2P");                                           // Start timer
  downwardPass(cells);                                          // Downward pass for L2L, L2P
  stop("L2L & L2P");                                            // Stop timer
  std::sort(bodies.begin(), bodies.end(), compare);

  //! Verify result
  char filename[128];
  strcpy(filename,"../../sandbox/direct/xi");
  strcat(filename,argv[2]);
  strcat(filename,".dat");
  file = fopen(filename,"r");
  double pDif = 0, pNrm = 0, FDif = 0, FNrm = 0;
  for (size_t b=0; b<bodies.size(); b++) {
    double p, ax, ay, az;
    err = fscanf(file,"%le ", &p);
    err = fscanf(file,"%le ", &ax);
    err = fscanf(file,"%le ", &ay);
    err = fscanf(file,"%le ", &az);
    p /= 16;
    ax /= 16*16*16;
    ay /= 16*16*16;
    az /= 16*16*16;
    pDif += (bodies[b].p - p) * (bodies[b].p - p);
    pNrm += p * p;
    FDif += (bodies[b].F[0] - ax) * (bodies[b].F[0] - ax) +
      (bodies[b].F[1] - ay) * (bodies[b].F[1] - ay) +
      (bodies[b].F[2] - az) * (bodies[b].F[2] - az);
    FNrm += ax * ax + ay * ay + az * az;
  }
  fclose(file);
  printf("--- %-16s ------------\n", "FMM vs. direct");         // Print message
  printf("%-20s : %8.5e s\n","Rel. L2 Error (p)", sqrt(pDif/pNrm));
  printf("%-20s : %8.5e s\n","Rel. L2 Error (F)", sqrt(FDif/FNrm));
  return 0*err;
}
