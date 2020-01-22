#include "kernel2.h"
#include "build_tree.h"
#include "timer.h"
//#include "traverse_eager.h"
#include "traverse_lazy.h"
#include <iostream>
using namespace exafmm;

int main(int argc, char ** argv) {
  const int numBodies = atoi(argv[1]);                          // Number of bodies
  const real_t cycle = 2 * M_PI;                                // Cycle of periodic boundary condition
  P = atoi(argv[2]);                                            // Order of expansions
  D = atof(argv[3]);                                            // Buffer size
  ncrit = 64;                                                   // Number of bodies per leaf cell
  theta = 0.5;                                                  // Multipole acceptance criterion
  images = 0;                                                   // 3^images * 3^images * 3^images periodic images

  //! Initialize bodies
  Bodies bodies(numBodies);                                     // Initialize bodies
  real_t average = 0;                                           // Average charge
  srand48(0);                                                   // Set seed for random number generator
  for (size_t b=0; b<bodies.size(); b++) {                      // Loop over bodies
    bodies[b].I = b;                                            //  Body index
    for (int d=0; d<2; d++) {                                   //  Loop over dimension
      bodies[b].X[d] = drand48() * 2 * M_PI - M_PI;             //   Initialize positions
    }                                                           //  End loop over dimension
    bodies[b].q = drand48() - .5;                               //  Initialize charge
    average += bodies[b].q;                                     //  Accumulate charge
    bodies[b].p = 0;                                            //  Clear potential
    for (int d=0; d<2; d++) bodies[b].F[d] = 0;                 //  Clear force
  }                                                             // End loop over bodies
  average /= bodies.size();                                     // Average charge
  for (size_t b=0; b<bodies.size(); b++) {                      // Loop over bodies
    bodies[b].q -= average;                                     // Charge neutral
  }                                                             // End loop over bodies

  //! Build tree
  Bodies bodies2 = bodies;
  Cells  cells = buildTree(bodies,cycle);                       // Build tree

  //! FMM evaluation
  start("FMM time");                                           // Start timer
  upwardPass(cells);                                            // Upward pass for P2M, M2M
  horizontalPass(cells, cells, cycle);                          // Horizontal pass for M2L, P2P
  downwardPass(cells);                                          // Downward pass for L2L, L2P
  stop("FMM time");                                           // Start timer
  Bodies jbodies = bodies2;
  joinBuffer(cells, jbodies, cycle);
  bodies = jbodies;

  // Direct N-Body
  const int numTargets = 10;                                    // Number of targets for checking answer
  int stride = bodies.size() / numTargets;                      // Stride of sampling
  for (int b=0; b<numTargets; b++) {                            // Loop over target samples
    bodies[b] = bodies[b*stride];                               //  Sample targets
  }                                                             // End loop over target samples
  bodies.resize(numTargets);                                    // Resize bodies
  jbodies = bodies2;                                            // Save bodies in jbodies
  bodies2 = bodies;                                             // Backup bodies
  for (size_t b=0; b<bodies.size(); b++) {                      // Loop over bodies
    bodies[b].p = 0;                                            //  Clear potential
    for (int d=0; d<2; d++) bodies[b].F[d] = 0;                 //  Clear force
  }                                                             // End loop over bodies
  direct(bodies, jbodies, cycle);                               // Direct N-Body

  //! Verify result
  double pDif = 0, pNrm = 0, FDif = 0, FNrm = 0;
  for (size_t b=0; b<bodies2.size(); b++) {                      // Loop over bodies & bodies2
    pDif += (bodies[b].p - bodies2[b].p) * (bodies[b].p - bodies2[b].p);// Difference of potential
    pNrm += bodies[b].p * bodies[b].p;                        //  Value of potential
    FDif += (bodies[b].F[0] - bodies2[b].F[0]) * (bodies[b].F[0] - bodies2[b].F[0])// Difference of force
      + (bodies[b].F[0] - bodies2[b].F[0]) * (bodies[b].F[0] - bodies2[b].F[0]);// Difference of force
    FNrm += bodies[b].F[0] * bodies[b].F[0] + bodies[b].F[1] * bodies[b].F[1];//  Value of force
  }                                                             // End loop over bodies & bodies2
  printf("%-20s : %8.5e\n","Rel. L2 Error (F)", sqrt(FDif/FNrm));// Print force error
  return 0;
}
