#include "build_tree.h"
#include "kernel0.h"
#include "timer.h"
#include "traverse.h"
#include <iostream>
using namespace exafmm;

int main(int argc, char ** argv) {
  const int numBodies = 3;                                  // Number of bodies
  P = 30;                                                       // Order of expansions
  D = 0.76;                                                     // Buffer size
  ncrit = 2;                                                   // Number of bodies per leaf cell
  theta = 0.0;                                                  // Multipole acceptance criterion

  printf("--- %-16s ------------\n", "FMM Profiling");          // Start profiling
  //! Initialize bodie
  start("Initialize bodies");                                   // Start timer
  Bodies bodies(numBodies);                                     // Initialize bodies
  real_t average = 0;                                           // Average charge
  srand48(0);                                                   // Set seed for random number generator
  for (size_t b=0; b<bodies.size(); b++) {                      // Loop over bodies
    bodies[b].I = b;                                            //  Body index
    bodies[b].listp.resize(numBodies);
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
  stop("Initialize bodies");                                    // Stop timer

  //! Build tree
  start("Build tree");                                          // Start timer
  Bodies bodies2 = bodies;
  Cells cells = buildTree(bodies);                              // Build tree
  stop("Build tree");                                           // Stop timer
  for (size_t b=0; b<bodies.size(); b++) {                      // Loop over bodies
    bodies[b].listp.resize(numBodies);
  }
  //! FMM evaluation
  start("P2M & M2M");                                           // Start timer
  upwardPass(cells);                                            // Upward pass for P2M, M2M
  stop("P2M & M2M");                                            // Stop timer
  start("M2L & P2P");                                           // Start timer
  horizontalPass(cells, cells);                                 // Horizontal pass for M2L, P2P
  stop("M2L & P2P");                                            // Stop timer
  start("L2L & L2P");                                           // Start timer
  downwardPass(cells);                                          // Downward pass for L2L, L2P
  stop("L2L & L2P");                                            // Stop timer
  Bodies jbodies = bodies2;
  joinBuffer(cells, jbodies);
  bodies = jbodies;

  //! Direct N-Body
  start("Direct N-Body");                                       // Start timer
  const int numTargets = std::min(100,int(bodies.size()));           // Number of targets for checking answer
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
  direct(bodies, jbodies);                                      // Direct N-Body
  stop("Direct N-Body");                                        // Stop timer

  //! Verify result
  double pDif = 0, pNrm = 0;
  for (size_t b=0; b<bodies2.size(); b++) {                      // Loop over bodies & bodies2
    pDif += (bodies[b].p - bodies2[b].p) * (bodies[b].p - bodies2[b].p);// Difference of potential
    //if (pDif > 1e-8) std::cout << bodies[b].I << " " << bodies[b].p << " " << bodies2[b].p << std::endl;
    pNrm += bodies[b].p * bodies[b].p;                           //  Value of potential
  }                                                             // End loop over bodies & bodies2
  printf("--- %-16s ------------\n", "FMM vs. direct");         // Print message
  printf("%-20s : %8.5e s\n","Rel. L2 Error (p)", sqrt(pDif/pNrm));// Print potential error
  return 0;
}
