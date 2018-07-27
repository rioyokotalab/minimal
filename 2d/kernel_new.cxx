#include "kernel.h"
using namespace exafmm;

int main(int argc, char ** argv) {
    P = atoi(argv[1]);

    Bodies jbodies(3);
    jbodies[0].X[0] = 6.5;
    jbodies[0].X[1] = 0;
    jbodies[0].q = 1;

    jbodies[1].X[0] = 6.0;
    jbodies[1].X[1] = 0;
    jbodies[1].q = 1;

    jbodies[2].X[0] = 5.5;
    jbodies[2].X[1] = 0;
    jbodies[2].q = 1;

    Cells cells(4);

    // P2M
    Cell *CJ = &cells[0];
    CJ->X[0] = 8;
    CJ->X[1] = 0;
    CJ->R = 2;
    CJ->BODY = &jbodies[0];
    CJ->NBODY = jbodies.size();
    CJ->M.resize(P, 0.0);
    P2M(CJ);

    //P2M CJ2
    Cell *CJ2 = &cells[1];
    CJ2->X[0] = 4;
    CJ2->X[1] = 0;
    CJ2->R = 2;
    CJ2->BODY = &jbodies[0];
    CJ2->NBODY = jbodies.size();
    CJ2->M.resize(P, 0.0);
    P2M(CJ2);

    // M2L
    Cell *CI = &cells[2];
    CI->X[0] = -8;
    CI->X[1] = 0;
    CI->R = 2;
    CI->L.resize(P, 0.0);

    //M2L (CI2,CJ) (CI,CJ2)
    Cell *CI2 = &cells[3];
    CI2->X[0] = -4;
    CI2->X[1] = 0;
    CI2->R = 2;
    CI2->L.resize(P, 0.0);

    M2L(CI, CJ);
    M2L(CI, CJ2);
    M2L(CI2, CJ);


    //L2P
    Bodies bodies(3);
    bodies[0].X[0] = - 6.5;
    bodies[0].X[1] = 0;
    bodies[0].q = 1;
    bodies[0].p = 0;

    bodies[1].X[0] = -6.0;
    bodies[1].X[1] = 0;
    bodies[1].q = 1;
    bodies[1].p = 0;

    bodies[2].X[0] = -5.5;
    bodies[2].X[1] = 0;
    bodies[2].q = 1;
    bodies[2].p = 0;

    for (int d = 0; d < 2; d++){
        bodies[0].F[d] = 0;
        bodies[1].F[d] = 0;
        bodies[2].F[d] = 0;
    }

    CI->BODY = &bodies[0];
    CI->NBODY = bodies.size();
    CI2->BODY = &bodies[0];
    CI2->NBODY = bodies.size();
    L2P(CI);
    L2P(CI2);

    // P2P
    P2P(CI2, CJ2);

    Bodies bodies2(3);
    for (size_t b = 0; b < bodies2.size(); b++) {
        bodies2[b] = bodies[b];
        bodies2[b].p = 0;
        for (int d = 0; d < 2; d++) bodies2[b].F[d] = 0;
   }

    CJ->NBODY = jbodies.size();
    CI->NBODY = bodies2.size();
    CI->BODY = &bodies2[0];
    P2PX(CI, CJ);

  // Verify results
  real_t pDif = 0, pNrm = 0, FDif = 0, FNrm = 0;
  for (size_t b=0; b<bodies.size(); b++) {
    pDif += (bodies[b].p - bodies2[b].p) * (bodies[b].p - bodies2[b].p);
    pNrm += bodies[b].p * bodies[b].p;
    FDif += (bodies[b].F[0] - bodies2[b].F[0]) * (bodies[b].F[0] - bodies2[b].F[0]) +
      (bodies[b].F[1] - bodies2[b].F[1]) * (bodies[b].F[1] - bodies2[b].F[1]);
    FNrm += bodies[b].F[0] * bodies[b].F[0] + bodies[b].F[1] * bodies[b].F[1];
  }

  printf("%-20s : %8.5e s\n","Rel. L2 Error (p)", sqrt(pDif/pNrm));
  printf("%-20s : %8.5e s\n","Rel. L2 Error (F)", sqrt(FDif/FNrm));
  return 0;
}
