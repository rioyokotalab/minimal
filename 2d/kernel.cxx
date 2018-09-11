#include "kernel2.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  P = atoi(argv[1]);
  D = 0.25;

  Bodies jbodies(9);
  jbodies[0].X[0] = 8.5;
  jbodies[0].X[1] = 1.5;
  jbodies[0].q = 1;

  jbodies[1].X[0] = 8.0;
  jbodies[1].X[1] = 1.5;
  jbodies[1].q = 1;

  jbodies[2].X[0] = 7.5;
  jbodies[2].X[1] = 1.5;
  jbodies[2].q = 1;

  jbodies[3].X[0] = 8.5;
  jbodies[3].X[1] = 2.0;
  jbodies[3].q = 1;

  jbodies[4].X[0] = 8.0;
  jbodies[4].X[1] = 2.0;
  jbodies[4].q = 1;

  jbodies[5].X[0] = 7.5;
  jbodies[5].X[1] = 2.0;
  jbodies[5].q = 1;

  jbodies[6].X[0] = 8.5;
  jbodies[6].X[1] = 2.5;
  jbodies[6].q = 1;

  jbodies[7].X[0] = 8.0;
  jbodies[7].X[1] = 2.5;
  jbodies[7].q = 1;

  jbodies[8].X[0] = 7.5;
  jbodies[8].X[1] = 2.5;
  jbodies[8].q = 1;

  Cells cells(10);

  // P2M
  Cell *CJ = &cells[0];
  CJ->X[0] = 10;
  CJ->X[1] = 0;
  CJ->R = 2;
  CJ->BODY = &jbodies[0];
  CJ->NBODY = jbodies.size();
  CJ->M.resize(P, 0.0);
  P2M(CJ);

  Cell *CJ2 = &cells[1];
  CJ2->X[0] = 6;
  CJ2->X[1] = 0;
  CJ2->R = 2;
  CJ2->BODY = &jbodies[0];
  CJ2->NBODY = jbodies.size();
  CJ2->M.resize(P, 0.0);
  P2M(CJ2);

  Cell *CJ3 = &cells[2];
  CJ3->X[0] = 10;
  CJ3->X[1] = 4;
  CJ3->R = 2;
  CJ3->BODY = &jbodies[0];
  CJ3->NBODY = jbodies.size();
  CJ3->M.resize(P, 0.0);
  P2M(CJ3);

  Cell *CJ4 = &cells[3];
  CJ4->X[0] = 6;
  CJ4->X[1] = 4;
  CJ4->R = 2;
  CJ4->BODY = &jbodies[0];
  CJ4->NBODY = jbodies.size();
  CJ4->M.resize(P, 0.0);
  P2M(CJ4);

  //M2M
  Cell *CJ5 =&cells[4];
  CJ5->CHILD = CJ;
  CJ5->NCHILD = 4;

  CJ5->X[0] = 8;
  CJ5->X[1] = 2;
  CJ5->R = 4;
  CJ5->M.resize(P,0.0);
  M2M(CJ5);


  // M2L
  Cell *CI = &cells[5];
  CI->X[0] = -10;
  CI->X[1] = 0;
  CI->R = 2;
  CI->L.resize(P, 0.0);

  Cell *CI2 = &cells[6];
  CI2->X[0] = -6;
  CI2->X[1] = 0;
  CI2->R = 2;
  CI2->L.resize(P, 0.0);

  Cell *CI3 = &cells[7];
  CI3->X[0] = -10;
  CI3->X[1] = 4;
  CI3->R = 2;
  CI3->L.resize(P, 0.0);

  Cell *CI4 = &cells[8];
  CI4->X[0] = -6;
  CI4->X[1] = 4;
  CI4->R = 2;
  CI4->L.resize(P, 0.0);

  Cell *CI5 = &cells[9];
  CI5->CHILD = CI;
  CI5->NCHILD = 4;
  CI5->X[0] = -8;
  CI5->X[1] = 2;
  CI5->L.resize(P,0.0);

  M2L(CI5,CJ5);

  // L2L
  L2L(CI5);

  //L2P
  Bodies bodies(9);
  bodies[0].X[0] = -8.5;
  bodies[0].X[1] = 1.5;
  bodies[0].q = 1;

  bodies[1].X[0] = -8.0;
  bodies[1].X[1] = 1.5;
  bodies[1].q = 1;

  bodies[2].X[0] = -7.5;
  bodies[2].X[1] = 1.5;
  bodies[2].q = 1;

  bodies[3].X[0] = -8.5;
  bodies[3].X[1] = 2.0;
  bodies[3].q = 1;

  bodies[4].X[0] = -8.0;
  bodies[4].X[1] = 2.0;
  bodies[4].q = 1;

  bodies[5].X[0] = -7.5;
  bodies[5].X[1] = 2.0;
  bodies[5].q = 1;

  bodies[6].X[0] = -8.5;
  bodies[6].X[1] = 2.5;
  bodies[6].q = 1;

  bodies[7].X[0] = -8.0;
  bodies[7].X[1] = 2.5;
  bodies[7].q = 1;

  bodies[8].X[0] = -7.5;
  bodies[8].X[1] = 2.5;
  bodies[8].q = 1;

  for (size_t b = 0; b < bodies.size(); b++) {
    bodies[b].p = 0;
    for (int d=0; d<2; d++){
      bodies[b].F[d] = 0;
    }
  }

  CI->BODY = &bodies[0];
  CI->NBODY = bodies.size();
  CI2->BODY = &bodies[0];
  CI2->NBODY = bodies.size();
  CI3->BODY = &bodies[0];
  CI3->NBODY = bodies.size();
  CI4->BODY = &bodies[0];
  CI4->NBODY = bodies.size();
  L2P(CI);
  L2P(CI2);
  L2P(CI3);
  L2P(CI4);

  // P2P
#if 0
  P2P(CI, CJ);
  P2P(CI, CJ2);
  P2P(CI, CJ3);
  P2P(CI, CJ4);
  P2P(CI2, CJ);
  P2P(CI2, CJ2);
  P2P(CI2, CJ3);
  P2P(CI2, CJ4);
  P2P(CI3, CJ);
  P2P(CI3, CJ2);
  P2P(CI3, CJ3);
  P2P(CI3, CJ4);
  P2P(CI4, CJ);
  P2P(CI4, CJ2);
  P2P(CI4, CJ3);
  P2P(CI4, CJ4);
#endif

  Bodies bodies2(bodies.size());
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
