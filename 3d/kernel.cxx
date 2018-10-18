#include "kernel.h"
#include "string.h"
#include <iostream>
using namespace exafmm;

int main(int argc, char ** argv) {
  P = atoi(argv[1]);
  initKernel();
  int err = 0, N = 1 << atoi(argv[2]);
  Bodies bodies(N);
  Bodies jbodies(N);
  FILE * file = fopen("../../sandbox/direct/xj.dat","r");
  for (int i=0; i<N; i++) {
    err = fscanf(file,"%f ", &bodies[i].X[0]);
    err = fscanf(file,"%f ", &bodies[i].X[1]);
    err = fscanf(file,"%f ", &bodies[i].X[2]);
    err = fscanf(file,"%f ", &jbodies[i].X[0]);
    err = fscanf(file,"%f ", &jbodies[i].X[1]);
    err = fscanf(file,"%f ", &jbodies[i].X[2]);
    err = fscanf(file,"%f ", &jbodies[i].q);
    bodies[i].p = 0;
    for (int d=0; d<3; d++) bodies[i].F[d] = 0;
  }
  fclose(file);
  Cells cells(2);
  // P2M
  Cell * Cj = &cells[0];
  for (int d=0; d<3; d++) Cj->X[d] = 5.5;
  Cj->R = 0.5;
  Cj->BODY = &jbodies[0];
  Cj->NBODY = jbodies.size();
  Cj->M.resize(NTERM, 0.0);
  P2M(Cj);

  // M2L
  Cell * Ci = &cells[1];
  for (int d=0; d<3; d++) Ci->X[d] = 0.5;
  Ci->R = 0.5;
  Ci->BODY = &bodies[0];
  Ci->NBODY = bodies.size();
  Ci->L.resize(NTERM, 0.0);
  M2L(Ci, Cj);

  // L2P
  L2P(Ci);

  // Verify results
  char filename[128];
  strcpy(filename,"../../sandbox/direct/xj");
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
  printf("%-20s : %8.5e s\n","Rel. L2 Error (p)", sqrt(pDif/pNrm));
  printf("%-20s : %8.5e s\n","Rel. L2 Error (F)", sqrt(FDif/FNrm));
  return 0*err;
}
