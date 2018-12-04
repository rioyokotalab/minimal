#ifndef kernel_h
#define kernel_h
#include "exafmm.h"

namespace exafmm {
  //!< Weight of smoothing function
  inline real_t weight(Body * B, Cell * C) {
    assert(C->R > D);
    real_t x = fmin(C->R - std::abs(B->X[0] - C->X[0]), D);
    real_t y = fmin(C->R - std::abs(B->X[1] - C->X[1]), D);
    if (R0 - std::abs(B->X[0] - X0[0]) < D) x = D;
    if (R0 - std::abs(B->X[1] - X0[1]) < D) y = D;
    assert(x >= -D);
    assert(y >= -D);
    x /= D;
    y /= D;
    real_t w = (2 + 3 * x - x * x * x) / 4;
    w *= (2 + 3 * y - y * y * y) / 4;
    assert(0 <= w && w <= 1);
    return w;
  }

  //!< P2P kernel between cells Ci and Cj
  void P2PX(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->BODY;                                       // Target body pointer
    Body * Bj = Cj->BODY;                                       // Source body pointer
    for (int i=0; i<Ci->NBODY; i++) {                           // Loop over target bodies
      real_t p = 0;                              //  Initialize potential, force
      for (int j=0; j<Cj->NBODY; j++) {                         //  Loop over source bodies
        for (int d=0; d<2; d++) dX[d] = Bi[i].X[d] - Bj[j].X[d];//   Calculate distance vector
        real_t R2 = norm(dX);                                   //   Calculate distance squared
        if (R2 != 0) {                                          //   If not the same point
          p += 1;                                            //    Potential
        }                                                       //   End if for same point
      }                                                         //  End loop over source points
      Bi[i].p += p;                                             //  Accumulate potential
    }                                                           // End loop over target bodies
  }

  //!< P2P kernel between cells Ci and Cj
  void P2P(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->BODY;                                       // Target body pointer
    Body * Bj = Cj->BODY;                                       // Source body pointer
    for (int i=0; i<Ci->NBODY; i++) {                           // Loop over target bodies
      real_t p = 0;                              //  Initialize potential, force
      real_t wi = weight(&Bi[i], Ci);
      for (int j=0; j<Cj->NBODY; j++) {                         //  Loop over source bodies
        real_t wj = weight(&Bj[j], Cj);
        for (int d=0; d<2; d++) dX[d] = Bi[i].X[d] - Bj[j].X[d];//   Calculate distance vector
        real_t R2 = norm(dX);                                   //   Calculate distance squared
        if (R2 != 0) {                                          //   If not the same point
          p += wj;                                       //    Potential
          Bi[i].listp[Bj[j].I] += wj;
        }                                                       //   End if for same point
      }                                                         //  End loop over source points
      Bi[i].p += p * wi;                                        //  Accumulate potential
    }                                                           // End loop over target bodies
  }

  //!< P2M kernel for cell C
  void P2M(Cell * C) {
    for (Body * B=C->BODY; B!=C->BODY+C->NBODY; B++) {          // Loop over bodies
      for (int d=0; d<2; d++) dX[d] = B->X[d] - C->X[d];        //  Get distance vector
      real_t w = weight(B, C);
      complex_t Z(dX[0],dX[1]), powZ(1.0, 0.0);                 //  Convert to complex plane
      C->M[0] += w;                                      //  Add constant term
    }                                                           // End loop
  }

  //!< M2M kernel for one parent cell Ci
  void M2M(Cell * Ci) {
    for (Cell * Cj=Ci->CHILD; Cj!=Ci->CHILD+Ci->NCHILD; Cj++) { // Loop over child cells
      for (int d=0; d<2; d++) dX[d] = Cj->X[d] - Ci->X[d];      //  Get distance vector
      Ci->M[0] += Cj->M[0];                                   //   Add constant term
    }                                                           // End loop
  }

  //!< M2L kernel between cells Ci and Cj
  void M2L(Cell * Ci, Cell * Cj) {
    for (int d=0; d<2; d++) dX[d] = Ci->X[d] - Cj->X[d];        // Get distance vector
    complex_t Z(dX[0],dX[1]), powZn(1.0, 0.0), powZnk(1.0, 0.0), invZ(powZn/Z);// Convert to complex plane
    Ci->L[0] += Cj->M[0];                             // Log term (for 0th order)
  }

  //!< L2L kernel for one parent cell Cj
  void L2L(Cell * Cj) {
    for (Cell * Ci=Cj->CHILD; Ci<Cj->CHILD+Cj->NCHILD; Ci++) {  // Loop over child cells
      for (int d=0; d<2; d++) dX[d] = Ci->X[d] - Cj->X[d];      //  Get distance vector
      Ci->L[0] += Cj->L[0];                                   //   Add constant term
    }                                                           // End loop
  }

  //!< L2P kernel for cell C
  void L2P(Cell * C) {
    for (Body * B=C->BODY; B!=C->BODY+C->NBODY; B++) {          // Loop over bodies
      real_t w = weight(B, C);
      B->p += std::real(C->L[0]) * w;                           //  Add constant term
    }                                                           // End loop
  }
}

#endif
