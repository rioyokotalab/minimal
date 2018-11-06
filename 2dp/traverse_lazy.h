#ifndef traverse_lazy_h
#define traverse_lazy_h
#include "types.h"

namespace exafmm {

  //! Recursive call to post-order tree traversal for upward pass
  void upwardPass(Cell * Ci) {
    for (Cell * Cj=Ci->CHILD; Cj!=Ci->CHILD+Ci->NCHILD; Cj++) { // Loop over child cells
#pragma omp task untied if(Cj->NBODY > 100)                     //  Start OpenMP task if large enough task
      upwardPass(Cj);                                           //  Recursive call for child cell
    }                                                           // End loop over child cells
#pragma omp taskwait                                            // Synchronize OpenMP tasks
    Ci->M.resize(P, 0.0);                                       // Allocate and initialize multipole coefs
    Ci->L.resize(P, 0.0);                                       // Allocate and initialize local coefs
    if (Ci->NCHILD == 0) P2M(Ci);                               // P2M kernel
    M2M(Ci);                                                    // M2M kernel
  }

  //! Upward pass interface
  void upwardPass(Cells & cells) {
#pragma omp parallel                                            // Start OpenMP
#pragma omp single nowait                                       // Start OpenMP single region with nowait
    upwardPass(&cells[0]);                                      // Pass root cell to recursive call
  }

  //! Recursive call to dual tree traversal for list construction
  void getList(Cell * Ci, Cell * Cj) {
    for (int d=0; d<2; d++) dX[d] = Ci->X[d] - Cj->X[d] - Xperiodic[d];// Distance vector from source to target
    real_t R2 = norm(dX) * theta * theta;                       // Scalar distance squared
    if (R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R)) {               // If distance is far enough
      Ci->listM2L.push_back(Cj);                                //  Add to M2L list
    } else if (Ci->NCHILD == 0 && Cj->NCHILD == 0) {            // Else if both cells are leafs
      Ci->listP2P.push_back(Cj);                                //  Add to P2P list
    } else if (Cj->NCHILD == 0 || (Ci->R >= Cj->R && Ci->NCHILD != 0)) {// Else if Cj is leaf or Ci is larger
      for (Cell * ci=Ci->CHILD; ci!=Ci->CHILD+Ci->NCHILD; ci++) {// Loop over Ci's children
        getList(ci, Cj);                                        //   Recursive call to target child cells
      }                                                         //  End loop over Ci's children
    } else {                                                    // Else if Ci is leaf or Cj is larger
      for (Cell * cj=Cj->CHILD; cj!=Cj->CHILD+Cj->NCHILD; cj++) {//  Loop over Cj's children
        getList(Ci, cj);                                        //   Recursive call to source child cells
      }                                                         //  End loop over Cj's children
    }                                                           // End if for leafs and Ci Cj size
  }

  //! Evaluate M2L, P2P kernels
  void evaluate(Cells & cells) {
#pragma omp parallel for schedule(dynamic)
    for (size_t i=0; i<cells.size(); i++) {                     // Loop over cells
      for (size_t j=0; j<cells[i].listM2L.size(); j++) {        //  Loop over M2L list
        M2L(&cells[i],cells[i].listM2L[j]);                     //   M2L kernel
      }                                                         //  End loop over M2L list
      for (size_t j=0; j<cells[i].listP2P.size(); j++) {        //  Loop over P2P list
        P2P(&cells[i],cells[i].listP2P[j]);                     //   P2P kernel
      }                                                         //  End loop over P2P list
    }                                                           // End loop over cells
  }

  //! Horizontal pass for periodic images
  void periodic(Cell * Ci0, Cell * Cj0, real_t cycle) {
    Cells pcells(9);                                            // Create cells
    for (size_t c=0; c<pcells.size(); c++) {                    // Loop over periodic cells
      pcells[c].M.resize(P, 0.0);                               //  Allocate & initialize M coefs
      pcells[c].L.resize(P, 0.0);                               //  Allocate & initialize L coefs
    }                                                           // End loop over periodic cells
    Cell * Ci = &pcells.back();                                 // Last cell is periodic parent cell
    *Ci = *Cj0;                                                 // Copy values from source root
    Ci->NCHILD = 8;                                             // Number of child cells for periodic center cell
    for (int level=0; level<images-1; level++) {                // Loop over sublevels of tree
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
        for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
          if (ix != 0 || iy != 0) {                             //    If periodic cell is not at center
            for (int cx=-1; cx<=1; cx++) {                      //     Loop over x periodic direction (child)
              for (int cy=-1; cy<=1; cy++) {                    //      Loop over y periodic direction (child)
                Xperiodic[0] = (ix * 3 + cx) * cycle;           //       Coordinate offset for x periodic direction
                Xperiodic[1] = (iy * 3 + cy) * cycle;           //       Coordinate offset for y periodic direction
                M2L(Ci0, Ci);                                   //       Perform M2L kernel
              }                                                 //      End loop over y periodic direction (child)
            }                                                   //     End loop over x periodic direction (child)
          }                                                     //    Endif for periodic center cell
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      Cell * Cj = &pcells.front();                              //  Iterator of periodic neighbor cells
      Ci->CHILD = Cj;                                           //  Child cells for periodic center cell
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
        for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
          if( ix != 0 || iy != 0) {                             //    If periodic cell is not at center
            Cj->X[0] = Ci->X[0] + ix * cycle;                   //     Set new x coordinate for periodic image
            Cj->X[1] = Ci->X[1] + iy * cycle;                   //     Set new y cooridnate for periodic image
            for (int n=0; n<P; n++) Cj->M[n] = Ci->M[n];        //     Copy multipoles to new periodic image
            Cj++;                                               //     Increment periodic cell iterator
          }                                                     //    Endif for periodic center cell
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      M2M(Ci);                                                  //  Evaluate periodic M2M kernels for this sublevel
      cycle *= 3;                                               //  Increase center cell size three times
    }                                                           // End loop over sublevels of tree
  }

  //! Horizontal pass interface
  void horizontalPass(Cells & icells, Cells & jcells, real_t cycle) {
    if (images == 0) {                                          // If non-periodic boundary condition
      for (int d=0; d<2; d++) Xperiodic[d] = 0;                 //  No periodic shift
      getList(&icells[0], &jcells[0]);                          //  Pass root cell to recursive call
      evaluate(icells);                                         //  Evaluate M2L & P2P kernels
    } else {                                                    // If periodic boundary condition
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
        for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
          Xperiodic[0] = ix * cycle;                            //    Coordinate shift for x periodic direction
          Xperiodic[1] = iy * cycle;                            //    Coordinate shift for y periodic direction
          getList(&icells[0], &jcells[0]);                      //    Pass root cell to recursive call
          evaluate(icells);                                     //    Evaluate M2L & P2P kernels
          for (size_t i=0; i<icells.size(); i++) {              //    Loop over target cells
            icells[i].listM2L.clear();                          //     Clear M2L interaction list
            icells[i].listP2P.clear();                          //     Clear P2P interaction list
          }                                                     //    End loop over target cells
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      periodic(&icells[0], &jcells[0], cycle);                  //  Horizontal pass for periodic images
    }                                                           // End if for periodic boundary condition
  }                                                             // End if for empty cell vectors

  //! Recursive call to pre-order tree traversal for downward pass
  void downwardPass(Cell * Cj) {
    L2L(Cj);                                                    // L2L kernel
    if (Cj->NCHILD == 0) L2P(Cj);                               // L2P kernel
    for (Cell * Ci=Cj->CHILD; Ci!=Cj->CHILD+Cj->NCHILD; Ci++) { // Loop over child cells
#pragma omp task untied if(Ci->NBODY > 100)                     //  Start OpenMP task if large enough task
      downwardPass(Ci);                                         //  Recursive call for child cell
    }                                                           // End loop over child cells
#pragma omp taskwait                                            // Synchronize OpenMP tasks
  }

  //! Downward pass interface
  void downwardPass(Cells & cells) {
#pragma omp parallel                                            // Start OpenMP
#pragma omp single nowait                                       // Start OpenMP single region with nowait
    downwardPass(&cells[0]);                                    // Pass root cell to recursive call
  }

  //! Direct summation
  void direct(Bodies & bodies, Bodies & jbodies, real_t cycle) {
    Cells cells(2);                                             // Define a pair of cells to pass to P2P kernel
    Cell * Ci = &cells[0];                                      // Allocate single target cell
    Cell * Cj = &cells[1];                                      // Allocate single source cell
    Ci->BODY = &bodies[0];                                      // Pointer of first target body
    Ci->NBODY = bodies.size();                                  // Number of target bodies
    Cj->BODY = &jbodies[0];                                     // Pointer of first source body
    Cj->NBODY = jbodies.size();                                 // Number of source bodies
    int prange = 0;                                             // Range of periodic images
    for (int i=0; i<images; i++) {                              // Loop over periodic image sublevels
      prange += int(powf(3.,i));                                //  Accumulate range of periodic images
    }                                                           // End loop over perioidc image sublevels
#pragma omp parallel for collapse(2)
    for (int ix=-prange; ix<=prange; ix++) {                    // Loop over x periodic direction
      for (int iy=-prange; iy<=prange; iy++) {                  //  Loop over y periodic direction
        Xperiodic[0] = ix * cycle;                              //   Coordinate shift for x periodic direction
        Xperiodic[1] = iy * cycle;                              //   Coordinate shift for y periodic direction
        P2PX(Ci, Cj);                                           //   Evaluate P2P kernel
      }                                                         //  End loop over y periodic direction
    }                                                           // End loop over x periodic direction
  }
}

#endif
