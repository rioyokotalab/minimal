#ifndef buildtree_h
#define buildtree_h
#include "exafmm.h"
#include <iostream>

namespace exafmm {
  //! Get bounding box of bodies
  void getBounds(Bodies & bodies) {
    real_t Xmin[2], Xmax[2];                                    // Min, max of domain
    for (int d=0; d<2; d++) Xmin[d] = Xmax[d] = bodies[0].X[d]; // Initialize Xmin, Xmax
    for (size_t b=0; b<bodies.size(); b++) {                    // Loop over range of bodies
      for (int d=0; d<2; d++) Xmin[d] = fmin(bodies[b].X[d], Xmin[d]);//  Update Xmin
      for (int d=0; d<2; d++) Xmax[d] = fmax(bodies[b].X[d], Xmax[d]);//  Update Xmax
    }                                                           // End loop over range of bodies
    for (int d=0; d<2; d++) X0[d] = (Xmax[d] + Xmin[d]) / 2;    // Calculate center of domain
    R0 = 0;                                                     // Initialize localRadius
    for (int d=0; d<2; d++) {                                   // Loop over dimensions
      R0 = fmax(X0[d] - Xmin[d], R0);                           //  Calculate min distance from center
      R0 = fmax(Xmax[d] - X0[d], R0);                           //  Calculate max distance from center
    }                                                           // End loop over dimensions
    R0 *= 1.00001;                                              // Add some leeway to radius
  }

  //! Build cells of tree adaptively using a top-down approach based on recursion
  void buildCells(Body * bodies, Body * buffer, int begin, int end, Cell * cell, Cells & cells,
                  real_t * X, real_t R, int level=0, bool direction=false) {
    //! Create a tree cell
    cell->BODY = bodies + begin;                                // Pointer of first body in cell
    if(direction) cell->BODY = buffer + begin;                  // Pointer of first body in cell
    cell->NBODY = end - begin;                                  // Number of bodies in cell
    cell->NCHILD = 0;                                           // Initialize counter for child cells
    for (int d=0; d<2; d++) cell->X[d] = X[d];                  // Center position of cell
    cell->R = R / (1 << level);                                 // Cell radius
    maxlevel = std::max(maxlevel, level);
    //! If cell is a leaf
    if (end - begin <= ncrit) {                                 // If number of bodies is less than threshold
      if (direction) {                                          //  If direction of data is from bodies to buffer
        for (int i=begin; i<end; i++) {                         //   Loop over bodies in cell
          for (int d=0; d<2; d++) buffer[i].X[d] = bodies[i].X[d];//  Copy bodies coordinates to buffer
          buffer[i].q = bodies[i].q;                            //    Copy bodies source to buffer
        }                                                       //   End loop over bodies in cell
      }                                                         //  End if for direction of data
      return;                                                   //  Return without recursion
    }                                                           // End if for number of bodies
    //! Count number of bodies in each quadrant
    int size[4] = {0,0,0,0};
    real_t x[2];                                                // Coordinates of bodies
    for (int i=begin; i<end; i++) {                             // Loop over bodies in cell
      for (int d=0; d<2; d++) x[d] = bodies[i].X[d];            //  Position of body
      int quadrant = (x[0] > X[0]) + ((x[1] > X[1]) << 1);      //  Which quadrant body belongs to
      size[quadrant]++;                                         //  Increment body count in quadrant
    }                                                           // End loop over bodies in cell
    //! Exclusive scan to get offsets
    int offset = begin;                                         // Offset of first quadrant
    int offsets[4], counter[4];                                 // Offsets and counter for each quadrant
    for (int i=0; i<4; i++) {                                   // Loop over elements
      offsets[i] = offset;                                      //  Set value
      offset += size[i];                                        //  Increment offset
      if (size[i]) cell->NCHILD++;                              //  Increment child cell counter
    }                                                           // End loop over elements
    //! Sort bodies by quadrant
    for (int i=0; i<4; i++) counter[i] = offsets[i];            // Copy offsets to counter
    for (int i=begin; i<end; i++) {                             // Loop over bodies
      for (int d=0; d<2; d++) x[d] = bodies[i].X[d];            //  Position of body
      int quadrant = (x[0] > X[0]) + ((x[1] > X[1]) << 1);      //  Which quadrant body belongs to`
      for (int d=0; d<2; d++) buffer[counter[quadrant]].X[d] = bodies[i].X[d];// Permute bodies coordinates out-of-place according to quadrant
      buffer[counter[quadrant]].q = bodies[i].q;                //  Permute bodies sources out-of-place according to quadrant
      counter[quadrant]++;                                      //  Increment body count in quadrant
    }                                                           // End loop over bodies
    //! Loop over children and recurse
    real_t Xchild[2];                                           // Coordinates of children
    cells.resize(cells.size()+cell->NCHILD);                    // Resize cell vector
    Cell * child = &cells.back() - cell->NCHILD + 1;            // Pointer for first child cell
    cell->CHILD = child;                                        // Point to first child cell
    int c = 0;                                                  // Counter for child cells
    for (int i=0; i<4; i++) {                                   // Loop over children
      for (int d=0; d<2; d++) Xchild[d] = X[d];                 //  Initialize center position of child cell
      real_t r = R / (1 << (level + 1));                        //  Radius of cells for child's level
      for (int d=0; d<2; d++) {                                 //  Loop over dimensions
        Xchild[d] += r * (((i & 1 << d) >> d) * 2 - 1);         //   Shift center position to that of child cell
      }                                                         //  End loop over dimensions
      if (size[i]) {                                            //  If child exists
        buildCells(buffer, bodies, offsets[i], offsets[i] + size[i],// Recursive call for each child
                   &child[c], cells, Xchild, R, level+1, !direction);
        c++;                                                    //   Increment child cell counter
      }                                                         //  End if for child
    }                                                           // End loop over children
  }

  //! Recursive call to dual tree traversal for list construction
  void getNeighbor(Cell * Ci, Cell * Cj) {
    real_t dX[2];
    for (int d=0; d<2; d++) dX[d] = Ci->X[d] - Cj->X[d];        // Distance vector from source to target
    real_t R2 = norm(dX) * theta * theta;                       // Scalar distance squared
    if (R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R)) {               // If distance is far enough
    } else if (Ci->NCHILD == 0 && Cj->NCHILD == 0) {            // Else if both cells are leafs
      Ci->listP2P.push_back(Cj);                                //  Add to P2P list
    } else if (Cj->NCHILD == 0 || (Ci->R >= Cj->R && Ci->NCHILD != 0)) {// Else if Cj is leaf or Ci is larger
      for (Cell * ci=Ci->CHILD; ci!=Ci->CHILD+Ci->NCHILD; ci++) {// Loop over Ci's children
        getNeighbor(ci, Cj);                                    //   Recursive call to target child cells
      }                                                         //  End loop over Ci's children
    } else {                                                    // Else if Ci is leaf or Cj is larger
      for (Cell * cj=Cj->CHILD; cj!=Cj->CHILD+Cj->NCHILD; cj++) {//  Loop over Cj's children
        getNeighbor(Ci, cj);                                    //   Recursive call to source child cells
      }                                                         //  End loop over Cj's children
    }                                                           // End if for leafs and Ci Cj size
  }

  void addBuffer(Cells & cells, Bodies & bodies) {
    for (size_t i=0; i<cells.size(); i++) {
      if (cells[i].NCHILD == 0) {
        cells[i].NBODY = 0;
        for (size_t j=0; j<cells[i].listP2P.size(); j++) {
          Cell * Cj = cells[i].listP2P[j];
          for (int b=0; b<Cj->NBODY; b++) {
            if ((std::abs(Cj->BODY[b].X[0] - cells[i].X[0]) - cells[i].R < D)
                && (std::abs(Cj->BODY[b].X[1] - cells[i].X[1]) - cells[i].R < D)) {
              Body body;
              body.I = Cj->BODY[b].I;
              for (int d=0; d<2; d++) body.X[d] = Cj->BODY[b].X[d];
              body.q = Cj->BODY[b].q;
              body.p = Cj->BODY[b].p;
              for (int d=0; d<2; d++) body.F[d] = Cj->BODY[b].F[d];
              real_t x = fmin(cells[i].R - std::abs(body.X[0] - cells[i].X[0]), D);
              real_t y = fmin(cells[i].R - std::abs(body.X[1] - cells[i].X[1]), D);
              assert(x > -D);
              assert(y > -D);
              bodies.push_back(body);
              cells[i].NBODY++;
            }
          }
        }
        cells[i].BODY = &bodies.back() - cells[i].NBODY + 1;
      } else {
        cells[i].NBODY = 0;
      }
      cells[i].listP2P.clear();
    }
  }

  Cells buildTree(Bodies & bodies) {
    getBounds(bodies);
    Cells cells(1), jcells(1);
    cells.reserve(bodies.size());
    jcells.reserve(bodies.size());
    Bodies jbodies = bodies;
    Bodies buffer = bodies;
    buildCells(&buffer[0], &bodies[0], 0, bodies.size(), &cells[0], cells, X0, R0);
#if 0
    for (size_t i=0; i<cells.size(); i++) {
      Cell Ci = cells[i];
      if (Ci.NCHILD == 0) {
        std::cout << i << " " << Ci.X[0] << " " << Ci.X[1] << " " << Ci.R <<  std::endl;
        for (int b=0; b<Ci.NBODY; b++) {
          Body Bi = Ci.BODY[b];
          std::cout << i << " " << b << " " << Bi.X[0] << " " << Bi.X[1] << std::endl;
        }
      }
    }
#endif
    buffer = jbodies;
    buildCells(&buffer[0], &jbodies[0], 0, jbodies.size(), &jcells[0], jcells, X0, R0);
    D *= R0 / (maxlevel + 1);
    getNeighbor(&cells[0], &jcells[0]);
    bodies.clear();
    bodies.reserve(buffer.size()*27);
    addBuffer(cells, bodies);
    return cells;
  }

  void joinBuffer(Cells & cells) {
    for (size_t i=0; i<cells.size(); i++) {
      Cell * Ci = &cells[i];
      if (Ci->NCHILD == 0) {
        for (size_t j=0; j<Ci->listP2P.size(); j++) {
          Cell * Cj = Ci->listP2P[j];
          for (int bi=0; bi<Ci->NBODY; bi++) {
            Body * Bi = &Ci->BODY[bi];
            for (int bj=0; bj<Cj->NBODY; bj++) {
              Body * Bj = &Cj->BODY[bj];
              if (Bi->I == Bj->I) {
                Bi->p += Bj->p;
                if (Bi->I<10) std::cout << Bi->I << std::endl;
                for (int d=0; d<2; d++) Bi->F[d] += Bj->F[d];
              }
            }
          }
        }
      }
      cells[i].listP2P.clear();
    }
  }

  void joinBuffer(Cells & jcells, Bodies & bodies) {
    getBounds(bodies);
    Cells cells(1);
    cells.reserve(bodies.size());
    Bodies buffer = bodies;
    buildCells(&bodies[0], &buffer[0], 0, bodies.size(), &cells[0], cells, X0, R0);
    getNeighbor(&cells[0], &jcells[0]);
    joinBuffer(cells);
  }
}

#endif
