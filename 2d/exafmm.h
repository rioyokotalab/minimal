#ifndef exafmm_h
#define exafmm_h
#include <cassert>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace exafmm {
  //! Basic type definitions
  typedef double real_t;                                        //!< Floating point type is single precision
  typedef std::complex<real_t> complex_t;                       //!< Complex type

  //! Structure of bodies
  struct Body {
    int I;                                                      //!< Index
    real_t X[2];                                                //!< Position
    real_t q;                                                   //!< Charge
    real_t p;                                                   //!< Potential
    real_t F[2];                                                //!< Force
  };
  typedef std::vector<Body> Bodies;                             //!< Vector of bodies

  //! Structure of cells
  struct Cell {
    int NCHILD;                                                 //!< Number of child cells
    int NBODY;                                                  //!< Number of descendant bodies
    Cell * CHILD;                                               //!< Pointer of first child cell
    Body * BODY;                                                //!< Pointer of first body
    real_t X[2];                                                //!< Cell center
    real_t R;                                                   //!< Cell radius
    std::vector<Cell*> listM2L;                                 //!< M2L interaction list
    std::vector<Cell*> listP2P;                                 //!< P2P interaction list
    std::vector<complex_t> M;                                   //!< Multipole expansion coefficients
    std::vector<complex_t> L;                                   //!< Local expansion coefficients
  };
  typedef std::vector<Cell> Cells;                              //!< Vector of cells

  //! Global variables
  int P;                                                        //!< Order of expansions
  int ncrit;                                                    //!< Number of bodies per leaf cell
  int maxlevel;                                                 //!< Maximum number of levels in tree
  real_t D;                                                     //!< Buffer size
  real_t theta;                                                 //!< Multipole acceptance criterion
  real_t R0;                                                    //!< Radius of root cell
  real_t X0[2];                                                 //!< Center of root cell
  real_t dX[2];                                                 //!< Distance vector
#pragma omp threadprivate(dX)                                   //!< Make global variables private

  //!< L2 norm of vector X
  inline real_t norm(real_t * X) {
    return X[0] * X[0] + X[1] * X[1];                           // L2 norm
  }
}

#endif
