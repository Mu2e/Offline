#ifndef _matrixUtils_h_
#define _matrixUtils_h_
//
// Joe Boudreau September 1997
//
// This routine is meant to remedy a defect in the CLHEP Matrix class,
// namely, the symmetric matrix diagonalization routine fails to give
// a sensible orthogonal transformation when diagonalizing even pretty
// simple matrices.
//
// This is the jacobi method that should always work with real symmetric
// matrices.
class HepMatrix;
class HepSymMatrix;
inline HepMatrix jacobi(HepSymMatrix * A);
#include "CDFFitting/inc/matrixUtils.icc"
#endif
