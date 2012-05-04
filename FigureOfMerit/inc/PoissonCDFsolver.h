#ifndef POISSON_CDF_SOLVER_H
#define POISSON_CDF_SOLVER_H

// ----------------------------------------------------------------------
//
// PoissonCDFsolver.h
//
// Functions to determine values of mu and n leading to given p-values
// based on the Poisson distribution 
//
// ----------------------------------------------------------------------

#include "splines/Spline.h"

#include <vector>

namespace mu2e {

  // determine, for a given mu, the first n such that CDF > a given p-value
  int poissonInverseCDF (double mu, double pValue);
  
  // determine, for a given p-value and n, the mu such that CDF at that n
  // is equal to that p-value
  double poissonMeanSolver (int n, double pValue);

  // determine, for a given p-value and (n1,n2), the mu such that the prob
  // slice from n1 to n2 is is equal to that p-value.  This takes the hint that
  // the solution is in the interval [mu_1,mu_2]
  double poissonMeanSliceSolver (int n1, int n2, double pValue, 
                                 double mu_1, double mu_2);

  // determine pValue corresponding to given mu and n (for completeness)
  double poissonCDF (double mu, int n);
  double poissonComplementaryCDF (double mu, int n); //1-CDF; more accurate
  double poissonProbabilitySlice (double mu, int n1, int n2);

  
} // end namespace mu2e

#endif  /* POISSON_CDF_SOLVER_H */
