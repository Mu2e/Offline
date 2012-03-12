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
  
  // determine pValue corresponding to given mu and n (for completeness)
  double poissonCDF (double mu, int n);
  double poissonComplementaryCDF (double mu, int n); //1-CDF; more accurate

  // determine Feldman/Cousins 90% CL sensitivity given a background mean b
  double feldmanCousins90pctSensitivity ( double b );

// Helpers

  class FeldmanCousins90pctSensitivity {
    public:
      FeldmanCousins90pctSensitivity();
      double operator() (double b) const;
    private:
      static std::vector<double> sensValues_4();
      static std::vector<double> sensValues_15();
      splines::Grid<1>   grid_4;
      splines::Grid<1>   grid_15;
      splines::Spline<1> sensitivity_4;
      splines::Spline<1> sensitivity_15;
  }; // FeldmanCousins90pctSensitivity
  
} // end namespace mu2e

#endif  /* POISSON_CDF_SOLVER_H */
