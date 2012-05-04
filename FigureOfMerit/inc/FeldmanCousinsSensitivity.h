#ifndef FELDMAN_COUSINS_SENSITIVITY_H
#define FELDMAN_COUSINS_SENSITIVITY_H

// ----------------------------------------------------------------------
//
// FeldmanCousinsSensitivity.h
//
// Functions to determine values of sensitivity for a given background b 
//
// ----------------------------------------------------------------------

#include "splines/Spline.h"

#include <vector>

namespace mu2e {


  // determine Feldman/Cousins 90% CL sensitivity given a background mean b
  double feldmanCousins90pctSensitivity ( double b );

// Helpers

  class FeldmanCousins90pctSensitivity {
    public:
      FeldmanCousins90pctSensitivity();
      double operator() (double b) const;
    private:
      static std::vector<double> sensValues_10();
      splines::Grid<1>   grid_10;
      splines::Spline<1> sensitivity_10;
  }; // FeldmanCousins90pctSensitivity
  
} // end namespace mu2e

#endif  /* POISSON_CDF_SOLVER_H */
