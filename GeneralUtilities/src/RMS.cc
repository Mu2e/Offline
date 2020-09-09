//
// Compute mean, RMS and some other stats about a set of numbers.
//
//
// Contact person Rob Kutschke
//

#include "GeneralUtilities/inc/RMS.hh"

#include <cmath>
#include <sstream>
#include <stdexcept>

double RMS::rms ( int npar) const {

  // Consider throwing here; if we do, update the operator<<.
  if ( n_ < npar+1 ){
    return 0.;
  }

  double tmp = sumsq_/n_ - mean()*mean();
  if ( tmp < 0. ) return 0;

  double fac= double(n_)/double(n_-npar);

  return sqrt(fac*tmp);
}

double RMS::rms0 ( int npar) const {

  // Consider throwing here; if we do, update the operator<<.
  if ( n_ <= npar ){
    return 0.;
  }

  double tmp = sumsq_/n_;
  if ( tmp < 0. ) return 0;

  double fac= double(n_)/double(n_-npar);

  return sqrt(fac*tmp);
}

double RMS::errorMean( int npar ) const{

  // Consider throwing here.
  if ( n_ == 0 ) {
    return 0;
  }

  return rms(npar)/sqrt(n_);
}

double RMS::errorRMS( int npar ) const{

  // Consider throwing here.
  if ( n_ == 0 ) {
    return 0;
  }

  return rms(npar)/sqrt(2.*n_);
}
