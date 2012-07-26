#ifndef GeneralUtilities_sqrtOrThrow_hh
#define GeneralUtilities_sqrtOrThrow_hh
//
//  Take the sqrt of its argument but protect against
//  roundoff error that can take the argument negative.
//
//  If the argument is only a little negative, assume
//  that this is round off error and set the answer
//  to zero.  If the argument is very negative, throw.
//
//  A simpler but less detailed version is found in safeSqrt.hh
//
// $Id: sqrtOrThrow.hh,v 1.1 2012/07/26 18:59:07 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/26 18:59:07 $
//
// Original author Rob Kutschke
//

#include <cmath>

namespace mu2e {

  // Some helper functions.
  // These are found in ../src/sqrtOrThrow.cc and are used to avoid code bloat
  // in the instantiated templates.
  void throwHelperForSqrtOrThrow( double, double );
  void throwHelperForSqrtOrThrow( float,  float  );

  template<typename T>
  T sqrtOrThrow ( T x, T eps ){
    T retval(0.);
    if ( x > 0. ) {
      retval = sqrt(x);
    }else if ( x < -eps ){
      throwHelperForSqrtOrThrow(x,eps);
    }
    return retval;
  }

}
#endif /* GeneralUtilities_sqrtOrThrow_hh */
