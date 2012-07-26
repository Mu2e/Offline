#ifndef GeneralUtilities_safeSqrt_hh
#define GeneralUtilities_safeSqrt_hh

//
//  Take the sqrt of its argument but protect against
//  roundoff error that can take the argument negative.
//
//
// $Id: safeSqrt.hh,v 1.1 2012/07/26 18:59:07 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/26 18:59:07 $
//
// Original author Rob Kutschke
//
//  A more sophisticated version is available in
//  sqrtOrThrow.hh, which will throw if the argument
//  is too negative.
//

#include <cmath>

namespace mu2e {

  template<typename T>
  inline T safeSqrt( T x ){
    return (x>0.) ? sqrt(x) : 0.;
  }

}
#endif /* GeneralUtilities_safeSqrt_hh */
