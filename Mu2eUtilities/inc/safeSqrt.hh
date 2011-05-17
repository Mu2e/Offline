#ifndef Mu2eUtilities_safeSqrt_hh
#define Mu2eUtilities_safeSqrt_hh

//
//  Take the sqrt of its argument but protect against
//  roundoff error that can take the argument negative.
//
//
// $Id: safeSqrt.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:36 $
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
#endif /* Mu2eUtilities_safeSqrt_hh */
