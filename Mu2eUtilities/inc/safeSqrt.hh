#ifndef Mu2eUtilities_safesqrt_HH
#define Mu2eUtilities_safesqrt_HH

//
//  Take the sqrt of its argument but protect against
//  roundoff error that can take the argument negative.
//
//
// $Id: safeSqrt.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
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
#endif
