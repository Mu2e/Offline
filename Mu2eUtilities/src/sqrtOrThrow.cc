//
//  Take the sqrt of its argument but protect against
//  roundoff error that can take the argument negative.
//
//  If the argument is only a little negative, assume
//  that this is round off error and set the answer
//  to zero.  If the argument is very negative, throw.
//
// $Id: sqrtOrThrow.cc,v 1.4 2011/05/18 02:27:19 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:19 $
//
// Original author Rob Kutschke
//

#include "Mu2eUtilities/inc/sqrtOrThrow.hh"

#include <cmath>

#include <cmath>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/Exception.h"

namespace mu2e {

  double sqrtOrThrow( double x, double eps){
    if ( x > 0. ) {
      return sqrt(x);
    }else if ( x > -eps ){
      return 0.;
    }else {
      throw cet::exception("RANGE")
        << "sqrtOrThrow has an input of: "
        << x;
    }
  }


  float  sqrtOrThrow( float x,  float eps ){
    if ( x > 0. ) {
      return sqrt(x);
    }else if ( x > -eps ){
      return 0.;
    }else {
      throw cet::exception("RANGE")
        << "sqrtOrThrow has an input of: "
        << x;
    }
  }

} // end namespace mu2e
