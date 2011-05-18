#ifndef Mu2eUtilities_sqrtOrThrow_hh
#define Mu2eUtilities_sqrtOrThrow_hh

//
//  Take the sqrt of its argument but protect against
//  roundoff error that can take the argument negative.
//
//  If the argument is only a little negative, assume
//  that this is round off error and set the answer
//  to zero.  If the argument is very negative, throw.
//
// $Id: sqrtOrThrow.hh,v 1.3 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:18 $
//
// Original author Rob Kutschke
//


namespace mu2e {


  double sqrtOrThrow( double x, double eps);
  float  sqrtOrThrow( float x,  float eps );

}
#endif /* Mu2eUtilities_sqrtOrThrow_hh */
