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
// $Id: sqrtOrThrow.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//


namespace mu2e {


  double sqrtOrThrow( double x, double eps);
  float  sqrtOrThrow( float x,  float eps );

}
#endif /* Mu2eUtilities_sqrtOrThrow_hh */
