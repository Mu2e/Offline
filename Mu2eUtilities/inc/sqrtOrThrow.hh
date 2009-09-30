#ifndef Mu2eUtilities_sqrtOrThrow_HH
#define Mu2eUtilities_sqrtOrThrow_HH

//
//  Take the sqrt of its argument but protect against
//  roundoff error that can take the argument negative.
//
//  If the argument is only a little negative, assume
//  that this is round off error and set the answer
//  to zero.  If the argument is very negative, throw.
//
// $Id: sqrtOrThrow.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//


namespace mu2e {


  double sqrtOrThrow( double x, double eps);
  float  sqrtOrThrow( float x,  float eps );

}
#endif
