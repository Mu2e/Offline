//
// Kinematics of 2 body decay.
//
// $Id: TwoBodyKinematics.cc,v 1.3 2011/05/18 20:09:10 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 20:09:10 $
//
// Original author Rob Kutschke

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "GeneralUtilities/inc/TwoBodyKinematics.hh"
#include "GeneralUtilities/inc/pow.hh"

TwoBodyKinematics::TwoBodyKinematics( double m0, double m1, double m2):
  m0_(m0),
  m1_(m1),
  m2_(m2){

  // Numerator of the expression for the momentum.
  double num = ( square(m0_) - square(m1_+m2) )*( square(m0_) - square(m1_-m2_) );
  if ( num < 0. ){
    std::ostringstream out;
    out << "TwoBodyKinematics: masses of daughters exceed that of mother: \n"
        << "   Mother:     " << m0_ << " "
        << "   Daughter 1: " << m1_ << " "
        << "   Daughter 2: " << m2_ << " "
        << "   Numerator:  " << num;
    throw std::invalid_argument( out.str() );
  }

  p_ = sqrt(num)/2./m0_;
  e1_ = sqrt( square(m1_) + square(p_) );
  e2_ = sqrt( square(m2_) + square(p_) );
}
