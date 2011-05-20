//
// Kinematics of 2 body decay.
//
// $Id: TwoBodyKinematics.cc,v 1.4 2011/05/20 22:39:28 wb Exp $
// $Author: wb $
// $Date: 2011/05/20 22:39:28 $
//
// Original author Rob Kutschke

#include "GeneralUtilities/inc/TwoBodyKinematics.hh"

#include "cetlib/pow.h"
#include <cmath>
#include <sstream>
#include <stdexcept>

using cet::square;
using std::sqrt;

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
