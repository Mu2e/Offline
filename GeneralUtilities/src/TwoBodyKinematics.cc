//
// Kinematics of 2 body decay.
//
// $Id: TwoBodyKinematics.cc,v 1.5 2011/05/25 18:31:42 wb Exp $
// $Author: wb $
// $Date: 2011/05/25 18:31:42 $
//
// Original author Rob Kutschke

#include "GeneralUtilities/inc/TwoBodyKinematics.hh"

#include "cetlib/pow.h"
#include <cmath>
#include <sstream>
#include <stdexcept>

using cet::diff_of_squares;
using cet::sum_of_squares;
using std::sqrt;

TwoBodyKinematics::TwoBodyKinematics( double m0, double m1, double m2):
  m0_(m0),
  m1_(m1),
  m2_(m2){

  // Numerator of the expression for the momentum.
  double num = diff_of_squares(m0_, m1_+m2 ) * diff_of_squares(m0_, m1_-m2_);
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
  e1_ = sqrt( sum_of_squares(m1_, p_) );
  e2_ = sqrt( sum_of_squares(m2_, p_) );
}
