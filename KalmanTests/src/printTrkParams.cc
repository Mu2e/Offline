//
// Free function to print the BaBar style track parameters.
//
// $Id: printTrkParams.cc,v 1.2 2011/05/18 02:27:16 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:16 $
//
// Original author Rob Kutschke
//

#include <ostream>

#include "KalmanTests/inc/printTrkParams.hh"

#include "CLHEP/Matrix/Vector.h"

using CLHEP::HepVector;

using std::ostream;
using std::endl;

namespace mu2e {

  void printTrkParams ( ostream& ost, const HepVector& par, bool doEndl ){
    ost << "( ";
    for ( size_t i=0; i<5; ++i ){
      if ( i > 0) ost << ", ";
      // HepVectors are 1...N style indexing.
      ost << par(i+1);
    }
    ost << " )";
    if ( doEndl ) ost << endl;
  }
}
