//
// Free function to print the BaBar style track parameters.
//
// $Id: printTrkParams.cc,v 1.1 2010/04/23 20:18:57 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/04/23 20:18:57 $
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
