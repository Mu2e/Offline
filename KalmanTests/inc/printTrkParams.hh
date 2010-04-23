#ifndef KalmanTests_printTrkParams_hh
#define KalmanTests_printTrkParams_hh
//
// Free function to print the BaBar style track parameters.
//
// $Id: printTrkParams.hh,v 1.1 2010/04/23 20:18:57 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/04/23 20:18:57 $
//
// Original author Rob Kutschke
//
#include <iostream>

namespace CLHEP{
  class HepVector;
}

namespace mu2e {
  void printTrkParams ( std::ostream& ost, const CLHEP::HepVector& par, bool doEndl=true );

  inline void printTrkParams ( const CLHEP::HepVector& par ){
    printTrkParams ( std::cout, par );
  }

}
#endif
