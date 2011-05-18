#ifndef KalmanTests_printTrkParams_hh
#define KalmanTests_printTrkParams_hh
//
// Free function to print the BaBar style track parameters.
//
// $Id: printTrkParams.hh,v 1.3 2011/05/18 02:27:16 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:16 $
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
#endif /* KalmanTests_printTrkParams_hh */
