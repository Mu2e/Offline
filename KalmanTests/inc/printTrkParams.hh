#ifndef KalmanTests_printTrkParams_hh
#define KalmanTests_printTrkParams_hh
//
// Free function to print the BaBar style track parameters.
//
// $Id: printTrkParams.hh,v 1.2 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:35 $
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
