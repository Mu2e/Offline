#ifndef Mu2eUtilities_CzarneckiSpectrum_hh
#define Mu2eUtilities_CzarneckiSpectrum_hh
//
// Read Czarnecki DIO spectrum from a table and merge it
// with the spectrum coming from the endopoint region formula

// $Id: CzarneckiSpectrum.hh,v 1.5 2013/07/12 17:17:38 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/12 17:17:38 $
//
// Original Author: Kyle Knoepfel
//

// Mu2e includes
#include "Mu2eUtilities/inc/DIOBase.hh"

namespace mu2e {

  class CzarneckiSpectrum : public DIOBase {

  public:
    
    CzarneckiSpectrum();

    double getWeight(double E) override;

  private:

    //    double FitCzarnecki(double E); maybe we'll use it later

    void readTable() override;

    double interpolate(double E, double e1, double p1,
                       double e2, double p2, double e3, double p3);  

    double interpolateE5(double E, const SpectrumValue& value);

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_CzarneckiSpectrum_hh */

