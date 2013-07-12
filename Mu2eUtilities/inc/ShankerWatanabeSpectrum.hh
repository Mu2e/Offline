#ifndef Mu2eUtilities_ShankerWatanabeSpectrum_hh
#define Mu2eUtilities_ShankerWatanabeSpectrum_hh
//
// Read Watanabe data about DIO spectrum from a table and merge it
// with the spectrum coming from the Shanker formula

// $Id: ShankerWatanabeSpectrum.hh,v 1.1 2013/07/12 17:17:38 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/12 17:17:38 $
//
//

// Mu2e includes
#include "Mu2eUtilities/inc/DIOBase.hh"

namespace mu2e {

  class ShankerWatanabeSpectrum : public DIOBase {

  public:

    ShankerWatanabeSpectrum();
    ~ShankerWatanabeSpectrum(){}

    double getWeight(double E) override;

  private:

    double _wanaEndPoint, _wanaEndPointVal, _norm;

    // built-in energy tolerance (in MeV)
    const double _tolerance = 0.0049;

    void readTable() override;

    double evaluateShanker (double E);
    double evaluateWatanabe(double E);

    double interpolate(double E, double e1, double p1,
                       double e2, double p2, double e3, double p3);


  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_ShankerWatanabeSpectrum_hh */

