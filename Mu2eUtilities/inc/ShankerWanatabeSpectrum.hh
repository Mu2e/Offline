#ifndef Mu2eUtilities_ShankerWanatabeSpectrum_hh
#define Mu2eUtilities_ShankerWanatabeSpectrum_hh
//
// Read Wanatabe data about DIO spectrum from a table and merge it
// with the spectrum coming from the Shanker formula

// $Id: ShankerWanatabeSpectrum.hh,v 1.8 2012/02/07 07:17:09 onoratog Exp $
// $Author: onoratog $
// $Date: 2012/02/07 07:17:09 $
//
//

// C++ includes
#include <utility>
#include <vector>

// Framework includes
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// Mu2e includes
#include "Mu2eUtilities/inc/DIOBase.hh"

//CLHEP includes
#include "CLHEP/Random/RandGeneral.h"

namespace mu2e {

  class ShankerWanatabeSpectrum {

  public:

    ShankerWanatabeSpectrum(int atomicZ);

    ~ShankerWanatabeSpectrum();

    double operator[](double E);

  private:

    int _znum;

    double _wanaEndPoint, _wanaEndPointVal, _norm;

    std::vector<std::pair<double, double> > _wanatable;

    void readWanatabeTable();

    double evaluateShanker(double E);

    double evaluateWanatabe(double E);

    double interpulate(double E, double e1, double p1,
                       double e2, double p2, double e3, double p3);


  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_ShankerWanatabeSpectrum_hh */

