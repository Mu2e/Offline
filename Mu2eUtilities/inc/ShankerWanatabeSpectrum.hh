#ifndef Mu2eUtilities_ShankerWanatabeSpectrum_hh
#define Mu2eUtilities_ShankerWanatabeSpectrum_hh
//
// Read Wanatabe data about DIO spectrum from a table and merge it
// with the spectrum coming from the Shanker formula

// $Id: ShankerWanatabeSpectrum.hh,v 1.5 2011/05/18 04:26:49 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/18 04:26:49 $
//
//

// C++ incldues
#include <vector>
#include <utility>

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

    int _Znum;

    double _WanaEndPoint, _WanaEndPointVal, _norm;

    std::vector<std::pair<double, double> > _wanatable;

    void ReadWanatabeTable();

    double EvaluateShanker(double E);

    double EvaluateWanatabe(double E);

    double Interpulate(double E, double e1, double p1,
                       double e2, double p2, double e3, double p3);


  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_ShankerWanatabeSpectrum_hh */

