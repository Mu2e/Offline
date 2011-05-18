#ifndef Mu2eUtilities_DIOShankerWanatabe_hh
#define Mu2eUtilities_DIOShankerWanatabe_hh
//
// Generate a momentum for the DIO electrons, using Wanatabe
// data, merged to Shanker's formula near the endpoint
//
// $Id: DIOShankerWanatabe.hh,v 1.9 2011/05/18 16:21:55 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/18 16:21:55 $
//
//

// C++ includes
#include <vector>

// Framework includes
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// Mu2e includes
#include "Mu2eUtilities/inc/DIOBase.hh"

//CLHEP includes
#include "CLHEP/Random/RandGeneral.h"

namespace mu2e {

  class DIOShankerWanatabe: public DIOBase {

  public:

    DIOShankerWanatabe(int atomicZ, double emin, double emax, double spectRes,
                       art::RandomNumberGenerator::base_engine_t & engine);

    ~DIOShankerWanatabe();

    double fire();

  private:

    int _Znum;

    double _emin, _emax, _res;

    int calculateNBins();

    int _nBinsSpectrum;

    CLHEP::RandGeneral _randEnergy;

    std::vector<double> shWaSpectrum();

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_DIOShankerWanatabe_hh */

