#ifndef DIOSHANKERWANATABE_HH
#define DIOSHANKERWANATABE_HH
//
// Generate a momentum for the DIO electrons, using Wanatabe
// data, merged to Shanker's formula near the endpoint
//
// $Id: DIOShankerWanatabe.hh,v 1.4 2011/05/11 19:08:49 onoratog Exp $
// $Author: onoratog $
// $Date: 2011/05/11 19:08:49 $
//
// 

// C++ incldues
#include <vector>

// Framework includes
#include "FWCore/Services/interface/RandomNumberGeneratorService.h"

// Mu2e includes
#include "Mu2eUtilities/inc/DIOBase.hh"

//CLHEP includes
#include "CLHEP/Random/RandGeneral.h"

namespace mu2e {

  class DIOShankerWanatabe: public DIOBase {
    
  public:

    DIOShankerWanatabe(int atomicZ, double emin, double emax, double spectRes,
                       edm::RandomNumberGeneratorService::base_engine_t & engine);
    
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

#endif

