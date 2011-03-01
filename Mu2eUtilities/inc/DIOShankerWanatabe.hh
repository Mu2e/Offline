#ifndef DIOSHANKERWANATABE_HH
#define DIOSHANKERWANATABE_HH
//
// Generate a momentum for the DIO electrons, using Wanatabe
// data, merged to Shanker's formula near the endpoint
//
// $Id: DIOShankerWanatabe.hh,v 1.2 2011/03/01 04:38:33 onoratog Exp $
// $Author: onoratog $
// $Date: 2011/03/01 04:38:33 $
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

    CLHEP::RandGeneral _randEnergy;

    std::vector<double> ShWaSpectrum();

    int _nBinsSpectrum;
  };

} // end of namespace mu2e

#endif

