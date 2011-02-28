#ifndef DIOSHANKERWANATABE_HH
#define DIOSHANKERWANATABE_HH
//
// Generate a momentum for the DIO electrons, using Wanatabe
// data, merged to Shanker's formula near the endpoint
//
// $Id: DIOShankerWanatabe.hh,v 1.1 2011/02/28 16:18:50 onoratog Exp $
// $Author: onoratog $
// $Date: 2011/02/28 16:18:50 $
//
// 

// C++ incldues
#include <vector>
#include <list>
#include <utility>

// Framework includes
#include "FWCore/Services/interface/RandomNumberGeneratorService.h"

// Mu2e includes
#include "Mu2eUtilities/inc/DIOBase.hh"

//CLHEP includes
#include "CLHEP/Random/RandGeneral.h"

namespace mu2e {

  class DIOShankerWanatabe: public DIOBase {
    
  public:

    DIOShankerWanatabe(int atomicZ, double emin, double emax,
                       edm::RandomNumberGeneratorService::base_engine_t & engine);

    ~DIOShankerWanatabe();

    double fire();

  private:

    int _Znum;

    double _emin, _emax;

    CLHEP::RandGeneral _randEnergy;

    std::vector<double> ShWaSpectrum();

    void AddShanker(std::pair<double, double> wEnd, std::list<std::pair<double, double> > & shankList);

    double EvalShanker(double E);

    double _startpoint;
    double _endpoint;
    int _nBinsSpectrum;
  };

} // end of namespace mu2e

#endif

