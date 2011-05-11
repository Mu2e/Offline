//
// Generate an energy value for the DIO electrons, using Wanatabe
// data, merged to Shanker's formula near the endpoint. The construction
// of the spectrum is made by ShankerWanatabeSpectrum class
//
// $Id: DIOShankerWanatabe.cc,v 1.6 2011/05/11 19:08:49 onoratog Exp $
// $Author: onoratog $
// $Date: 2011/05/11 19:08:49 $
//
// 

// C++ includes
#include <iostream>

// Framework includes
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eUtilities/inc/DIOShankerWanatabe.hh"
#include "Mu2eUtilities/inc/ShankerWanatabeSpectrum.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"

using namespace std;

namespace mu2e {

  DIOShankerWanatabe::DIOShankerWanatabe(int atomicZ, double emin, double emax, double spectRes,
                                         edm::RandomNumberGeneratorService::base_engine_t& engine):
  //atomic number of the foil material
    _Znum ( atomicZ ),
  //limits on energy generation
    _emin ( emin ),
    _emax ( emax ),
    _res ( spectRes ),
    _nBinsSpectrum ( calculateNBins() ),
    _randEnergy ( engine, &(shWaSpectrum()[0]), _nBinsSpectrum )
  {
    if (_Znum!=13) {
      throw cms::Exception("GEOM")
        << "Foil material different from Alluminum";
    }
  }

  DIOShankerWanatabe::~DIOShankerWanatabe()
  {
  }

  double DIOShankerWanatabe::fire() {

    return _emin + (_emax-_emin)*_randEnergy.fire();
 
  }

  int DIOShankerWanatabe::calculateNBins() {

    double step = _emin;
    int size = 0;
    while (step <= _emax) {
      step += _res;
      size++;
    }

    return size;

  }

  vector<double> DIOShankerWanatabe::shWaSpectrum() {

    vector<double> spectrum;
    ShankerWanatabeSpectrum WSspec(_Znum);

    double step = _emin;

    while (step <= _emax) {
      spectrum.push_back(WSspec[step]);
      step += _res;
    }

    _nBinsSpectrum = spectrum.size();

    return spectrum;
  }
  
}

