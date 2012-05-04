//
// Generate an energy value for the DIO electrons, using a generic 
// custom distribution as momentum spectrum
// The construction
// of the spectrum is made by specialized classes
//
// $Id: ReadDIOSpectrum.cc,v 1.4 2012/05/04 20:12:16 onoratog Exp $
// $Author: onoratog $
// $Date: 2012/05/04 20:12:16 $
//
// Original author: Gianni Onorato
//

// C++ includes
#include <iostream>

// Framework includes
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eUtilities/inc/ReadDIOSpectrum.hh"
#include "Mu2eUtilities/inc/CzarneckiSpectrum.hh"
#include "Mu2eUtilities/inc/ShankerWanatabeSpectrum.hh"
#include "MCDataProducts/inc/PDGCode.hh"

using namespace std;

namespace mu2e {

  ReadDIOSpectrum::ReadDIOSpectrum(int atomicZ, double mumass, double emass, double emin, double emax, 
                                   double spectRes, string spectrum,
                                   art::RandomNumberGenerator::base_engine_t& engine):
  //atomic number of the foil material
    _znum ( atomicZ ),
    _mumass ( mumass ),
    _emass ( emass ),
  //limits on energy generation
    _emin ( emin ),
    _emax ( emax ),
    _res ( spectRes ),
    _nBinsSpectrum ( calculateNBins() ),
    _spectrum ( spectrum ),
    _randEnergy ( engine, &(readSpectrum()[0]), _nBinsSpectrum )
  {
    if (_znum!=13) {
      throw cet::exception("GEOM")
        << "Foil material different from Alluminum";
    }
  }

  ReadDIOSpectrum::~ReadDIOSpectrum()
  {
  }




  double ReadDIOSpectrum::fire() {

    double e = _randEnergy.fire();
    
    return _emin + (_emax-_emin)*e;

  }

  int ReadDIOSpectrum::calculateNBins() {

    double step = _emin;
    int size = 0;
    while (step <= _emax) {
      step += _res;
      size++;
    }

    return size;

  }

  vector<double> ReadDIOSpectrum::readSpectrum() {

    vector<double> spectrum;

    if (_spectrum == "ShankerWanatabe") {
      ShankerWanatabeSpectrum WSspec(_znum, _mumass, _emass);
      
      double step = _emin;
      
      while (step <= _emax) {
        spectrum.push_back(WSspec(step));
        step += _res;
      }
    } else if (_spectrum == "Czarnecki") {

      CzarneckiSpectrum CZspec(_znum);
      
      double step = _emin;
      
      while (step <= _emax) {
        spectrum.push_back(CZspec(step));
        step += _res;

      }
    } else {
      throw cet::exception("MODEL")
        << "Wrong or not allowed DIO energy spectrum";
    }

    //    vector<double>::iterator it = spectrum.begin();
    //while (it!=spectrum.end()) {
    //  cout << *it << endl;
    //  it++;
    //}


    _nBinsSpectrum = spectrum.size();

    return spectrum;
  }
  
}

