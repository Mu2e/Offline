// Mu2e includes
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "Mu2eUtilities/inc/MueXSpectrum.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

// Framework includes
#include "cetlib/pow.h"

// CLHEP includes
#include "CLHEP/Random/RandFlat.h"

// C++ includes
#include <iostream>
#include <cmath>

//GSL includes
#include "stdio.h"
#include "gsl/gsl_integration.h"

using namespace std;
  
namespace mu2e {
 
  MueXSpectrum::MueXSpectrum(double maxEnergy, double bin, int RadCorrected) :
    _bin          (bin         ),
    _spectrumType (RadCorrected) 
  {
    GlobalConstantsHandle<ParticleDataTable> pdt;

    _par.me    = pdt->particle(PDGCode::e_minus ).ref().mass().value();
    _par.eMax  = maxEnergy;
    _nbins     = maxEnergy/_bin;
    double de = _bin;
    if (_nbins*_bin < maxEnergy){
      _nbins += 1;
      de = _par.eMax-_bin*(_nbins-1);
    }
     _integral = evalIntegral(de); 

  }
    
  double MueXSpectrum::f(double E, void *p) { 
    double eMax  = ((MueXSpectrum::Params_t*) p)->eMax;
    double me    = ((MueXSpectrum::Params_t*) p)->me;
   
    double f     = 1 //TODO
    if (f < 0) f = 0;
    return f;
  }

  double MueXSpectrum::getCorrectedMueXSpectrum(double e) const {
    return MueXSpectrum::_f(e,(void*) &_par);
  }

  double MueXSpectrum::getWeight(double E) const {
    
    double weight(0.);
  
    int    bin = E/_bin ; 

    if (bin < _nbins-1) {
      weight = _bin* getCorrectedMueXSpectrum(E);
    }
    else {
      weight =( 1.-_integral);
    }
    
    return weight;
  }

//TODO  
  double MueXSpectrum::evalIntegral(double de){
   return 1.
  }
}
