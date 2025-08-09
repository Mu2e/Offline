// Spectrum based on arXiv: 1110.2874
// Sophie Middleton, 2021
// Mu2e includes
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/Mu2eUtilities/inc/Mu2eXSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"

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

  Mu2eXSpectrum::Mu2eXSpectrum(double maxEnergy, double bin, int RadCorrected):
    _bin (bin),
    _spectrumType (RadCorrected)
  {
    GlobalConstantsHandle<ParticleDataList> pdt;
    _par.me    = pdt->particle(PDGCode::e_minus ).mass();
    _par.eMax  = maxEnergy;
    _par.mmu = 105.6584; //mu mass MeV
    _par.Emu = 105.194; //MeV
    _par.BR = 5e-5; //Branching ratio
    //_par.Gamma = 2.99561e-16;
    _par.mN = 25133; //Mass of Al in MeV
    // Following from arXiv:1110.2874:
    _par.a0 = 3.289e-10;
    _par.a1 = 3.137e-7;
    _par.a2 = 1.027e-4;
    _par.a3 = 1.438e-3;
    _par.a4 = 2.419e-3;
    _par.a5 = 1.215e-1;
    _nbins     = maxEnergy/_bin;
    double de = _bin;
    if (_nbins*_bin < maxEnergy){
      _nbins += 1;
      de = _par.eMax-_bin*(_nbins-1);
    }
    _integral = evalIntegral(de);
  }

  double Mu2eXSpectrum::f_mu2eX_gen(double E, void *p) { //For E>100MeV Only
    double mmu   = ((Mu2eXSpectrum::Params_t*) p)->mmu;
    double Emu   = ((Mu2eXSpectrum::Params_t*) p)->Emu;
    double BR    = ((Mu2eXSpectrum::Params_t*) p)->BR;
    double mN    = ((Mu2eXSpectrum::Params_t*) p)->mN;
    double a0    = ((Mu2eXSpectrum::Params_t*) p)->a0;
    double a1    = ((Mu2eXSpectrum::Params_t*) p)->a1;
    double a2    = ((Mu2eXSpectrum::Params_t*) p)->a2;
    double a3    = ((Mu2eXSpectrum::Params_t*) p)->a3;
    double a4    = ((Mu2eXSpectrum::Params_t*) p)->a4;
    double a5    = ((Mu2eXSpectrum::Params_t*) p)->a5;
    double delta = ((Emu - E - pow(E,2)/(2*mN))/mmu);
    double f     = BR*(1/mmu)*(a0*pow(delta,1) + a1*pow(delta,2) + a2*pow(delta,3) + a3*pow(delta,4) + a4*pow(delta,5) + a5*pow(delta,6));
    if (f < 0) f = 0;
    return f; //F = 1/Gamma * d(Gamma)/dEe
  }

  double Mu2eXSpectrum::getCorrectedMu2eXSpectrum(double e) const {
    return Mu2eXSpectrum::f_mu2eX_gen(e,(void*) &_par);
  }

  /* For compatibility - perhaps we dont need this? */
  double Mu2eXSpectrum::getWeight(double E) const {

    double weight(0.);
    int    bin = E/_bin ;
    if (bin < _nbins-1) {
      weight = _bin* getCorrectedMu2eXSpectrum(E);
    }
    else {
      weight =( 1.-_integral);
    }

    return weight;
  }


  double Mu2eXSpectrum::evalIntegral(double de){
    gsl_function F;
    F.function = &f_mu2eX_gen;
    F.params   = &_par;

    size_t limit  = 1000;
    double epsabs = 0.001;
    double epsrel = 0.001;

    gsl_integration_workspace * ws = gsl_integration_workspace_alloc(10000);

    double result, abserr;

    double emin = 100; //eqn only valid for > 100 MeV
    double emax = _par.eMax-de;

    gsl_integration_qags(&F,
      emin, emax,
      epsabs,
      epsrel,
      limit,
      ws,
      &result,
      &abserr);

    gsl_integration_workspace_free(ws);
    return result;
  }
}
