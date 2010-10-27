#ifndef EJECTEDPROTONGUN_HH
#define EJECTEDPROTONGUN_HH

//
// Simulate the protons that come from the stopping target when muons capture
// on an Al nucleus.  Use the MECO distribution for the kinetic energy of the
// protons.  
//
// $Id: EjectedProtonGun.hh,v 1.6 2010/10/27 16:42:56 onoratog Exp $
// $Author: onoratog $ 
// $Date: 2010/10/27 16:42:56 $
//
//

// Framework includes
#include "EventGenerator/inc/GeneratorBase.hh"

// Mu2e includes
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

// CLHEP includes
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"

// Forward declarations outside of namespace mu2e
class TH1D;
namespace edm {
  class Run;
}

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class EjectedProtonGun: public GeneratorBase{

  public:
    EjectedProtonGun( edm::Run& run, const SimpleConfig& config );
    virtual ~EjectedProtonGun();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    // Start: parameters that can be configured from the config file.

    double _mean;    // Mean number of protons per event
    double _elow;    // Range of proton energy.
    double _ehi;     // 
    double _czmin;   // Range of cos(polar angle)
    double _czmax;  
    double _phimin;  // Range of azimuth
    double _phimax;
    int    _nbins;   // number of bins in proton energy pdf

    double _mass; //Particle mass

    double _p; //Particle momentum

    // Limits on the generated time.
    double _tmin;
    double _tmax;

    // Histogram control.
    bool _doHistograms;

    // end: parameters that can be configured from the config file.

    //Random generators
    CLHEP::RandPoissonQ _randPoissonQ;
    RandomUnitSphere    _randomUnitSphere;
    CLHEP::RandGeneral _shape;


    TH1D* _hMultiplicity;
    TH1D* _hKE;
    TH1D* _hKEZoom;
    TH1D* _hMomentumMeV;
    TH1D* _hzPosition;
    TH1D* _hcz;
    TH1D* _hphi;
    TH1D* _htime;


    //Functions used to calculate the energy spectrum of the proton
    std::vector<double> binnedEnergySpectrum();
    double energySpectrum( double e );

  };

} // end namespace mu2e,

#endif
