#ifndef EJECTEDNEUTRONGUN_HH
#define EJECTEDNEUTRONGUN_HH

//
// Simulate the neutrons that come from the stopping target when muons capture
// on an Al nucleus.  Use the MARS distribution for the kinetic energy of the
// neutron.  
//
// $Id: EjectedNeutronGun.hh,v 1.2 2011/05/17 06:00:47 onoratog Exp $
// $Author: onoratog $
// $Date: 2011/05/17 06:00:47 $
//
//

// C++ includes
#include <memory>

// Mu2e includes
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "EventGenerator/inc/GeneratorBase.hh"
#include "EventGenerator/inc/FoilParticleGenerator.hh"

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

  class EjectedNeutronGun: public GeneratorBase{

  public:
    EjectedNeutronGun( edm::Run& run, const SimpleConfig& config );
    virtual ~EjectedNeutronGun();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    // Start: parameters that can be configured from the config file.

    double _mean;    // Mean number of neutrons per event
    double _elow;    // Range of neutrons energy.
    double _ehi;     // 
    double _czmin;   // Range of cos(polar angle)
    double _czmax;  
    double _phimin;  // Range of azimuth
    double _phimax;
    bool _PStoDSDelay;
    bool _pPulseDelay;
    int    _nbins;   // number of bins in neutron energy pdf

    // Class object to generate position and time of the particle
    std::auto_ptr<FoilParticleGenerator> _fGenerator;

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


    //Functions used to calculate the energy spectrum of the neutron
    std::vector<double> binnedEnergySpectrum();

  };

} // end namespace mu2e,

#endif
