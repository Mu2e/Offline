#ifndef EventGenerator_EjectedNeutronGun_hh
#define EventGenerator_EjectedNeutronGun_hh

//
// Simulate the neutrons that come from the stopping target when muons capture
// on an Al nucleus.  Use the MARS distribution for the kinetic energy of the
// neutron.
//
// $Id: EjectedNeutronGun.hh,v 1.9 2011/08/12 19:33:57 vbiliyar Exp $
// $Author: vbiliyar $
// $Date: 2011/08/12 19:33:57 $
//
//

// C++ includes
#include <memory>

// Mu2e includes
#include "EventGenerator/inc/FoilParticleGenerator.hh"
#include "EventGenerator/inc/GeneratorBase.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

// CLHEP includes
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandPoissonQ.h"

// Forward declarations outside of namespace mu2e
class TH1D;
namespace art {
  class Run;
}

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class EjectedNeutronGun: public GeneratorBase{

  public:
    EjectedNeutronGun( art::Run& run, const SimpleConfig& config );
    virtual ~EjectedNeutronGun();

    virtual void generate( GenParticleCollection&  );

  private:

    // Start: parameters that can be configured from the config file.

    enum SpectrumType { mars , docdb1619 , last_enum };
    
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
    
    // Function to check that a valid neutron spectrum model is being used
    SpectrumType checkSpectrumModel(const int i);
    SpectrumType _spectrumModel;

    // end: parameters that can be configured from the config file.

    //Random generators
    CLHEP::RandPoissonQ _randPoissonQ;
    RandomUnitSphere    _randomUnitSphere;
    CLHEP::RandGeneral _shape;

    int _nToSkip;

    TH1D* _hMultiplicity;
    TH1D* _hKE;
    TH1D* _hKEZoom;
    TH1D* _hMomentumMeV;
    TH1D* _hzPosition;
    TH1D* _hcz;
    TH1D* _hphi;
    TH1D* _htime;
    TH1D* _hmudelay;
    TH1D* _hpulsedelay;
    
    
    
    // Continuous spectrum function
    double energySpectrum(const double e);

    //Functions used to calculate the energy spectrum of the neutron
    std::vector<double> binnedEnergySpectrum();

  };

} // end namespace mu2e,

#endif /* EventGenerator_EjectedNeutronGun_hh */
