#ifndef EventGenerator_EjectedPhotonGun_hh
#define EventGenerator_EjectedPhotonGun_hh

//
// Simulate the photons coming from the stopping target when muons are captured
// by an Al nucleus.
//
// $Id: EjectedPhotonGun.hh,v 1.7 2013/05/31 18:07:29 gandr Exp $
// $Author: &
// $Date: &
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
class TH2D;
namespace art {
  class Run;
}

namespace mu2e {

  // Forward reference.
  class SimpleConfig;
  class DetectorSystem;

  class EjectedPhotonGun: public GeneratorBase{

  public:
    EjectedPhotonGun( art::Run& run, const SimpleConfig& config );
    virtual ~EjectedPhotonGun();

    virtual void generate( GenParticleCollection&  );

  private:

    // Start: parameters that can be configured from the config file.

    double _mean;    // Mean number of protons per event
    double _elow;    // Range of proton energy.
    double _ehi;     //
    double _czmin;   // Range of cos(polar angle)
    double _czmax;
    double _phimin;  // Range of azimuth
    double _phimax;

    // Class object to generate position and time of the particle
    std::unique_ptr<FoilParticleGenerator> _fGenerator;

    double _p; //Particle momentum

    // Limits on the generated time.
    double _tmin;
    double _tmax;

    // Histogram control.
    bool _doHistograms;

    bool _PStoDSDelay;
    bool _pPulseDelay;
    double _pPulseShift;

    // Activate the folding procedure on generation time. Default is on
    bool _timeFolding;

    // Select the position, type and time type for the generation
    std::string _foilGen;
    std::string _posGen;
    std::string _timeGen;

    // end: parameters that can be configured from the config file.

    //Random generators
    CLHEP::RandPoissonQ _randPoissonQ;
    RandomUnitSphere    _randomUnitSphere;
    CLHEP::RandFlat _flatmomentum;
    std::string _STfname;
    int _nToSkip;

    const DetectorSystem *_detSys;

    TH1D* _hMultiplicity;
    TH1D* _hKE;
    TH1D* _hKEZoom;
    TH1D* _hMomentumMeV;
    TH1D* _hradius;
    TH1D* _hzPosition;
    TH1D* _hcz;
    TH1D* _hphi;
    TH1D* _htime;
    TH1D* _hmudelay;
    TH1D* _hpulsedelay;
    TH2D* _hyx;
    TH2D* _hrz;

  };

} // end namespace mu2e,

#endif /* EventGenerator_EjectedPhotonGun_hh */
