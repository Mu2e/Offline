#ifndef EventGenerator_EjectedNeutronGun_hh
#define EventGenerator_EjectedNeutronGun_hh

//
// Simulate the neutrons that come from the stopping target when muons capture
// on an Al nucleus.  Use the MARS distribution for the kinetic energy of the
// neutron.
//
// $Id: EjectedNeutronGun.hh,v 1.17 2013/05/31 18:07:29 gandr Exp $
// $Author: gandr $
// $Date: 2013/05/31 18:07:29 $
//
//

// C++ includes
#include <memory>
#include <string>

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
  class DetectorSystem;

  class EjectedNeutronGun: public GeneratorBase{

  public:
    EjectedNeutronGun( art::Run& run, const SimpleConfig& config );
    virtual ~EjectedNeutronGun();

    virtual void generate( GenParticleCollection&  );

  private:

    // Start: parameters that can be configured from the config file.

    // Which model of the energy spectrum will be used?
    enum SpectrumType { mars , docdb1619 , last_enum };

    double _mean;    // Mean number of neutrons per event
    double _kelow;   // Range of kinetic energy for the neutrons.
    double _kehi;    //
    double _czmin;   // Range of cos(polar angle)
    double _czmax;
    double _phimin;  // Range of azimuth
    double _phimax;
    bool   _PStoDSDelay;
    bool   _pPulseDelay;
    double _pPulseShift;

    // Activate the folding procedure on generation time. Default is on
    bool _timeFolding;

    // Select the position, type and time type for the generation
    std::string _foilGen;
    std::string _posGen;
    std::string _timeGen;

    // Number of bins in neutron energy pdf
    int  _nbins;

    // Limits on the generated time.
    double _tmin;
    double _tmax;

    // Histogram control.
    bool _doHistograms;

    // Function to check that a valid neutron spectrum model is being used
    SpectrumType checkSpectrumModel(const int i);
    SpectrumType _spectrumModel;

    // end: parameters that can be configured from the config file.

    // Class object to generate position and time of the particle
    std::unique_ptr<FoilParticleGenerator> _fGenerator;

    //Particle mass
    double _mass;

    std::string _filetoread;

    //Random generators
    CLHEP::RandPoissonQ _randPoissonQ;
    RandomUnitSphere    _randomUnitSphere;
    CLHEP::RandGeneral _shape;
    std::string _STfname;
    int _nToSkip;

    const DetectorSystem *_detSys;

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
