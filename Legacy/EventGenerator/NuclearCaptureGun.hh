#ifndef EventGenerator_NuclearCaptureGun_hh
#define EventGenerator_NuclearCaptureGun_hh

//
// Simulate the complete process of the nuclear capture of muons by aluminum atoms
// which results in protons, neutrons and photons
//
//
// $Id: NuclearCaptureGun.hh,v 1.12 2013/05/31 18:07:29 gandr Exp $
// $Author: gandr $
// $Date: 2013/05/31 18:07:29 $
//
// Original author Gianni Onorato

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
  class DetectorSystem;

  class NuclearCaptureGun: public GeneratorBase{

  public:
    NuclearCaptureGun( art::Run& run, const SimpleConfig& config );
    virtual ~NuclearCaptureGun();

    virtual void generate( GenParticleCollection&  );

  private:

    // Start: parameters that can be configured from the config file.

    double _mean; // mean number of nuclear captures per event
    double _protonMean; // mean number of ejected protons per nuclear capture
    double _neutronMean; // mean number of ejected neutrons per nuclear capture
    double _photonMean; // mean number of ejected photons per nuclear capture
    double _protonElow; // Range of proton energy
    double _protonEhi;
    double _neutronElow; //Range of neutron energy
    double _neutronEhi;
    double _photonElow; //Range of photon energy
    double _photonEhi;
    double _czmin; //Range of cos(polar angle) of ejection
    double _czmax;
    double _phimin; //Range of azimuth
    double _phimax;
    bool _PStoDSDelay;
    bool _pPulseDelay;
    double _pPulseShift;

    // Activate the folding procedure on generation time. Default is on
    bool _timeFolding;

    int _nProtonBins; //number of bins for proton energy spectrum
    int _nNeutronBins; //number of bins for neutrons energy spectrum
    int _nPhotonBins; //number of bins for photon energy spectrum

    std::string _STfname;
    int _nToSkip;

    // Class object to generate position and time of the particle
    std::unique_ptr<FoilParticleGenerator> _fGenerator;

    double _pMass, _nMass; //Particle masses

    double _p; //Particle momentum

    // Limits on the generated time.
    double _tmin;
    double _tmax;

    // Histogram control.
    bool _doHistograms;

    // end: parameters that can be configured from the config file.

    //Random generators
    CLHEP::RandPoissonQ _randPoissonQ, _randPoissonP, _randPoissonN, _randPoissonG;
    RandomUnitSphere    _randomUnitSphere;
    CLHEP::RandGeneral _shapeP, _shapeN, _shapeG;

    const DetectorSystem *_detSys;

    TH1D* _hNuclearCaptureMultiplicity;
    TH1D* _hProtonMultiplicity;
    TH1D* _hNeutronMultiplicity;
    TH1D* _hPhotonMultiplicity;
    TH1D* _hProtonKE;
    TH1D* _hNeutronKE;
    TH1D* _hPhotonKE;
    TH1D* _hProtonKEZoom;
    TH1D* _hNeutronKEZoom;
    TH1D* _hPhotonKEZoom;
    TH1D* _hProtonMomentumMeV;
    TH1D* _hNeutronMomentumMeV;
    TH1D* _hPhotonMomentumMeV;
    TH1D* _hProtonzPosition;
    TH1D* _hNeutronzPosition;
    TH1D* _hPhotonzPosition;
    TH1D* _hProtonCz;
    TH1D* _hNeutronCz;
    TH1D* _hPhotonCz;
    TH1D* _hProtonPhi;
    TH1D* _hNeutronPhi;
    TH1D* _hPhotonPhi;
    TH1D* _hProtonTime;
    TH1D* _hNeutronTime;
    TH1D* _hPhotonTime;
    TH1D* _hProtonMudelay;
    TH1D* _hProtonPulsedelay;
    TH1D* _hNeutronMudelay;
    TH1D* _hNeutronPulsedelay;
    TH1D* _hPhotonMudelay;
    TH1D* _hPhotonPulsedelay;

    //Functions used to calculate the energy spectrum of the proton
    std::vector<double> binnedEnergySpectrumProton();
    std::vector<double> binnedEnergySpectrumNeutron();
    std::vector<double> binnedEnergySpectrumPhoton();
    double protonEnergySpectrum( double e );
    int evaluateNeutronBins();

  };

} // end namespace mu2e,

#endif /* EventGenerator_NuclearCaptureGun_hh */
