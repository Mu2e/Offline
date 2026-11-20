#ifndef EventGenerator_CosmicDYB_hh
#define EventGenerator_CosmicDYB_hh
//
// Muon generator, uses Daya Bay libraries
//
//
// Original author Yury Kolomensky
//

#include <vector>

#include "Offline/EventGenerator/inc/GeneratorBase.hh"
#include "Offline/EventGenerator/inc/DYBGenerator.hh"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"

// Forward declarations.
class TH1D;
class TH2D;
namespace art{
  class Run;
}

namespace mu2e {

  // Forward declaratations.
  class DYBGenerator;
  class SimpleConfig;

  class CosmicDYB: public GeneratorBase {

  public:
    CosmicDYB(CLHEP::HepRandomEngine& engine, art::Run& run, const SimpleConfig& config);
    virtual ~CosmicDYB();

    virtual void generate( GenParticleCollection&  );

    // Muons per second crossing the production surfaces (Hz), from the flux parametrization.
    double rate() const { return _totalRate; }

    // Total area of the production surfaces (m^2) and the muon energy range (GeV).
    double area() const { return _area; }
    double eMin() const { return _muEMin; }
    double eMax() const { return _muEMax; }

    // Muons generated, and the live time they represent (s), since the last resetCounts().
    // Fixed multiplicity (cosmicDYB.mean < 0): muons / rate. Poisson multiplicity: events * mean / rate.
    unsigned long long nMuons() const { return _nMuons; }
    double liveTime() const;
    void resetCounts() { _nEvents = 0; _nMuons = 0; }

  private:

    int   _verbose;
    bool  _doHistograms;
    TH2D *_hStartXZ;
    TH1D *_hStartY;
    TH1D *_hStartPlane;
    TH1D *_hStartE;
    TH1D *_hStartTheta;
    TH1D *_hStartPhi;

    // Muon mass
    double _mMu;

    // Mean multiplicity. If negative, use -_mean as a fixed number
    double _mean;
    double _muEMin;   // min and max values of muon energy (MeV)
    double _muEMax;
    double _muThMin; // min and max zenith angles
    double _muThMax;
    double _muPhiMin;
    double _muPhiMax;

    // half lengths to generate events (mm)
    double _dx;
    double _dy;
    double _dz;

    // Number of lookup bins for DYBGenerator.
    long _ne;
    long _nth;

    // Time range (in ns) over which to generate events.
    double _tmin;
    double _tmax;
    double _dt;

    // Random number distributions.
    CLHEP::RandFlat     _randFlat;
    CLHEP::RandPoissonQ _randPoissonQ;

    CLHEP::Hep3Vector        _productionCenterInMu2e;
    DYBGenerator::Direction  _direction;
    bool                     _dontProjectToSurface;
    bool                     _checkedProductionPlanes;

    std::vector<boost::shared_ptr<DYBGenerator> > _generators;
    std::vector<double>                           _boxFraction;

    double             _totalRate = 0;  // Hz
    double             _area = 0;       // m^2
    unsigned long long _nEvents = 0;
    unsigned long long _nMuons = 0;
  };  // CosmicDYB

}  // namespace mu2e

#endif /* EventGenerator_CosmicDYB_hh */
