#ifndef EventGenerator_CosmicDYB_hh
#define EventGenerator_CosmicDYB_hh
//
// Muon generator, uses Daya Bay libraries
//
// $Id: CosmicDYB.hh,v 1.13 2014/03/22 21:40:43 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/03/22 21:40:43 $
//
// Original author Yury Kolomensky
//

#include <vector>

#include "EventGenerator/inc/GeneratorBase.hh"
#include "EventGenerator/inc/DYBGenerator.hh"

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
  };  // CosmicDYB

}  // namespace mu2e

#endif /* EventGenerator_CosmicDYB_hh */
