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
  class SimpleConfig;

  class CosmicDYB: public GeneratorBase {

  public:
    CosmicDYB( art::Run& run, const SimpleConfig& config );
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
    double _muCosThMin; // min and max zenith angles
    double _muCosThMax;

    // half area to generate events (cm)
    double _dx;
    double _dy;
    double _dz;
    double _y0;  // reference depth

    // Dimensions of the 2d working space for hrndg2.
    long _ne;
    long _nth;

    // Time range (in ns) over which to generate events.
    double _tmin;
    double _tmax;
    double _dt;
    double _constTime;

    // Random number distributions.
    CLHEP::RandFlat     _randFlat;
    CLHEP::RandPoissonQ _randPoissonQ;

    // Working space for hrndg2 ( working space will be on the heap).
    std::vector<double> _workingSpace;

    // production will be created only once.
    bool _createdProductionPlane;

    enum RefPointChoice {UNDEFINED, TRACKER, EXTMONFNAL, CALO, CUSTOMIZED};
    enum DirectionChoice {ALL, POSITIVE_X, NEGATIVE_X, POSITIVE_Z, NEGATIVE_Z};
    RefPointChoice    _choice;
    DirectionChoice   _directionChoice;
    CLHEP::Hep3Vector _cosmicReferencePointInMu2e;
    bool _vertical;
    bool _dontProjectToSurface;

  };  // CosmicDYB

}  // namespace mu2e

#endif /* EventGenerator_CosmicDYB_hh */
