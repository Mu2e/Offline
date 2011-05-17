#ifndef PARTICLEGUN_HH
#define PARTICLEGUN_HH
//
// Shoots a single particle gun and puts its output into a generated event.
//
// $Id: ParticleGun.hh,v 1.5 2011/05/17 15:35:59 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:35:59 $
//
// Original author Rob Kutschke
//
// The position is given in the Mu2e coordinate system.
// 

// Mu2e includes
#include "EventGenerator/inc/GeneratorBase.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"

// External includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"

// Forward references.
namespace art{
  class Run;
}
class TH1F;

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class ParticleGun: public GeneratorBase{

  public:
    ParticleGun( art::Run const& run, const SimpleConfig& config );
    virtual ~ParticleGun();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    // Start: Information from the run time configuration.

    // Number of particles per event.
    // If positive, mean of a Poisson distribution.
    // If negative, then exactly that number of particles per event.
    double _mean;

    // PDG particle id code of the particle to be generated.
    PDGCode::type _pdgId;

    // Angular range over which particles will be generated.
    double _czmin;
    double _czmax;
    double _phimin;
    double _phimax;

    // Momentum range of the particle.  Units are MeV.
    double _pmin;
    double _pmax;

    // Time range over which the particle will be produced; units are ns.
    double _tmin;
    double _tmax;

    // Particles will be produced in a box, specified by 
    // a point in the Tracker coordinate system and 
    // the half lengths of the box.  Units are mm.
    // The point (0,0,0) is at the origin of the Mu2e coordinate system.
    CLHEP::Hep3Vector _point;
    CLHEP::Hep3Vector _halfLength;

    // Enable histograms
    bool _doHistograms;

    // End: Information from the run time configuration.

    // Random number distributions.
    CLHEP::RandFlat     _randFlat;
    CLHEP::RandPoissonQ _randPoissonQ;
    RandomUnitSphere    _randomUnitSphere;

    // Derived information.

    // Mass of the particle to be generated.
    double _mass;

    // Ranges of momentum and time.
    double _dp, _dt;

    // Histogram information.
    TH1F* _hMultiplicity;
    TH1F* _hMomentum;
    TH1F* _hCz;
    TH1F* _hX0;
    TH1F* _hY0;
    TH1F* _hZ0;
    TH1F* _hT0;

  };

} // end namespace mu2e,

#endif
