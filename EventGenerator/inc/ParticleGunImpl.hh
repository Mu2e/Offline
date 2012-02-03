#ifndef EventGenerator_ParticleGunImpl_hh
#define EventGenerator_ParticleGunImpl_hh
//
// Shoots a single particle gun and puts its output into a generated event.
// This class implements the common code for "particle gun-like" generators.
//
// $Id: ParticleGunImpl.hh,v 1.1 2012/02/03 05:08:06 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/03 05:08:06 $
//
// Original author Rob Kutschke, re-factored for use in multiple generators by Andrei Gaponenko.
//
// The position is given in the Mu2e coordinate system.
//

// Mu2e includes
#include "EventGenerator/inc/GeneratorBase.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <string>

class TH1F;

namespace mu2e {

  // Forward reference.
  class RandomUnitSphereParams;

  class ParticleGunImpl: public GeneratorBase{
    
  public:
    ParticleGunImpl(double meanMultiplicity,
		    PDGCode::type pdgId,
		    double pmin,
		    double pmax,
		    const RandomUnitSphereParams& angles,
		    double tmin,
		    double tmax,

		    const CLHEP::Hep3Vector& point,
		    const CLHEP::Hep3Vector& halfLength,
		    
		    // empty string "" means don't histogram
		    const std::string& histoDir,

		    bool verbose);

    virtual void generate( GenParticleCollection&  );
    
  private:
    // Disallow copying because of the owned pointer members
    ParticleGunImpl(const ParticleGunImpl&);
    ParticleGunImpl& operator=(const ParticleGunImpl&);

    // Start: Information from the run time configuration.

    // Number of particles per event.
    // If positive, mean of a Poisson distribution.
    // If negative, then exactly that number of particles per event.
    double _mean;

    // PDG particle id code of the particle to be generated.
    PDGCode::type _pdgId;

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

    // enable output
    bool _verbose;
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

#endif /* EventGenerator_ParticleGunImpl_hh */
