#ifndef EventGenerator_PrimaryProtonGunImpl_hh
#define EventGenerator_PrimaryProtonGunImpl_hh
//
// Generate a proton with the primary proton energy
//
//
//
// The coordinate system used in this class is:
//
// 1) The origin is at the upstream face of the production target.
//    (upstream in the sense of the proton beam).
//
// 2) The positive z direction is along the cylinder axis of the production target,
//    in the usual Mu2e sense (from PS toward DS); ie the ideal beam direction
//    is in the -z direction.  The target is slightly tilted from the Mu2e z axis.
//    The x and y axes are in the plane of the upstream face of the production target;
//    the sense of these axes is defined by the rotation matrix that was used to position
//    the production target inside the PS vacuum volume; see the section titled
//    "Production Target" inside Mu2e-doc-938.
//
// 3) The transformation from this coordinate system to the G4 world coordinate
//    system is done in Mu2eG4/src/PrimaryGeneratorAction.cc .
//


// Mu2e includes
#include "EventGenerator/inc/GeneratorBase.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"

// Framework Includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/TupleAs.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Units/PhysicalConstants.h"

// C++ includes
#include <string>


// Forward references outside of namespace mu2e
namespace art {
  class Run;
}

namespace mu2e {
        
  class PrimaryProtonGunImpl: public GeneratorBase {
      
  public:
      
    struct PrimaryProtonGunConfig {
          
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      using Hep3Vector_t = CLHEP::Hep3Vector(double, double, double);
          
      fhicl::OptionalAtom<double> proton_momentum {Name("proton_momentum"), Comment("Momentum of generated proton in MeV. \nThe default is set dynamically in the constructor using the Global Constants Service, PhysicsParams.")};
          
      fhicl::TupleAs<Hep3Vector_t> beamDisplacementOnTarget {Name("beamDisplacementOnTarget"), Comment("Offset of production point relative to upstream face of production target; in mm.") };
          
      fhicl::Atom<double> beamRotationTheta {Name("beamRotationTheta"), Comment("Rotation of beam direction in Theta wrt to target angle; in deg"), 0.};
      fhicl::Atom<double> beamRotationPhi {Name("beamRotationPhi"), Comment("Rotation of beam direction in Phi wrt to target angle; in deg"), 0.};
      fhicl::Atom<double> beamRotationPsi {Name("beamRotationPsi"), Comment("Rotation of beam direction in Psi wrt to target angle; in deg"), 0.};
          
      fhicl::Atom<double> beamSpotSigma {Name("beamSpotSigma"), Comment("Beamspot is 2D gaussian with this sigma in x and y.") };
          
      fhicl::Atom<double> czmin {Name("czmin"), Comment("Limit on generated z over a unit sphere."), -1.};
      fhicl::Atom<double> czmax {Name("czmax"), Comment("Limit on generated z over a unit sphere."), 1.};
          
      fhicl::Atom<double> phimin {Name("phimin"), Comment("Limit on generated phi over a unit sphere."), 0.};
      fhicl::Atom<double> phimax {Name("phimax"), Comment("Limit on generated phi over a unit sphere."), CLHEP::twopi};
          
      fhicl::Atom<double> tmin {Name("tmin"), Comment("Time of generation is flat within these limits. Time in ns."), 0.};
      fhicl::Atom<double> tmax {Name("tmax"), Comment("Time of generation is flat within these limits. Time in ns."), 100.};
          
      fhicl::Atom<std::string> shape {Name("shape"), Comment("Shape of beam spot as a function of x and y."), "gaus"};
          
      fhicl::Atom<double> rmax {Name("rmax"), Comment("Maximum radius of beam spot for flat distribution"), 100.};
      fhicl::Atom<double> mean {Name("mean"), Comment("Poisson mean; negative for non-random abs(mean)"), -1};
    };
      
    explicit PrimaryProtonGunImpl(CLHEP::HepRandomEngine& engine, const PrimaryProtonGunConfig& config);
      
    virtual void generate( GenParticleCollection& );

  private:
      
    CLHEP::HepRotation _gunRotation; // rotates target frame to Mu2e frame
    CLHEP::Hep3Vector _gunOrigin;

    double _proton_mass;

    // Start parameters from the run time configuration.
    PrimaryProtonGunConfig _config;
      
    // Momentum of the generated proton; in MeV.
    double _p;

    // Offset of production point relative to the origin described above; in mm.
    CLHEP::Hep3Vector _beamDisplacementOnTarget;

    // Rotation of beam direction wrt to target angle; in deg.
    double _beamRotationTheta;
    double _beamRotationPhi;
    double _beamRotationPsi;

    // Time of generation is a flat distribution within these limits. Time in ns.
    double _tmin;
    double _tmax;

    // Shape
    std::string _shape;

    // radius max, for flat distribution
    double _rmax;
            
    CLHEP::RandPoissonQ _randPoissonQ;
    CLHEP::RandFlat     _randFlat;
    CLHEP::RandGaussQ   _randGaussQ;
    RandomUnitSphere    _randomUnitSphere;
      
    virtual void generateOne( GenParticleCollection& );
      
  };

} // end namespace mu2e,

#endif /* EventGenerator_PrimaryProtonGunImpl_hh */
