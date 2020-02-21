#ifndef EventGenerator_PrimaryProtonGunImpl_hh
#define EventGenerator_PrimaryProtonGunImpl_hh
//
// Generate a proton with the primary proton energy
//
// $Id: PrimaryProtonGunImpl.hh,v 1.17 2013/12/13 21:35:07 gandr Exp $
// $Author: gandr $
// $Date: 2013/12/13 21:35:07 $
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
#include "EventGenerator/inc/PrimaryProtonGunConfig.hh"

// Framework Includes

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Units/PhysicalConstants.h"


// Forward references outside of namespace mu2e
namespace art {
  class Run;
}

namespace mu2e {
    
  class PhysicsParams;
    
  class PrimaryProtonGunImpl: public GeneratorBase {
      
  public:
      
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
