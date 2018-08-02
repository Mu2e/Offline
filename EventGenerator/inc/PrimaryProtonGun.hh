#ifndef EventGenerator_PrimaryProtonGun_hh
#define EventGenerator_PrimaryProtonGun_hh
//
// Generate a proton with the primary proton energy
//
// $Id: PrimaryProtonGun.hh,v 1.17 2013/12/13 21:35:07 gandr Exp $
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

// Framework Includes
#include "art/Framework/Principal/Run.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandPoissonQ.h"

// Forward references outside of namespace mu2e
class TH1D;
namespace art {
  class Run;
}

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class PrimaryProtonGun: public GeneratorBase{
  public:

    PrimaryProtonGun(CLHEP::HepRandomEngine& engine, art::Run& run, SimpleConfig const& config);
    ~PrimaryProtonGun() = default;

    virtual void generate( GenParticleCollection&  );

  private:
    CLHEP::HepRotation _gunRotation; // rotates target frame to Mu2e frame
    CLHEP::Hep3Vector _gunOrigin;

    double _proton_mass;

    // Start parameters from the run time configuration.

    // Momentum of the generated proton; in MeV.
    double _p;

    // Offset of production point relative to the origin described above; in mm.
    CLHEP::Hep3Vector _beamDisplacementOnTarget;

    // Rotation of beam direction wrt to target angle; in deg.
    double _beamRotationTheta;
    double _beamRotationPhi;
    double _beamRotationPsi;

    // Beamspot is a 2D gaussian with this sigma in both x and y.
    double _beamSpotSigma;

    // Limits on the generated direction; generated over a unit sphere within these limits.
    double _czmin;
    double _czmax;
    double _phimin;
    double _phimax;

    // Time of generation is a flat distribution within these limits. Time in ns.
    double _tmin;
    double _tmax;

    // Shape
    std::string _shape;

    // radius max, for flat distribution
    double _rmax;

    double _mean;  // poisson mean; negative for non-random abs(mean)
    CLHEP::RandPoissonQ _randPoissonQ;
    CLHEP::RandFlat     _randFlat;
    CLHEP::RandGaussQ   _randGaussQ;
    RandomUnitSphere    _randomUnitSphere;


    // Make histograms or not.
    bool _doHistograms;

    // End parameters from the run time configuration.

    // Histograms.
    TH1D* _hKE;
    TH1D* _hKEZoom;
    TH1D* _hmomentum;
    TH1D* _hposx;
    TH1D* _hposy;
    TH1D* _hposz;
    TH1D* _hcosTheta;
    TH1D* _htime;

    virtual void generateOne( GenParticleCollection&  );
  };

} // end namespace mu2e,

#endif /* EventGenerator_PrimaryProtonGun_hh */
