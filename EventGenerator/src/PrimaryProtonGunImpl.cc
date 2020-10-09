//
// Generate one proton with the primary proton beam energy and
// incident on the upstream face of the production target.
// See the header file for details.
//
//
// Original author Rob Kutschke
//

// C++ includes.

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "EventGenerator/inc/PrimaryProtonGunImpl.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"

// CLHEP includes.
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

namespace mu2e {

    PrimaryProtonGunImpl::PrimaryProtonGunImpl(CLHEP::HepRandomEngine& engine, const PrimaryProtonGunConfig& config):
    
    _gunRotation(GeomHandle<ProductionTarget>()->protonBeamRotation()),
    _gunOrigin(GeomHandle<ProductionTarget>()->targetPositionByVersion()
               + _gunRotation*CLHEP::Hep3Vector(0., 0., GeomHandle<ProductionTarget>()->targetHalfLengthByVersion())),

    _proton_mass(GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::p_plus).ref().mass().value()),

    // Parameters from the run time configuration.
    _config(config),
    _p(GlobalConstantsHandle<PhysicsParams>()->getProtonMomentum()),
    _beamDisplacementOnTarget(config.beamDisplacementOnTarget()),
    _beamRotationTheta(config.beamRotationTheta()),
    _beamRotationPhi(config.beamRotationPhi()),
    _beamRotationPsi(config.beamRotationPsi()),
    _tmin(config.tmin()),
    _tmax(config.tmax()),
    _shape(config.shape()),
    _rmax(config.rmax()),
    
    // For all distributions, use the engine managed by the RandomNumberGenerator.
    _randPoissonQ{engine, std::abs(config.mean())},
    _randFlat{engine},
    _randGaussQ{engine, 0., config.beamSpotSigma()},
    _randomUnitSphere{engine, config.czmin(), config.czmax(), config.phimin(), config.phimax()}
    {
        _config.proton_momentum(_p);
        
    }

  void PrimaryProtonGunImpl::generate( GenParticleCollection& genParts ){
    long n = _config.mean() < 0 ? static_cast<long>(-_config.mean()): _randPoissonQ.fire();
    for ( int j =0; j<n; ++j ){
      generateOne(genParts);
    }
  }

  void PrimaryProtonGunImpl::generateOne( GenParticleCollection& genParts ){

    // Simulate the size of the beam spot.
    double dx = 0;
    double dy = 0;
    double dr = 0;
    double phi = 0;
    if (_shape == std::string("gaus")) {
      dx = _randGaussQ.fire();
      dy = _randGaussQ.fire();
    }
    else if (_shape == std::string("flat")) {
      // even in the circle
      dr = _rmax * sqrt(_randFlat.fire());
      phi = 2. * M_PI * _randFlat.fire();
      dx = dr*cos(phi);
      dy = dr*sin(phi);
    }

    // Generated position.
    CLHEP::Hep3Vector pos( _beamDisplacementOnTarget.x() + dx,
                           _beamDisplacementOnTarget.y() + dy,
                           _beamDisplacementOnTarget.z() );

    // Distribution of generation time is flat.  This needs to be improved.
    const double time = _tmin + (_tmax-_tmin) * _randFlat.fire();

    // Energy and kinetic energy.
    double e = sqrt( _p*_p + _proton_mass*_proton_mass);


    // Generated 4 momentum.
    CLHEP::HepLorentzVector mom( _randomUnitSphere.fire(_p), e );

    CLHEP::HepRotation rot;
    rot.setTheta( _beamRotationTheta * CLHEP::degree );
    rot.setPhi( _beamRotationPhi * CLHEP::degree );
    rot.setPsi( _beamRotationPsi * CLHEP::degree );
    mom = rot * mom;


    // Add the proton to the list of generated particles.
    genParts.push_back( GenParticle( PDGCode::p_plus, GenId::primaryProtonGun,
                                     // Convert position to Mu2e coordinates
                                     _gunRotation*pos + _gunOrigin,
                                     // Convert momentum to Mu2e coordinates
                                     _gunRotation*mom,
                                     time));
      
  }//generateOne
    

}
