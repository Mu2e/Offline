//
// Generate one proton with the primary proton beam energy and
// incident on the upstream face of the production target.
// See the header file for details.
//
// $Id: PrimaryProtonGun.cc,v 1.24 2013/12/13 21:35:07 gandr Exp $
// $Author: gandr $
// $Date: 2013/12/13 21:35:07 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Framework includes
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "EventGenerator/inc/PrimaryProtonGun.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"

// CLHEP includes.
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

//ROOT Includes
#include "TH1D.h"

using namespace std;

namespace mu2e {

  PrimaryProtonGun::PrimaryProtonGun(CLHEP::HepRandomEngine& engine, art::Run const& run, SimpleConfig const& config):

    _gunRotation(GeomHandle<ProductionTarget>()->protonBeamRotation()),
    _gunOrigin(GeomHandle<ProductionTarget>()->position()
               + _gunRotation*CLHEP::Hep3Vector(0., 0., GeomHandle<ProductionTarget>()->halfLength())),

    _proton_mass(GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::p_plus).ref().mass().value()),

    // Parameters from the run time configuration.
    _p(config.getDouble("primaryProtonGun.p", GlobalConstantsHandle<PhysicsParams>()->getProtonMomentum() )),
    _beamDisplacementOnTarget(config.getHep3Vector("beamDisplacementOnTarget")),
    _beamRotationTheta(config.getDouble("beamRotationTheta", 0)),
    _beamRotationPhi(config.getDouble("beamRotationPhi", 0)),
    _beamRotationPsi(config.getDouble("beamRotationPsi", 0)),
    _beamSpotSigma(config.getDouble("primaryProtonGun.beamSpotSigma")),
    _czmin(config.getDouble("primaryProtonGun.czmin", -1.)),
    _czmax(config.getDouble("primaryProtonGun.czmax", -1.)),
    _phimin(config.getDouble("primaryProtonGun.phimin", 0. )),
    _phimax(config.getDouble("primaryProtonGun.phimax", CLHEP::twopi)),
    _tmin(config.getDouble("primaryProtonGun.tmin",   0.)),
    _tmax(config.getDouble("primaryProtonGun.tmax", 100.)),
    _shape(config.getString("primaryProtonGun.shape", "gaus")),
    _rmax(config.getDouble("primaryProtonGun.rmax", 100.)),
    _mean(config.getDouble("primaryProtonGun.mean", -1.)),
    // For all distributions, use the engine managed by the RandomNumberGenerator.
    _randPoissonQ{engine, std::abs(_mean)},
    _randFlat{engine},
    _randGaussQ{engine, 0., _beamSpotSigma},
    _randomUnitSphere{engine, _czmin, _czmax, _phimin, _phimax}
    {}

  void PrimaryProtonGun::generate( GenParticleCollection& genParts ){
    long n = _mean < 0 ? static_cast<long>(-_mean): _randPoissonQ.fire();
    for ( int j =0; j<n; ++j ){
      generateOne(genParts);
    }
  }

  void PrimaryProtonGun::generateOne( GenParticleCollection& genParts ){

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

  }

}
