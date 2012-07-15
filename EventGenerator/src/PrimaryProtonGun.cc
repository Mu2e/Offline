//
// Generate one proton with the primary proton beam energy and
// incident on the upstream face of the production target.
// See the header file for details.
//
// $Id: PrimaryProtonGun.cc,v 1.20 2012/07/15 22:06:17 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:06:17 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Framework includes
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "EventGenerator/inc/PrimaryProtonGun.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "TargetGeom/inc/Target.hh"
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"

// CLHEP includes.
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

//ROOT Includes
#include "TH1D.h"

using namespace std;

namespace mu2e {

  // The kinetic energy of proton is 8 Gev; E = T + m, p = sqrt(E^2 - m^2)
  static const double pBeam = 8888.6*CLHEP::MeV;

  PrimaryProtonGun::PrimaryProtonGun( art::Run& run, const SimpleConfig& config ):

    // The base class;
    GeneratorBase(),

    _proton_mass(GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::p_plus).ref().mass().value()),

    // Parameters from the run time configuration.
    _p(config.getDouble("primaryProtonGun.p", pBeam )),
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
    _doHistograms(config.getBool("primaryProtonGun.doHistograms", true)){

    if ( _doHistograms ){
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "PrimaryProtonGun" );
      _hmomentum = tfdir.make<TH1D>( "hmomentum", "Primary Proton Momentum, MeV",        10, 7000., 9000.);
      _hKE       = tfdir.make<TH1D>( "hKE",       "Primary Proton Kinetic Energy, MeV", 200, 7000., 9000.);
    }

  }

  PrimaryProtonGun::~PrimaryProtonGun(){
  }

  void PrimaryProtonGun::generate( GenParticleCollection& genParts ){

    // For all distributions, use the engine managed by the RandomNumberGenerator.
    static CLHEP::RandFlat   randFlat        ( getEngine() );
    static CLHEP::RandGaussQ randGaussQ      ( getEngine(), 0., _beamSpotSigma );
    static RandomUnitSphere  randomUnitSphere( getEngine(), _czmin, _czmax, _phimin, _phimax);

    // Simulate the size of the beam spot.
    double dx = 0;
    double dy = 0;
    double dr = 0;
    double phi = 0;

    if (_shape == std::string("gaus")) {
      dx = randGaussQ.fire();
      dy = randGaussQ.fire();
    }
    else if (_shape == std::string("flat")) {
      // even in the circle
      dr = _rmax*sqrt( randFlat.fire() ); 
      phi = 2.*M_PI*randFlat.fire(); 
      dx = dr*cos(phi);
      dy = dr*sin(phi);
    }

    // Generated position.
    CLHEP::Hep3Vector pos( _beamDisplacementOnTarget.x() + dx,
                           _beamDisplacementOnTarget.y() + dy,
                           _beamDisplacementOnTarget.z() );

    // Distribution of generation time is flat.  This needs to be improved.
    const double time = _tmin + (_tmax-_tmin)*randFlat.fire();

    // Energy and kinetic energy.
    double e = sqrt( _p*_p + _proton_mass*_proton_mass);
    double ekine = e - _proton_mass;

    // Generated 4 momentum.
    CLHEP::HepLorentzVector mom(  randomUnitSphere.fire(_p), e );

    CLHEP::HepRotation rot;
    rot.setTheta( _beamRotationTheta * CLHEP::degree );
    rot.setPhi( _beamRotationPhi * CLHEP::degree );
    rot.setPsi( _beamRotationPsi * CLHEP::degree );
    mom = rot * mom;
 
    // Add the proton to the list of generated particles.
    genParts.push_back( GenParticle( PDGCode::p_plus, GenId::primaryProtonGun, pos, mom, time));

    if ( _doHistograms ){
      _hKE->Fill(ekine);
      _hmomentum->Fill(_p);
    }

  }

}
