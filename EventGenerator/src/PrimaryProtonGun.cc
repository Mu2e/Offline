//
// Generate one proton with the primary proton beam energy and
// incident on the upstream face of the production target.
// See the header file for details.
//
// $Id: PrimaryProtonGun.cc,v 1.12 2011/05/18 16:11:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 16:11:17 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Framework includes
#include "art/Framework/Core/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"

// Mu2e includes
#include "EventGenerator/inc/PrimaryProtonGun.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "TargetGeom/inc/Target.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

// CLHEP includes.
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/RandGaussQ.h"

//ROOT Includes
#include "TH1D.h"

using namespace std;

namespace mu2e {

  // Mass of the proton.
  // Once we have the HepPDT package installed, get this number from there.
  static const double m = 938.272*CLHEP::MeV;

  // The kinetic energy of proton is 8 Gev; E = T + m, p = sqrt(E^2 - m^2)
  static const double pBeam = 8888.6*CLHEP::MeV;

  PrimaryProtonGun::PrimaryProtonGun( art::Run& run, const SimpleConfig& config ):

    // The base class;
    GeneratorBase(),

    // Parameters from the run time configuration.
    _p(config.getDouble("primaryProtonGun.p", pBeam )),
    _beamDisplacementOnTarget(config.getHep3Vector("beamDisplacementOnTarget")),
    _beamSpotSigma(config.getDouble("primaryProtonGun.beamSpotSigma")),
    _czmin(config.getDouble("primaryProtonGun.czmin", -1.)),
    _czmax(config.getDouble("primaryProtonGun.czmax", -1.)),
    _phimin(config.getDouble("primaryProtonGun.phimin", 0. )),
    _phimax(config.getDouble("primaryProtonGun.phimax", CLHEP::twopi)),
    _tmin(config.getDouble("primaryProtonGun.tmin",   0.)),
    _tmax(config.getDouble("primaryProtonGun.tmax", 100.)),
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

  void PrimaryProtonGun::generate( ToyGenParticleCollection& genParts ){

    // For all distributions, use the engine managed by the RandomNumberGeneratorService.
    static CLHEP::RandFlat   randFlat        ( getEngine() );
    static CLHEP::RandGaussQ randGaussQ      ( getEngine(), 0., _beamSpotSigma );
    static RandomUnitSphere  randomUnitSphere( getEngine(), _czmin, _czmax, _phimin, _phimax);

    // Simulate the size of the beam spot.
    double dx = randGaussQ.fire();
    double dy = randGaussQ.fire();

    // Generated position.
    CLHEP::Hep3Vector pos( _beamDisplacementOnTarget.x() + dx,
                           _beamDisplacementOnTarget.y() + dy,
                           _beamDisplacementOnTarget.z() );

    // Distribution of generation time is flat.  This needs to be improved.
    const double time = _tmin + (_tmax-_tmin)*randFlat.fire();

    // Energy and kinetic energy.
    double e = sqrt( _p*_p +m*m);
    double ekine = e - m;

    // Generated 4 momentum.
    CLHEP::HepLorentzVector mom(  randomUnitSphere.fire(_p), e );

    // Add the proton to the list of generated particles.
    genParts.push_back( ToyGenParticle( PDGCode::p_plus, GenId::primaryProtonGun, pos, mom, time));

    if ( _doHistograms ){
      _hKE->Fill(ekine);
      _hmomentum->Fill(_p);
    }

  }

}
