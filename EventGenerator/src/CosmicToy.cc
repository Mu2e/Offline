//
// A really, really, stupid model of cosmic rays.
// The purpose is to provide an example of the interface.
//
// $Id: CosmicToy.cc,v 1.12 2011/05/18 22:01:46 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 22:01:46 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Framework includes.
#include "art/Framework/Core/Run.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "EventGenerator/inc/CosmicToy.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"
#include "TargetGeom/inc/Target.hh"

// CLHEP includes
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

// ROOT includes
#include "TH1D.h"

using namespace std;

namespace mu2e {

  // Mass of the muon, in GeV.
  // Once we have the HepPDT package installed, get the mass from there.
  static const double m = 0.1056584;

  CosmicToy::CosmicToy( art::Run& run, const SimpleConfig& config ):

    // Base class.
    GeneratorBase(),

    // From run time configuration.
    _mean(config.getDouble("cosmictoy.mean",2.)),
    _doHistograms(config.getBool("cosmictoy.doHistograms",true)),

    // Histograms.
    _hMultiplicity(0),
    _hMomentum(0),
    _hAngle(0),

    // Random number distributions; getEngine comes from the base class.
    _randFlat( getEngine() ),
    _randPoissonQ( getEngine(), std::abs(_mean) ){

    // Sanity check.
    if ( std::abs(_mean) > 99999. ) {
      throw cet::exception("RANGE")
        << "CosmicToy has been asked to produce a crazily large number of electrons."
        << _mean
        << "\n";
    }

    // Access conditions data.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    ConditionsHandle<DAQParams>         daqPar("ignored");

    // Start time for generation is a little before the start time
    // of the DAQ system.
    double offset = 100.;

    // Start and end times for generation.
    _tmin = (daqPar->t0 > offset)? daqPar->t0-offset : 0.;
    _tmax = accPar->deBuncherPeriod;
    _dt   = _tmax - _tmin;

    // Book histograms.
    if ( _doHistograms ) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "CosmicToy" );
      _hMultiplicity = tfdir.make<TH1D>( "hMultiplicity", "Cosmic Toy Multiplicity", 20, 0, 20);
      _hMomentum     = tfdir.make<TH1D>( "hMomentum",     "Cosmic Toy Momentum",   100, 0, 10000.);
      _hAngle        = tfdir.make<TH1D>( "hAngle",        "Cosmic Toy Angle from Zenith", 100, 0, 0.5);
    }

  }

  CosmicToy::~CosmicToy(){
  }

  void CosmicToy::generate( GenParticleCollection& genParts ){

    // Pick a number of muons from a Poisson distribution.
    // Choose the number of electrons to generate this event.
    long n = _mean < 0 ? static_cast<long>(-_mean) : _randPoissonQ.fire();
    if ( _doHistograms ){
      _hMultiplicity->Fill(n);
    }

    for ( int i=0; i<n; ++i ){

      // We are after the high energy tail.  Units are MeV.
      double p = 1500.;

      // Look at a small angle around the zenith, going downard.
      double theta = 0.1*_randFlat.fire();

      // Cosine and sin of polar angle wrt y axis.
      double cy = cos(theta);
      double sy = sin(theta);

      double phi = 2.*M_PI*_randFlat.fire();

      double e = sqrt(p*p +m*m);
      CLHEP::HepLorentzVector mom(p*sy*cos(phi), -p*cy, p*sy*sin(phi), e);

      // Footprint of this toy model is 10m on each side.
      double halfLength = 5000.;

      // Position in a reference plane that is just above the ground.
      // We can worry later about the exact meaning of the height.
      // The G4 interface code ( PrimaryGeneratorAction) will put it
      // at the right height.
      double x = (1.-2.*_randFlat.fire())*halfLength;
      double y = 0.;
      double z = (1.-2.*_randFlat.fire())*halfLength;
      CLHEP::Hep3Vector pos( x, y, z );

      double time = _tmin + _dt*_randFlat.fire();

      // Pick a random charge.
      PDGCode::type pid = (_randFlat.fire() >0.5) ? PDGCode::mu_minus : PDGCode::mu_plus;

      // Add the cosmic ray muon to the list of generated particles.
      genParts.push_back( GenParticle( pid, GenId::cosmicToy, pos, mom, time));

      if ( _doHistograms ){
        static CLHEP::Hep3Vector vertical( 0., -1., 0.);
        double angle = mom.angle(vertical);

        _hMomentum->Fill(p);
        _hAngle->Fill(angle);
      }

    }
  }

}
