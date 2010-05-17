//
// A really, really, stupid model of cosmic rays.
// The purpose is to provide an example of the interface.
//
// $Id: CosmicToy.cc,v 1.5 2010/05/17 21:47:33 genser Exp $
// $Author: genser $
// $Date: 2010/05/17 21:47:33 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Framework includes.
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"

// Mu2e includes.
#include "EventGenerator/inc/CosmicToy.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TargetGeom/inc/Target.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/DAQParams.hh"

// From CLHEP
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

// From Root.
#include "TH1D.h"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
using CLHEP::RandPoisson;
using CLHEP::RandFlat;


namespace mu2e {

  // Mass of the muon, in GeV.
  // Once we have the HepPDT package installed, get the mass from there.
  static const double m = 0.1056584;

  CosmicToy::CosmicToy( edm::Run& run, const SimpleConfig& config ):
    GeneratorBase(){

    _mean = config.getDouble("cosmictoy.mean",0.);

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
    edm::Service<edm::TFileService> tfs;
    _cosmicMultiplicity = tfs->make<TH1D>( "cosmicMultiplicity", "Cosmic Multiplicity", 20, 0, 20);

  }

  CosmicToy::~CosmicToy(){
  }

  void CosmicToy::generate( ToyGenParticleCollection& genParts ){

    if ( _mean <= -99999. ) return;

    // Should get the numbers here from the config file or from the
    // geometry manager.

    // Pick a number of muons from a Poisson distribution.
    long n;
    if (_mean<=0) n=(long)-_mean;
    else          n = CLHEP::RandPoisson::shoot(_mean);

    _cosmicMultiplicity->Fill(n);

    for ( int i=0; i<n; ++i ){

      // We are after the high energy tail.  Units are MeV.
      double p = 1500.;

      // Look at a small angle around the zenith, going downard.
      double theta = 0.1*CLHEP::RandFlat::shoot();

      // Cosine and sin of polar angle wrt y axis.
      double cy = cos(theta);
      double sy = sin(theta);

      double phi = 2.*M_PI*CLHEP::RandFlat::shoot();
     
      double e = sqrt(p*p +m*m);
      CLHEP::HepLorentzVector mom(p*sy*cos(phi), -p*cy, p*sy*sin(phi), e);

      // Footprint of this toy model is 10m on each side.
      double halfLength = 5000.;

      // Position in a reference plane that is just above the ground.
      // We can worry later about the exact meaning of the height.
      // The G4 interface code ( PrimaryGeneratorAction) will put it
      // at the right height.
      double x = (1.-2.*CLHEP::RandFlat::shoot())*halfLength;
      double y = 0.;
      double z = (1.-2.*CLHEP::RandFlat::shoot())*halfLength;
      CLHEP::Hep3Vector pos( x, y, z );

      double time = _tmin + _dt*CLHEP::RandFlat::shoot();

      // Pick a random charge.
      PDGCode::type pid = (CLHEP::RandFlat::shoot() >0.5) ? PDGCode::mu_minus : PDGCode::mu_plus;

      // Add the cosmic to  the list.
      genParts.push_back( ToyGenParticle( pid, GenId::cosmicToy, pos, mom, time));
    }
  }

}
