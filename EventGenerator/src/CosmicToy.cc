//
// A really, really, stupid model of cosmic rays.
// The purpose is to provide an example of the interface.
//
// $Id: CosmicToy.cc,v 1.2 2009/10/16 04:20:52 shanahan Exp $
// $Author: shanahan $
// $Date: 2009/10/16 04:20:52 $
//
// Original author Rob Kutschke
//
#include <iostream>

#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"

#include "EventGenerator/inc/CosmicToy.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TargetGeom/inc/Target.hh"

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

  // Need a home for these.  Can we use the root version?
  static const int  pdg_electron = 11;
  static const int  pdg_muon     = 13;
  static const double m = 0.1056584;

  CosmicToy::CosmicToy( edm::Run& run, const SimpleConfig& config ):
    GeneratorBase(){

    _mean = config.getDouble("cosmictoy.mean",0.);
    cout << "Cosmic mean: " << endl;

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
    else          n = RandPoisson::shoot(_mean);

    _cosmicMultiplicity->Fill(n);

    for ( int i=0; i<n; ++i ){

      // We are after the high energy tail.  Units are MeV.
      double p = 1500.;

      // Look at a small angle around the zenith, going downard.
      double theta = 0.1*RandFlat::shoot();

      // Cosine and sin of polar angle wrt y axis.
      double cy = cos(theta);
      double sy = sin(theta);

      double phi = 2.*M_PI*RandFlat::shoot();
     
      double e = sqrt(p*p +m*m);
      HepLorentzVector mom(p*sy*cos(phi), -p*cy, p*sy*sin(phi), e);

      // Toy footprint is 5m in diameter.
      double range = 5000.;

      // Position in a reference plane that is just above the ground.
      // We can worry later about the exact meaning of the height.
      // The G4 interface code ( PrimaryGeneratorAction) will put it
      // at the right height.
      double x = (1.-2.*RandFlat::shoot())*range;
      double y = 0.;
      double z = (1.-2.*RandFlat::shoot())*range;
      Hep3Vector pos( x, y, z );

      double time = 1694*+RandFlat::shoot();

      // Pick a random charge.
      int pid = (RandFlat::shoot() >0.5) ? pdg_muon : -pdg_muon;

      // Add the cosmic to  the list.
      genParts.push_back( ToyGenParticle( pid, GenId::cosmicToy, pos, mom, time));
    }
  }

}
