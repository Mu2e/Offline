//
//
// A really, really, stupid model of photons from pi-
// capture on the targets.
//
// $Id: PiCapture.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
// 

#include <iostream>

#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"

#include "EventGenerator/inc/PiCapture.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"
#include "Mu2eUtilities/inc/sqrtOrThrow.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TargetGeom/inc/Target.hh"

// From CLHEP
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

// From ROOT
#include "TH1D.h"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
using CLHEP::RandFlat;
using CLHEP::RandPoisson;


namespace mu2e {

  // Need a home for these.  Can we use the root version?
  static const int  pdg_electron = 11;
  static const int  pdg_muon     = 13;
  static const int  pdg_gamma    = 22;

  // Also need a home for this - the cycle time of the debuncher.
  static const double tcycle = 1694.;

  PiCapture::PiCapture( edm::Run& run, const SimpleConfig& config ):
    GeneratorBase(),
    _piCaptureMultiplicity(0),
    _randomUnitSphere(){

    _mean = config.getDouble("picapture.mean",0.);
    cout << "Pi capture mean: " << _mean << endl;

    // Book histograms.
    edm::Service<edm::TFileService> tfs;
    _piCaptureMultiplicity = tfs->make<TH1D>( "piCaptureMultiplicity", "Cosmic Multiplicity", 20, 0, 20);

  }

  PiCapture::~PiCapture(){
  }


  void PiCapture::generate( ToyGenParticleCollection& genParticles ){

    if ( _mean <= 0. ) return;

    // Get access to the geometry system.
    GeomHandle<Target> target;

    int nFoils = target->nFoils();

    // Pick a number of muons from a Poisson distribution.
    // This number is not meant to be real - its just to exercise the system.
    long n = RandPoisson::shoot(_mean);

    _piCaptureMultiplicity->Fill(n);

    for ( int i=0; i<n; ++i ){

      // Pick a foil.
      int ifoil = static_cast<int>(nFoils*RandFlat::shoot());
      TargetFoil const& foil = target->foil(ifoil);

      // Foil properties.
      CLHEP::Hep3Vector const& center = foil.center();
      const double r1 = foil.rIn();
      const double dr = foil.rOut() - r1;
      
      // A random point within the foil.
      const double r   = r1 + dr*RandFlat::shoot();
      const double dz  = 2.*(0.5-RandFlat::shoot())*foil.halfThickness();
      const double phi = 2.*M_PI*RandFlat::shoot();
      Hep3Vector pos( center.x()+r*cos(phi), 
		      center.y()+r*sin(phi), 
		      center.z()+dz );

      // Pick a random energy.  Replace this with a real distribution.
      const double e = 50. + 100.0*RandFlat::shoot();

      // Momentum 3 vector: uniform on a unit sphere.
      Hep3Vector p = e*_randomUnitSphere.shoot();

      // This should reflect the model of the tails of the beam
      // and the out-of-time protons.
      const double time = 1694*RandFlat::shoot();

      // Make the 4 vector.
      HepLorentzVector mom( p.x(), p.y(), p.z(), e);
      
      // Add the electron to  the list.
      genParticles.push_back( ToyGenParticle( pdg_gamma, GenId::pionCapture, pos, mom, time));
    }

  }

}
