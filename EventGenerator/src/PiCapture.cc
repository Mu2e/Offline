
// Based on Ivano Sarra's model described in mu2e Doc 665-v2
//
// $Id: PiCapture.cc,v 1.2 2009/10/16 04:20:52 shanahan Exp $
// $Author: shanahan $ 
// $Date: 2009/10/16 04:20:52 $
//
// Original author Rob Kutschke/P. Shanahan
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

  static const double emax = 138.2; // 

  PiCapture::PiCapture( edm::Run& run, const SimpleConfig& config ):
    GeneratorBase(),
    _piCaptureMultiplicity(0),
    _piCaptureEPhot(0),
    _randomUnitSphere(),
    _funcGen(0){

    _mean = config.getDouble("picapture.mean",0.);
    cout << "Pi capture mean: " << _mean << endl;
    _elow = config.getDouble("picapture.elow",38.2);
    cout << "Pi capture E Low: " << _elow << endl;
    _ehi = config.getDouble("picapture.ehi",emax);
    cout << "Pi capture E High: " << _ehi << endl;
    _nbins = config.getInt("picapture.nbins",1000);
    cout << "Pi capture N Bins: " << _nbins << endl;

    // set up the generator function
    if (_nbins>0) _bindE = (_ehi - _elow) / _nbins;
    else {
       // I'm sure this isn't the right way to do this...
       std::cout<<"Rubbish picapture.nbins = "<<_nbins<<std::endl;
       assert(0);
    }

    
    double YFunc[_nbins];
    for (int ib=0; ib<_nbins; ib++) {

       double x = _elow+(ib+0.5) * _bindE;
       YFunc[ib] = EPhotFunc(x);
    }
    _funcGen = new RandGeneral(YFunc,_nbins);

    // Book histograms.
    edm::Service<edm::TFileService> tfs;
    _piCaptureMultiplicity = tfs->make<TH1D>( "piCaptureMultiplicity", "Pi Capture Multiplicity", 20, 0, 20);
    _piCaptureEPhot = tfs->make<TH1D>( "piCaptureEPhot", "Pi Capture Photon Enegy", 200, 0, 200);
    _piCaptureEPhotZ = tfs->make<TH1D>( "piCaptureEPhotZ", "Pi Capture Photon Enegy (zoom)", 200, _elow, _ehi);

  }

  PiCapture::~PiCapture(){
  }


  void PiCapture::generate( ToyGenParticleCollection& genParticles ){

    // return for absurdly high fixed rate
    if (_mean<-99999) return;

    // Get access to the geometry system.
    GeomHandle<Target> target;

    int nFoils = target->nFoils();

    // Pick a number of muons from a Poisson distribution, or a fixed number
    // This number is not meant to be real - its just to exercise the system.
    long n;
    if (_mean<=0) n=(long)-_mean;
    else         n = RandPoisson::shoot(_mean);

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

      // Pick a random energy.  
      const double e = _elow + _funcGen->shoot() * ( _ehi - _elow );
      _piCaptureEPhot->Fill(e);
      _piCaptureEPhotZ->Fill(e);

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

  const double PiCapture::EPhotFunc(const double x)
  {
      // Parameters from doc 665-v2
      static const double emax=138.2;
      static const double alpha = 2.691;
      static const double gamma = 1.051;
      static const double tau = 8.043;
      static const double c0 = 2.741;
      static const double c1 = -0.005;

      return pow(emax-x,alpha) * exp(-(emax-gamma*x)/tau) * (c0 + c1*x);


  }// PiCapture::EPhotFunc

}
