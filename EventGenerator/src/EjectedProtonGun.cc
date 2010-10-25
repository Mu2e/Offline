//
//
// Simulate the protons that come from the stopping target when muons capture
// on an Al nucleus.  Use the MECO distribution for the kinetic energy of the
// protons.  Production is uniform across the targets and uniform in time;
// this model needs to be improved.
//
// $Id: EjectedProtonGun.cc,v 1.9 2010/10/25 21:12:44 onoratog Exp $ 
// $Author: onoratog $
// $Date: 2010/10/25 21:12:44 $
//
// Original author Rob Kutschke, heavily modified by R. Bernstein
// 
// 

// C++ incldues.
#include <iostream>

// Framework includes
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"

// Mu2e includes
#include "EventGenerator/inc/EjectedProtonGun.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "TargetGeom/inc/Target.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"

// General Utilities
#include "GeneralUtilities/inc/pow.hh"

// Other external includes.
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

//ROOT Includes
#include "TH1D.h"


using namespace std;

namespace mu2e {

  EjectedProtonGun::EjectedProtonGun( edm::Run& run, const SimpleConfig& config ):
    
    // Base class.
    GeneratorBase(),
    
    // Configurable parameters
    _mean(config.getDouble("ejectedProtonGun.mean",1.)),
    _elow(config.getDouble("ejectedProtonGun.elow",0.)),
    _ehi(config.getDouble("ejectedProtonGun.ehi",.100)),
    _czmin(config.getDouble("ejectedProtonGun.czmin",  0.3)),
    _czmax(config.getDouble("ejectedProtonGun.czmax",  0.6)),
    _phimin(config.getDouble("ejectedProtonGun.phimin", 0. )),
    _phimax(config.getDouble("ejectedProtonGun.phimax", CLHEP::twopi )),
    _nbins(config.getInt("ejectedProtonGun.nbins",1000)),
    _doHistograms(config.getBool("ejectedProtonGun.doHistograms",true)),

    // Initialize random number distributions; getEngine comes from the base class.
    _randPoissonQ( getEngine(), std::abs(_mean) ),

    // Histogram pointers
    _hMultiplicity(),
    _hKE(),
    _hKEZoom(),
    _hMomentumMeV(),
    _hzPosition(),
    _hcz(),
    _hphi(),
    _htime(){

    // About the ConditionsService:
    // The argument to the constructor is ignored for now.  It will be a
    // data base key.  There is a second argument that I have let take its
    // default value of "current"; it will be used to specify a version number.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    ConditionsHandle<DAQParams>         daqPar("ignored");
    
    // Default values for the start and end of the live window.
    // Can be overriden by the run-time config; see below.
    _tmin = daqPar->t0;
    _tmax = accPar->deBuncherPeriod;
    
    _tmin = config.getDouble("ejectedProtonGun.tmin",  _tmin );
    _tmax = config.getDouble("ejectedProtonGun.tmax",  _tmax );

    // Book histograms.
    if ( _doHistograms ){
      edm::Service<edm::TFileService> tfs;
      edm::TFileDirectory tfdir  = tfs->mkdir( "EjectedProtonGun" );
      _hMultiplicity = tfdir.make<TH1D>( "hMultiplicity", "Proton Multiplicity",                20,     0,     20  );
      _hKE           = tfdir.make<TH1D>( "hKE",           "Proton Kinetic Energy",              50, _elow,   _ehi  );
      _hMomentumMeV  = tfdir.make<TH1D>( "hMomentumMeV",  "Proton Momentum in MeV",             50, _elow,   _ehi  );
      _hKEZoom       = tfdir.make<TH1D>( "hEZoom",        "Proton Kinetic Energy (zoom)",      200, _elow,   _ehi  );
      _hzPosition    = tfdir.make<TH1D>( "hzPosition",    "Proton z Position (Tracker Coord)", 200, -6600., -5600. );
      _hcz           = tfdir.make<TH1D>( "hcz",           "Proton cos(theta)",                 100,    -1.,     1. );
      _hphi          = tfdir.make<TH1D>( "hphi",          "Proton azimuth",                    100,  -M_PI,  M_PI  );
      _htime         = tfdir.make<TH1D>( "htime",         "Proton time ",                      200,      0,  2000. );
    }

    FoilParticleGenerator generator(PDGCode::p_plus, GenId::ejectedProtonGun);
    fGenerator = generator;
    fGenerator._elow = _elow;
    fGenerator._ehi = _ehi;
    fGenerator._nbins = _nbins;
    fGenerator._czmin = _czmin;
    fGenerator._czmax = _czmax;
    fGenerator._phimin = _phimin;
    fGenerator._phimax = _phimax;
    fGenerator._tmin = _tmin;
    fGenerator._tmax = _tmax;

    fGenerator.setRandomEngine (GeneratorBase::getEngine() ); 

  }
  
  EjectedProtonGun::~EjectedProtonGun(){
  }
  
  void EjectedProtonGun::generate( ToyGenParticleCollection& genParts ){

    // Choose the number of protons to generate this event.
    long n = _mean < 0 ? static_cast<long>(-_mean): _randPoissonQ.fire();
    if ( _doHistograms ) { 
      _hMultiplicity->Fill(n); 
    }

    size_t genPartsEntries = genParts.size();    
    
    fGenerator.generateFromFoil( genParts, n );
    
    
    // Fill histograms.
    if ( _doHistograms) {
      for( size_t i=genPartsEntries; i<n+genPartsEntries ; i++) {     
        
        ToyGenParticle& particle = genParts[i];  
        
        _hKE->Fill( particle._momentum.e() );
        _hKEZoom->Fill( particle._momentum.e() );
        _hMomentumMeV->Fill( particle._momentum.mag() );
        _hzPosition->Fill( particle._position.z() );
        _hcz->Fill( particle._momentum.vect().cosTheta() );
        _hphi->Fill( particle._momentum.vect().phi() );
        _htime->Fill(particle._time);
      }
    }
  } // end generate
  
} // namespace mu2e
