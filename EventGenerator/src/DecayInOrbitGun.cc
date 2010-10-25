//
// Generate some number of DIO electrons.
//
// $Id: DecayInOrbitGun.cc,v 1.11 2010/10/25 19:50:21 onoratog Exp $ 
// $Author: onoratog $
// $Date: 2010/10/25 19:50:21 $
//
// Original author Rob Kutschke
// 
//
// 

// C++ includes.
#include <iostream>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"

// Mu2e includes
#include "EventGenerator/inc/DecayInOrbitGun.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"


//ROOT includes
#include "TH1D.h"

// CLHEP includes.
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;

//using CLHEP::Hep3Vector;
//using CLHEP::HepLorentzVector;
//using CLHEP::RandFlat;
//using CLHEP::twopi;

namespace mu2e {

  // Need a Conditions entity to hold info about conversions:
  // endpoints and lifetimes for different materials etc
  // Grab them from Andrew's minimc package?
  static const double conversionEnergyAluminum = 104.96;

  DecayInOrbitGun::DecayInOrbitGun( edm::Run& run, const SimpleConfig& config ):

    // Base class
    GeneratorBase(),

    // Information from config file.
    _mean(config.getDouble("decayinorbitGun.mean",1.)),
    _elow(config.getDouble("decayinorbitGun.elow",100.)),
    _ehi(config.getDouble("decayinorbitGun.ehi",conversionEnergyAluminum)),
    _nbins(config.getInt("decayinorbitGun.nbins",1000)),
    _czmin(config.getDouble("decayinorbitGun.czmin",  0.3)),
    _czmax(config.getDouble("decayinorbitGun.czmax",  0.6)),
    _phimin(config.getDouble("decayinorbitGun.phimin", 0. )),
    _phimax(config.getDouble("decayinorbitGun.phimax", CLHEP::twopi )),
    _doHistograms(config.getBool("decayinorbitGun.doHistograms", true)),

    // Random number distributions; getEngine comes from the base class.
    _randPoissonQ( getEngine(), std::abs(_mean) ),
  
    // Histograms.
    _hMultiplicity(0),
    _hEElec(0),
    _hEElecZ(0),
    _hzPosition(0),
    _hcz(0),
    _hphi(0),
    _ht()  {

    // Sanity check.
    if ( std::abs(_mean) > 99999. ) {
      throw cms::Exception("RANGE") 
        << "DecayInOrbit Gun has been asked to produce a crazily large number of electrons."
        << _mean
        << "\n";
    }

    // About the ConditionsService:
    // The argument to the constructor is ignored for now.  It will be a
    // data base key.  There is a second argument that I have let take its
    // default value of "current"; it will be used to specify a version number.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    ConditionsHandle<DAQParams>         daqPar("ignored");
    
    _tmin   = config.getDouble("decayinorbitGun.tmin",  daqPar->t0 );
    _tmax   = config.getDouble("decayinorbitGun.tmax",  accPar->deBuncherPeriod );

    // Make ROOT subdirectory to hold diagnostic histograms; book those histograms.
    if ( _doHistograms ){
      edm::Service<edm::TFileService> tfs;
      edm::TFileDirectory tfdir = tfs->mkdir( "DecayInOrbit" );
      _hMultiplicity = tfdir.make<TH1D>( "hMultiplicity", "DIO Multiplicity",                20,     0.,   20.    );
      _hEElec        = tfdir.make<TH1D>( "hEElec",        "DIO Electron Energy",             10,  _elow,    0.105 );
      _hEElecZ       = tfdir.make<TH1D>( "hEElecZ",       "DIO Electron Energy (zoom)",     200,  _elow, _ehi     );
      _hzPosition    = tfdir.make<TH1D>( "hzPosition",    "DIO z Position (Tracker Coord)", 200, -6600., -5600.   );
      _hcz           = tfdir.make<TH1D>( "hcz",           "DIO cos(theta)",                 100,    -1.,     1.   );
      _hphi          = tfdir.make<TH1D>( "hphi",          "DIO azimuth",                    100,  -M_PI,   M_PI   );
      _ht            = tfdir.make<TH1D>( "ht",            "DIO time ", 200, 0, 2000. );
    }

    
    FoilParticleGenerator generator(PDGCode::e_minus, GenId::dio1);
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
    fGenerator._conversionEnergyAluminum = conversionEnergyAluminum;

    fGenerator.setRandomEngine (GeneratorBase::getEngine() ); 
    
  }

  DecayInOrbitGun::~DecayInOrbitGun(){
  }

  void DecayInOrbitGun::generate( ToyGenParticleCollection& genParts ){

    // Choose the number of electrons to generate this event.
    long n = _mean < 0 ? static_cast<long>(-_mean): _randPoissonQ.fire();
      if ( _doHistograms ) { 
        _hMultiplicity->Fill(n); 
      }

      size_t genPartsEntries = genParts.size();    

      fGenerator.generateFromFoil( genParts, n );
    
      if( _doHistograms ){
        for( size_t i=genPartsEntries; i<n+genPartsEntries ; i++) {     

          ToyGenParticle& particle = genParts[i];  

          _hEElec    ->Fill( particle._momentum.e() );
          _hEElecZ   ->Fill( particle._momentum.e() );
          _hzPosition->Fill( particle._position.z() );
          _hcz       ->Fill( particle._momentum.vect().cosTheta() );
          _hphi      ->Fill( particle._momentum.vect().phi() );
          _ht        ->Fill( particle._time );
        
        }
      }

  } // DecayInOrbitGun::generate

}
