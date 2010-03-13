//
// Muon generator, uses Daya Bay libraries
//
// $Id: CosmicDYB.cc,v 1.3 2010/03/13 00:29:01 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/03/13 00:29:01 $
//
// Original author Yury Kolomensky
//

// C++ includes.
#include <iostream>
#include <math.h>


// Framework includes.
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"

// Mu2e includes.
#include "EventGenerator/inc/CosmicDYB.hh"
#include "EventGenerator/src/rm48.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TargetGeom/inc/Target.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "RandomNumberService/inc/RandomNumberService.hh"


// From CLHEP
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

// From Root.
#include "TH1D.h"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

// declare DYB functions here
extern "C" {
  void hrndg2_(double[100][4000],const long*,double*,double*,const long*,double*,double*,double*,double*,double*,float*);
}

namespace mu2e {

  static const long _ne = 4000;
  static const long _nth = 100;

  // Mass of the muon, in GeV.
  // Once we have the HepPDT package installed, get the mass from there.
  static const double m = 105.6584; // MeV

  CosmicDYB::CosmicDYB( edm::Run& run, const SimpleConfig& config ):
    GeneratorBase(){

    _mean = config.getDouble("cosmicDYB.mean",0.);
    _muEMin = config.getDouble("cosmicDYB.muEMin",0.);
    _muEMax = config.getDouble("cosmicDYB.muEMax",1.e5);

    // magic number from DYB ? 
    _muCosThMin = config.getDouble("cosmicDYB.muCosThMin",0.00366518);
    _muCosThMax = config.getDouble("cosmicDYB.muCosThMax",1.0);

    // size of the area to illuminate (in mm)
    _dx = config.getDouble("cosmicDYB.dx",5000);
    _dz = config.getDouble("cosmicDYB.dz",5000);
    _y0 = config.getDouble("cosmicDYB.y0",0);


    edm::LogInfo log("COSMIC");
    log << "cosmicDYB.mean = " << _mean << "\n"
	<< "cosmicDYB.muEMin = " << _muEMin 
	<< ", cosmicDYB.muEMax = " << _muEMax << "\n"
	<< "cosmicDYB.muCosThMin = " << _muCosThMin 
	<< ", cosmicDYB.muCosThMax = " << _muCosThMax << "\n"
	<< "cosmicDYB.dx = " << _dx 
	<< ", cosmicDYB.dz = " << _dz
	<< ", cosmicDYB.y0 = " << _y0 << "\n";

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
    _cosmicMultiplicityH = tfs->make<TH1D>( "cosmicMultiplicityH", "Cosmic Multiplicity", 20, -0.5, 19.5);

    // log of muon energy (GeV)
    _cosmicMomH = tfs->make<TH1D>( "cosmicMomH", "log (Momentum, GeV)", 60, -3, 3);

    // charge
    _cosmicChargeH = tfs->make<TH1D>( "cosmicChargeH", "Muon Charge", 2, -2, 2.);

    // Initialize fake RM48 that is used by DYB code.
    edm::Service<edm::RandomNumberGenerator> rng;
    static CLHEP::RandFlat flat(rng->getEngine());
    static bool init(true);
    if ( init ){
      init = false;
      setRm48Distribution(flat);
    }

    // initialize DYB generator
    float par = 1.;
    // convert to GeV
    _muEMin *= 1e-3;
    _muEMax *= 1e-3;
    
    double dim_sum,E,cosTh;
    ::hrndg2_(_buffer,&_ne,&_muEMin,&_muEMax,&_nth,&_muCosThMin,&_muCosThMax,
	      &dim_sum,&E,&cosTh,&par);

    // rate is per cm^2. The constants are 2*pi times the area
    double tRate = dim_sum*M_PI*0.08*_dx*_dz;
    log << "Total cosmic rate = " << tRate << " Hz\n";
    
  }

  CosmicDYB::~CosmicDYB(){
  }

  void CosmicDYB::generate( ToyGenParticleCollection& genParts ){

    if ( _mean <= -99999. ) return;

    // Should get the numbers here from the config file or from the
    // geometry manager.

    // Pick a number of muons from a Poisson distribution.
    long n;
    if (_mean<0) {
      n=(long)-_mean;
    } else {
      n = RandPoisson::shoot(_mean);
    }

    _cosmicMultiplicityH->Fill(n);

    for ( int i=0; i<n; ++i ){

      float par = 111.;  // double precision
      double dim_sum,E,cosTh;
      ::hrndg2_(_buffer,&_ne,&_muEMin,&_muEMax,&_nth,&_muCosThMin,&_muCosThMax,
		&dim_sum,&E,&cosTh,&par);
      
      // energy is in GeV, convert to MeV
      E *= 1e3;

      double p = 0;
      if ( E<= m ) {
	E = m;
      } else {
	p = sqrt(E*E-m*m);
      }

      _cosmicMomH->Fill(log10(p)-3);

      // Cosine and sin of polar angle wrt y axis.
      double cy = cosTh;
      double sy = sin(acos(cosTh));

      double phi = 2.*M_PI*RandFlat::shoot();
     
      HepLorentzVector mom(p*sy*cos(phi), -p*cy, p*sy*sin(phi), E);

      // Position in a reference plane that is just above the ground.
      // We can worry later about the exact meaning of the height.
      // The G4 interface code ( PrimaryGeneratorAction) will put it
      // at the right height.
      // units of (x,y,z)
      double x = (1.-2.*RandFlat::shoot())*_dx;
      double y = _y0;
      double z = (1.-2.*RandFlat::shoot())*_dz;
      Hep3Vector pos( x, y, z );

      double time = _tmin + _dt*RandFlat::shoot();

      // Pick a random charge.
      // implement a rough charge asymmetry
      double logP = log10(p)-3.;
      double asym = 1.15;
      if ( logP > 0 ) {
	if ( logP > 1 ) {
	  asym = 1.3;
	} else {
	  asym = 1.15+0.15*logP;
	}
      }


      PDGCode::type pid = (RandFlat::shoot() > asym/(1+asym) ) ? PDGCode::mu_minus : PDGCode::mu_plus;

      _cosmicChargeH->Fill(-pid/abs(pid));

      // Add the cosmic to  the list.
      genParts.push_back( ToyGenParticle( pid, GenId::cosmicDYB, pos, mom, time));
    }
  }

}
