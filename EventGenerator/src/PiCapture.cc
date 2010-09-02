//
// Generate photons from pi- capture on Al nuclei.
// Based on Ivano Sarra's model described in mu2e Doc 665-v2
//
// $Id: PiCapture.cc,v 1.11 2010/09/02 18:26:02 rhbob Exp $
// $Author: rhbob $ 
// $Date: 2010/09/02 18:26:02 $
//
// Original author Rob Kutschke/P. Shanahan
// 

#include <iostream>

// Framework includes
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"

// Mu2e includes
#include "EventGenerator/inc/PiCapture.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
//#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TargetGeom/inc/Target.hh"

// CLHEP includes
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

// ROOT includes
#include "TH1D.h"

using namespace std;

namespace mu2e {

  // Also need a home for this - the cycle time of the debuncher.
  static const double tcycle = 1694.;

  static const double emax = 138.2; // 

  PiCapture::PiCapture( edm::Run& run, const SimpleConfig& config ):
    
    // Base class
    GeneratorBase(),

    // Parameters from run time configuration
    _mean(config.getDouble("picapture.mean", -1.)),
    _elow(config.getDouble("picapture.elow", 38.2)),
    _ehi(config.getDouble("picapture.ehi",   emax)),
    _nbins(config.getInt("picapture.nbins",  1000)),
    _doHistograms(config.getBool("picapture.doHistograms",true)),

    // Histograms
    _hEPhot(0),
    _hEPhotZ(0),

    // Random number distributions; getEngine is found in the base class.
    _randFlat( getEngine() ),
    _randPoissonQ( getEngine(), std::abs(_mean) ),
    _randExponential( getEngine(), 1.),
    _randomUnitSphere( getEngine() ),
    _spectrum( getEngine(), &(binnedEnergySpectrum()[0]),_nbins){

    // Sanity checks
    if ( std::abs(_mean) > 99999. ) {
      throw cms::Exception("RANGE") 
        << "PiCapture has been asked to produce a crazily large number of photons ."
        << _mean
        << "\n";
    }
    if (_nbins <= 0) {
      throw cms::Exception("RANGE") 
        << "Nonsense picapture.nbins requested="
        << _nbins
        << "\n";
    }

    // Book histograms.
    if ( _doHistograms ){
      edm::Service<edm::TFileService> tfs;
      edm::TFileDirectory tfdir = tfs->mkdir( "PiCapture" );

      _hMultiplicity = tfdir.make<TH1D>( "hMultiplicity", "DIO Multiplicity",    20,     0.,      20. );

      _hEPhot      = tfdir.make<TH1D>( "hEPhot",  "PiCapture E(Photon)",              200,     0.,   200. );
      _hEPhotZ     = tfdir.make<TH1D>( "hEPhotZ", "PiCapture E(Photon)(zoom)",        200,  _elow,   _ehi );
      _hzPos       = tfdir.make<TH1D>( "hzPos",   "PiCap z Position (Tracker Coord)", 200, -6600., -5600. );
      _hcz         = tfdir.make<TH1D>( "hcz",     "PiCap cos(theta)",                 100,    -1.,     1. );
      _hphi        = tfdir.make<TH1D>( "hphi",    "PiCapture azimuth",                100,  -M_PI,   M_PI );
      _ht          = tfdir.make<TH1D>( "ht",     "PiCapture time ",                  100,      0,  2000. );
      _hFoilNumber = tfdir.make<TH1D>( "hFoilNumber", "Foil Number", 20,0.,20.);

    }
     //
      // set up exponential RNG for foils.  this code is a complete hack which is why I'm not dignifying it
    // by putting something in the config file.

    foilMean = 1./0.693;

    //
    // 24.5 is from Rick Coleman telling me the stopped pi's fall a factor of two over the 17 foils  

  } // end PiCapture::PiCapture

  PiCapture::~PiCapture(){
  }

  void PiCapture::generate( ToyGenParticleCollection& genParticles ){

    // Choose the number of photons to generate this event.
    long n = _mean < 0 ? static_cast<long>(-_mean): _randPoissonQ.fire();
    if ( _doHistograms ) { 
      _hMultiplicity->Fill(n); 
    }

    // Get access to the geometry system.
    GeomHandle<Target> target;
    nFoils = target->nFoils();

    for ( long i=0; i<n; ++i ){

      // Pick a foil.  Exponentials can go past the number of foils, so have to do a trick (and recall foils numbered
      // 0--16 so do while <17
      //int ifoil = static_cast<int>(nFoils*_randFlat.fire());
      uint32_t ifoil;
      do{
        ifoil = static_cast<int>(nFoils*_randExponential.fire());
      }while(ifoil > nFoils-1);
      //cout << "ifoil = " << ifoil << endl;
      _hFoilNumber->Fill(static_cast<double>(ifoil));
      TargetFoil const& foil = target->foil(ifoil);

      // Foil properties.
      CLHEP::Hep3Vector const& center = foil.center();
      const double r1 = foil.rIn();
      const double dr = foil.rOut() - r1;

      // A random point within the foil.
      const double r   = r1 + dr*_randFlat.fire();
      const double dz  = 2.*(0.5-_randFlat.fire())*foil.halfThickness();
      const double phi = 2.*M_PI*_randFlat.fire();
      CLHEP::Hep3Vector pos( center.x()+r*cos(phi), 
                             center.y()+r*sin(phi), 
                             center.z()+dz );

      // Pick a random photon energy from the spectrum.
      const double e = _elow + _spectrum.fire() * ( _ehi - _elow );

      // This should reflect the model of the tails of the beam
      // and the out-of-time protons.
      const double time = 1694*_randFlat.fire();

      // Make the 4 vector.
      CLHEP::HepLorentzVector mom( _randomUnitSphere.fire(e), e);
      
      // Add the photon to the list of generated particles.
      genParticles.push_back( ToyGenParticle( PDGCode::gamma, GenId::pionCapture, pos, mom, time));

      if ( _doHistograms ){
        _hEPhot->Fill(e);
        _hEPhotZ->Fill(e);
        _hzPos->Fill(pos.z());
        _hcz->Fill(mom.vect().cosTheta());
        _hphi->Fill(mom.vect().phi());
        _ht->Fill(time);
      }

    } // end loop over photons to generate

  }

  // Photon energy spectrum as a continuous function.
  const double PiCapture::energySpectrum(const double x)
  {
    // Parameters from doc 665-v2
    static const double emax  = 138.2;
    static const double alpha =   2.691;
    static const double gamma =   1.051;
    static const double tau   =   8.043;
    static const double c0    =   2.741;
    static const double c1    =  -0.005;

    return pow(emax-x,alpha) * exp(-(emax-gamma*x)/tau) * (c0 + c1*x);

  } // PiCapture::energySpectrum

  // Compute a binned representation of the photon energy spectrum.
  std::vector<double> PiCapture::binnedEnergySpectrum(){

    // Sanity check.
    if (_nbins <= 0) {
      throw cms::Exception("RANGE") 
        << "Nonsense PiCaptureGun.nbins requested="
        << _nbins
        << "\n";
    }

    // Bin width.
    double dE = (_ehi - _elow) / _nbins;

    // Vector to hold the binned representation of the energy spectrum.
    std::vector<double> spectrum;
    spectrum.reserve(_nbins);
    
    for (int ib=0; ib<_nbins; ib++) {
      double x = _elow+(ib+0.5) * dE;
      spectrum.push_back(energySpectrum(x));
    }

    return spectrum;
  } // PiCapture::binnedEnergySpectrum

} // namespace mu2e
