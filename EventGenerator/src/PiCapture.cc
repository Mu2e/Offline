//
// Generate photons from pi- capture on Al nuclei.
// Based on Ivano Sarra's model described in mu2e Doc 665-v2
// add internal conversion, 11/2011 rhb
//
// $Id: PiCapture.cc,v 1.41 2014/01/27 22:20:17 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/01/27 22:20:17 $
//
// Original author Rob Kutschke/P. Shanahan
//

#include <iostream>

// Framework includes
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "EventGenerator/inc/PiCapture.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "StoppingTargetGeom/inc/zBinningForFoils.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// CLHEP includes
#include "CLHEP/Random/RandPoisson.h"



// ROOT includes
#include "TH1D.h"

using namespace std;

namespace mu2e {

  // Also need a home for this - the cycle time of the debuncher.
  static const double tcycle = 1694.;

  static const double emax = 138.2; //

  static const double internalRatio = 0.0069;

  PiCapture::PiCapture( art::Run& run, const SimpleConfig& config ):

    // Base class
    GeneratorBase(),

    // Parameters from run time configuration
    _mean(config.getDouble("picapture.mean", -1.)),
    _elow(config.getDouble("picapture.elow", 38.2)),
    _ehi(config.getDouble("picapture.ehi",   emax)),
    _probInternalConversion(config.getDouble("internalpicapture.mean", 0.0069)), //.0069 from Kroll and Wada
    _PStoDSDelay(config.getBool("picapture.PStoDSDelay", false)),
    _pPulseDelay(config.getBool("picapture.pPulseDelay", true)),
    _pPulseShift(config.getDouble("picapture.pPulseShift", 0)),
    _timeFolding(config.getBool("FoilParticleGenerator.foldingTimeOption", true)),
    _foilGen(config.getString("picapture.foilGen", "expoVolWeightFoil")),
    _posGen(config.getString("picapture.posGen", "flatPos")),
    _timeGen(config.getString("picapture.timeGen", "limitedExpoTime")),
    _nbins(config.getInt("picapture.nbins",  1000)),
    _doHistograms(config.getBool("picapture.doHistograms",true)),

    // Random number distributions; getEngine is found in the base class.
    _randPoissonQ( getEngine(), std::abs(_mean) ),

    _STfname(config.getString("FoilParticleGenerator.STfilename","ExampleDataFiles/StoppedMuons/stoppedMuons_02.txt")),
    _detSys(),
    // Histograms
    _hMultiplicity(0),
    _hEPhot(0),
    _hEPhotZ(0),
    _hEElect(0),
    _hEElectZ(0),
    _hzPos(0),
    _hcz(0),
    _hphi(0),
    _ht(0),
    _hmudelay(),
    _hpulsedelay(),
    _hFoilNumber(0){

    // Sanity checks
    if ( std::abs(_mean) > 99999. ) {
      throw cet::exception("RANGE")
        << "PiCapture has been asked to produce a crazily large number of photons ."
        << _mean
        << "\n";
    }
    if (_nbins <= 0) {
      throw cet::exception("RANGE")
        << "Nonsense picapture.nbins requested="
        << _nbins
        << "\n";
    }

    _detSys = &*GeomHandle<DetectorSystem>();

    // Book histograms.
   
    if ( _doHistograms ) {
      // Compute a binning that ensures that the stopping target foils are at bin centers.
      GeomHandle<StoppingTarget> target;
      Binning bins = zBinningForFoils(*target,7);

      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "PiCapture" );

      _hMultiplicity = tfdir.make<TH1D>( "hMultiplicity", "DIO Multiplicity",    20,     0.,      20. );

      _hEPhot      = tfdir.make<TH1D>( "hEPhot",  "PiCapture E(Photon)",              200,     0.,   200. );
      _hEElect      = tfdir.make<TH1D>( "hEElect",  "PiCapture Internal E(Electron)",  200,     0.,   200. );
      _hEPhotZ     = tfdir.make<TH1D>( "hEPhotZ", "PiCapture E(Photon)(zoom)",        200,  _elow,   _ehi );
      _hEElectZ     = tfdir.make<TH1D>( "hEElectZ", "PiCapture E(Electron)(zoom)",    200,  _elow,   _ehi );
      _hInternalFraction     = tfdir.make<TH1D>( "hInternalFraction", "Internal Conversion Fractional Distribution",50,0.,1.);
      _hzPos       = tfdir.make<TH1D>( "hzPos",   "PiCap z Position (Tracker Coord)", bins.nbins(), bins.low(), bins.high() );
      _hcz         = tfdir.make<TH1D>( "hcz",     "PiCap cos(theta)",                 100,    -1.,     1. );
      _hphi        = tfdir.make<TH1D>( "hphi",    "PiCapture azimuth",                100,  -M_PI,   M_PI );
      _ht          = tfdir.make<TH1D>( "ht",     "PiCapture time ",                  210,   -200.,  2000. );
      _hmudelay    = tfdir.make<TH1D>( "hmudelay",      "Production delay due to muons arriving at ST;(ns)", 600, 0., 3000. );
      _hpulsedelay = tfdir.make<TH1D>( "hpdelay",       "Production delay due to the proton pulse;(ns)", 60, -150., 150. );  
      _hFoilNumber = tfdir.make<TH1D>( "hFoilNumber", "Foil Number", 20,0.,20.);
      
    }
    // UNUSED CODE.
    // set up exponential RNG for foils.  this code is a complete hack which is why I'm not dignifying it
    // by putting something in the config file.
    
    //    foilMean = 1./0.693;
    
    //
    // 24.5 is from Rick Coleman telling me the stopped pi's fall a factor of two over the 17 foils


    _fGenerator = unique_ptr<FoilParticleGenerator>
      (new FoilParticleGenerator( getEngine(), 0, tcycle,
                                  FoilParticleGenerator::findFoilGenByName(_foilGen),
                                  FoilParticleGenerator::findPosGenByName(_posGen),
                                  FoilParticleGenerator::findTimeGenByName(_timeGen),
                                  _PStoDSDelay,
                                  _pPulseDelay,
                                  _pPulseShift,
                                  _STfname));

    _piCaptureInfo = new PiCaptureEffects(_probInternalConversion, _elow, _ehi, _nbins);

  } // end PiCapture::PiCapture
  
  
  PiCapture::~PiCapture(){
  }
  


  void PiCapture::generate( GenParticleCollection& genParticles ){

    // Choose the number of photons to generate this event.
    long n = _mean < 0 ? static_cast<long>(-_mean): _randPoissonQ.fire();
    if ( _doHistograms ) {
      _hMultiplicity->Fill(n);
    }
    for ( long i=0; i<n; ++i ){

      _piCaptureInfo->defineOutput();

      CLHEP::Hep3Vector pos(0,0,0);
      double time;
      _fGenerator->generatePositionAndTime(pos, time, _timeFolding);

      _hFoilNumber->Fill(static_cast<double>(_fGenerator->iFoil()));

      const CLHEP::Hep3Vector detPos(_detSys->toDetector(pos));

      if (_piCaptureInfo->doPhoton()) {
        genParticles.push_back( _piCaptureInfo->outputGamma(pos, time));
        if (_doHistograms) {
          _hEPhot->Fill(genParticles.back().momentum().e());
          _hEPhotZ->Fill(genParticles.back().momentum().e());
          _hcz->Fill(genParticles.back().momentum().vect().cosTheta());
          _hphi->Fill(genParticles.back().momentum().vect().phi());
        }
      }
      if (_piCaptureInfo->doElectron()){
        genParticles.push_back( _piCaptureInfo->outputElec(pos, time));
        if (_doHistograms) {
          _hEElect->Fill(genParticles.back().momentum().e());
          _hEElectZ->Fill(genParticles.back().momentum().e());
          _hcz->Fill(genParticles.back().momentum().vect().cosTheta());
          _hphi->Fill(genParticles.back().momentum().vect().phi());
        }
      }
      if (_piCaptureInfo->doPositron()){
        genParticles.push_back( _piCaptureInfo->outputPosit(pos, time));
        if (_doHistograms) {
          _hEElect->Fill(genParticles.back().momentum().e());
          _hEElectZ->Fill(genParticles.back().momentum().e());
          _hcz->Fill(genParticles.back().momentum().vect().cosTheta());
        _hphi->Fill(genParticles.back().momentum().vect().phi());
        }
      }
      
      if (_doHistograms) {
        
        _hInternalFraction->Fill(_piCaptureInfo->internalFraction());
        _hzPos->Fill(detPos.z());
        _ht->Fill(time);
        _hmudelay   ->Fill(_fGenerator->muDelay());
        _hpulsedelay->Fill(_fGenerator->pulseDelay());
      }
    }
  } // end PiCapture::generate

} // namespace mu2e
