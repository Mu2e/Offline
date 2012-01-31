//
// Generate photons from pi- capture on Al nuclei.
// Based on Ivano Sarra's model described in mu2e Doc 665-v2
// add internal conversion, 11/2011 rhb
//
// $Id: PiCapture.cc,v 1.32 2012/01/31 05:34:19 onoratog Exp $
// $Author: onoratog $
// $Date: 2012/01/31 05:34:19 $
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
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "EventGenerator/inc/PiCapture.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"

// CLHEP includes
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandFlat.h"


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
    _nbins(config.getInt("picapture.nbins",  1000)),
    _doHistograms(config.getBool("picapture.doHistograms",true)),

    // Random number distributions; getEngine is found in the base class.
    _randPoissonQ( getEngine(), std::abs(_mean) ),
    _randomUnitSphere( getEngine() ),
    _randFlat( getEngine() ),
    _spectrum( getEngine(), &(binnedEnergySpectrum()[0]),_nbins),
    _internalFractionalSpectrum( getEngine(), &(internalFractionalBinnedSpectrum()[0]),_nbins),

    _STfname(config.getString("FoilParticleGenerator.STfilename","ExampleDataFiles/StoppedMuons/stoppedMuons_02.txt")),
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

    ConditionsHandle<ParticleDataTable> pdt("ignored");
    // pick up particle mass
    const HepPDT::ParticleData& e_data = pdt->particle(PDGCode::e_minus).ref();
    _electMass = e_data.mass().value();

    // Book histograms.
    if ( _doHistograms ){
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "PiCapture" );

      _hMultiplicity = tfdir.make<TH1D>( "hMultiplicity", "DIO Multiplicity",    20,     0.,      20. );

      _hEPhot      = tfdir.make<TH1D>( "hEPhot",  "PiCapture E(Photon)",              200,     0.,   200. );
      _hEElect      = tfdir.make<TH1D>( "hEElect",  "PiCapture Internal E(Electron)",  200,     0.,   200. );
      _hEPhotZ     = tfdir.make<TH1D>( "hEPhotZ", "PiCapture E(Photon)(zoom)",        200,  _elow,   _ehi );
      _hEElectZ     = tfdir.make<TH1D>( "hEElectZ", "PiCapture E(Electron)(zoom)",    200,  _elow,   _ehi );
      _hInternalFraction     = tfdir.make<TH1D>( "hInternalFraction", "Internal Conversion Fractional Distribution",50,0.,1.);
      _hzPos       = tfdir.make<TH1D>( "hzPos",   "PiCap z Position (Tracker Coord)", 200, -6600., -5600. );
      _hcz         = tfdir.make<TH1D>( "hcz",     "PiCap cos(theta)",                 100,    -1.,     1. );
      _hphi        = tfdir.make<TH1D>( "hphi",    "PiCapture azimuth",                100,  -M_PI,   M_PI );
      _ht          = tfdir.make<TH1D>( "ht",     "PiCapture time ",                  210,   -200.,  2000. );
      _hmudelay    = tfdir.make<TH1D>( "hmudelay",      "Production delay due to muons arriving at ST;(ns)", 600, 0., 3000. );
      _hpulsedelay = tfdir.make<TH1D>( "hpdelay",       "Production delay due to the proton pulse;(ns)", 60, 0., 300. );  
      _hFoilNumber = tfdir.make<TH1D>( "hFoilNumber", "Foil Number", 20,0.,20.);

    }
    // UNUSED CODE.
    // set up exponential RNG for foils.  this code is a complete hack which is why I'm not dignifying it
    // by putting something in the config file.

    //    foilMean = 1./0.693;

    //
    // 24.5 is from Rick Coleman telling me the stopped pi's fall a factor of two over the 17 foils

    _fGenerator = auto_ptr<FoilParticleGenerator>( new FoilParticleGenerator(getEngine(), 0 ,tcycle,
                                                                             FoilParticleGenerator::expoVolWeightFoil,
                                                                             FoilParticleGenerator::flatPos,
                                                                             FoilParticleGenerator::limitedExpoTime,
									     _PStoDSDelay,
                                                                             _pPulseDelay,
									     _pPulseShift,
									     _STfname));

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
      //Pick up position and momentum
      CLHEP::Hep3Vector pos(0,0,0);
      double time;
      _fGenerator->generatePositionAndTime(pos, time);

      _hFoilNumber->Fill(static_cast<double>(_fGenerator->iFoil()));

      // Pick a random photon energy from the spectrum.
      const double e = _elow + _spectrum.fire() * ( _ehi - _elow );

      // Make the 4 vector.
      CLHEP::HepLorentzVector mom( _randomUnitSphere.fire(e), e);

      bool internal = true;
      //
      // pull out a random number to determine whether regular RPC or internal conversion process
      if (_randFlat.fire() > _probInternalConversion) internal = false;
      if (!internal) {
        genParticles.push_back( GenParticle( PDGCode::gamma, GenId::pionCapture, pos, mom, time));
      }
      if (internal){
        // Pick a random photon energy from the spectrum.
        double electMom;
        double positMom;
        const double fract = _internalFractionalSpectrum.fire();
        _eElect = e*fract;
        _ePosit = e*(1. - fract);
        //
        // make the momentum four-vector.  the formula in Rossi assumed e>> m_e.  RPCs of interest ~100 MeV, so 
        // just fudge this so there's never a negative sqrt.  Assumes CPT.  
        electMom= safeSqrt((_eElect)*(_eElect) - _electMass*_electMass);
        positMom = safeSqrt((_ePosit)*(_ePosit) - _electMass*_electMass);
        CLHEP::HepLorentzVector momElectVec( _randomUnitSphere.fire(electMom), electMom);
        CLHEP::HepLorentzVector momPositVec( _randomUnitSphere.fire(positMom), positMom);
        //
        // only stack particles with momenta > 0 or someone segfaults
        if (electMom > 0.) {genParticles.push_back( GenParticle( PDGCode::e_minus, GenId::internalRPC, pos, momElectVec, time));}
        if (positMom > 0.) {genParticles.push_back( GenParticle( PDGCode::e_plus,  GenId::internalRPC, pos, momPositVec, time));}
        if (_doHistograms){
          _hEElect->Fill(_eElect);
          _hEElect->Fill(_ePosit);
          _hEElectZ->Fill(_eElect);
          _hEElectZ->Fill(_ePosit);
          _hInternalFraction->Fill(fract);
        }
    
      }
      if ( _doHistograms ){
        _hEPhot->Fill(e);
        _hEPhotZ->Fill(e);    
        _hzPos->Fill(pos.z());
        _hcz->Fill(mom.vect().cosTheta());
        _hphi->Fill(mom.vect().phi());
        _ht->Fill(time);
        _hmudelay   ->Fill(_fGenerator->muDelay());
        _hpulsedelay->Fill(_fGenerator->pulseDelay());
      }
    }
  } // end loop over photons to generate


  // Photon energy spectrum as a continuous function.
  double PiCapture::energySpectrum(const double x)
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
      throw cet::exception("RANGE")
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



  // Photon energy spectrum as a continuous function.
  double PiCapture::internalFractionalSpectrum(const double x)
  {
    // the usual distribution for positron/electron fraction, normalized to unity
    return (9./7.)*(1. - (4/3)*x*(1 - x));
  } // PiCapture::internalFractionalSpectrum

  // Compute a binned representation of the photon energy spectrum.
  // This runs from 0 to 1; we'll take the photon energy and multiply.
  std::vector<double> PiCapture::internalFractionalBinnedSpectrum(){

    // Sanity check.
    if (_nbins <= 0) {
      throw cet::exception("RANGE")
        << "Nonsense InternalPiCaptureGun.nbins requested="
        << _nbins
        << "\n";
    }

    // Bin width.
    double dy = 1./_nbins;

    // Vector to hold the binned representation of the energy spectrum.
    std::vector<double> spectrum;
    spectrum.reserve(_nbins);

    for (int ib=0; ib<_nbins; ib++) {
      double x = (ib+0.5) * dy;
      spectrum.push_back(internalFractionalSpectrum(x));
    }

    return spectrum;
  } // PiCapture::internalFractionalBinnedSpectrum




} // namespace mu2e
