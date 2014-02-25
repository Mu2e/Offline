//
// Generate some number of DIO electrons.
//
// $Id: DecayInOrbitGun.cc,v 1.62 2014/02/25 17:14:10 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/02/25 17:14:10 $
//
// Original author Rob Kutschke
//

// CLHEP includes
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Framework includes
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib/pow.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/PhysicsParams.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "EventGenerator/inc/DecayInOrbitGun.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "StoppingTargetGeom/inc/zBinningForFoils.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/CzarneckiSpectrum.hh"
#include "Mu2eUtilities/inc/ShankerWatanabeSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"

// ROOT includes
#include "TH1D.h"

// C++ includes
#include <iostream>

using namespace std;

namespace mu2e {

  DecayInOrbitGun::DecayInOrbitGun( art::Run& run, const SimpleConfig& config ):

    // Base class
    GeneratorBase(),

    // Information from config file.
    _mean(config.getDouble("decayinorbitGun.mean",1.)),
    _elow(config.getDouble("decayinorbitGun.elow",0.)),
    _ehi(config.getDouble("decayinorbitGun.ehi",GlobalConstantsHandle<PhysicsParams>()->getEndpointEnergy())),
    _nbins(config.getInt("decayinorbitGun.nbins",1000)),
    _czmin(config.getDouble("decayinorbitGun.czmin", -1.0)),
    _czmax(config.getDouble("decayinorbitGun.czmax",  1.0)),
    _phimin(config.getDouble("decayinorbitGun.phimin", 0. )),
    _phimax(config.getDouble("decayinorbitGun.phimax", CLHEP::twopi )),
    _pStodSDelay(config.getBool("decayinorbitGun.PStoDSDelay", true)),
    _pPulseDelay(config.getBool("decayinorbitGun.pPulseDelay", false)),
    _pPulseShift(config.getDouble("decayinorbitGun.pPulseShift", 0)),
    _timeFolding(config.getBool("FoilParticleGenerator.foldingTimeOption", true)),
    _foilGen(config.getString("decayinorbitGun.foilGen", "muonFileInputFoil")),
    _posGen(config.getString("decayinorbitGun.posGen", "muonFileInputPos")),
    _timeGen(config.getString("decayinorbitGun.timeGen", "negExp")),
    _doHistograms(config.getBool("decayinorbitGun.doHistograms", true)),
    _spectrumResolution(config.getDouble("decayinorbitGun.spectrumResolution", 0.1)),
    _energySpectrum(config.getString("decayinorbitGun.energySpectrum", "Czarnecki")),
    _stFname(config.getString("FoilParticleGenerator.STfilename","ExampleDataFiles/StoppedMuons/stoppedMuons_02.txt")),
    _nToSkip (config.getInt("decayinorbitGun.nToSkip",0)),

    _detSys(),

    // Histograms.
    _hMultiplicity(0),
    _hEElec(0),
    _hEElecZ(0),
    _hradius(0),
    _hzPosition(0),
    _hcz(0),
    _hphi(0),
    _ht(),
    _hmudelay(),
    _hpulsedelay()  {

    //pick up particle mass
    GlobalConstantsHandle<ParticleDataTable> pdt;
    const HepPDT::ParticleData& e_data = pdt->particle(PDGCode::e_minus).ref();
    _mass = e_data.mass().value();

    // Sanity checks.
    if ( std::abs(_mean) > 99999. ) {
      throw cet::exception("RANGE")
        << "DecayInOrbit Gun has been asked to produce a crazily large number of electrons."
        << _mean
        << "\n";
    }

    GlobalConstantsHandle<PhysicsParams> phy;
    if (_ehi > (phy->getEndpointEnergy()+1)) { // 1 MeV of tolerance for energy range
      throw cet::exception("RANGE")
        << "Generation energy range must be within 0 and " << phy->getEndpointEnergy() << " (plus 1 MeV of tolerance)"
        << '\n';
    }

    if (_elow <= _mass) {
      _elow = _mass + 0.001;
      cout << "Lower bound of the DIO spectrum must be higher than electron mass. Set to electron mass + 0.001 MeV" << endl;
    }
    if ((_energySpectrum == "pol5" || _energySpectrum == "pol58" || _energySpectrum == "flat") && _ehi <=_mass) {
      throw cet::exception("RANGE")
        << "Using a not tabulated spectrum provides that the minimum energy must be higher than "
        << "electron mass"
        << '\n';
    }

    // About the ConditionsService:
    // The argument to the constructor is ignored for now.  It will be a
    // data base key.  There is a second argument that I have let take its
    // default value of "current"; it will be used to specify a version number.
    ConditionsHandle<AcceleratorParams> accPar("ignored");

    _tmin   = config.getDouble("decayinorbitGun.tmin",  0. );
    _tmax   = config.getDouble("decayinorbitGun.tmax",  accPar->deBuncherPeriod );

    _detSys = &*GeomHandle<DetectorSystem>();

    // Make ROOT subdirectory to hold diagnostic histograms; book those histograms.

    if ( _doHistograms ){
      // Compute a binning that ensures that the stopping target foils are at bin centers.
      GeomHandle<StoppingTarget> target;
      Binning bins = zBinningForFoils(*target,7);
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "DecayInOrbit" );
      _hMultiplicity = tfdir.make<TH1D>( "hMultiplicity", "DIO Multiplicity",               100,     0.,   200.   );
      _hEElec        = tfdir.make<TH1D>( "hEElec",        "DIO Electron Energy",            100,     0.,   105.   );
      _hEElecZ       = tfdir.make<TH1D>( "hEElecZ",       "DIO Electron Energy (zoom)",     200, _elow,   _ehi     );
      _hradius       = tfdir.make<TH1D>( "hradius",       "DIO radius(Tracker Coord)",      150,     0.,   150.   );
      _hzPosition    = tfdir.make<TH1D>( "hzPosition",    "DIO z Position (Tracker Coord)", bins.nbins(), bins.low(), bins.high() );
      _hcz           = tfdir.make<TH1D>( "hcz",           "DIO cos(theta)",                 100,    -1.,     1.   );
      _hphi          = tfdir.make<TH1D>( "hphi",          "DIO azimuth",                    100,  -M_PI,   M_PI   );
      _ht            = tfdir.make<TH1D>( "ht",            "DIO time ", 210, -200., 3000. );
      _hmudelay      = tfdir.make<TH1D>( "hmudelay",      "Production delay due to muons arriving at ST;(ns)", 300, 0., 2000. );
      _hpulsedelay   = tfdir.make<TH1D>( "hpdelay",       "Production delay due to the proton pulse;(ns)", 60, -150., 150. );
    }

    _fGenerator = unique_ptr<FoilParticleGenerator>
      (new FoilParticleGenerator( getEngine(), _tmin, _tmax,
                                  FoilParticleGenerator::findFoilGenByName(_foilGen),
                                  FoilParticleGenerator::findPosGenByName(_posGen),
                                  FoilParticleGenerator::findTimeGenByName(_timeGen),
                                  _pStodSDelay,
                                  _pPulseDelay,
                                  _pPulseShift,
                                  _stFname,
                                  _nToSkip));

    // Pick correct model
    if      ( _energySpectrum == "ShankerWatanabe" ) { 
      _dioGenId = GenId::dioShankerWatanabe;       
      _dioSpectrum.initialize< ShankerWatanabeSpectrum >( _elow, _ehi, _spectrumResolution );
    }
    else if ( _energySpectrum == "Czarnecki"       ) { 
      _dioGenId = GenId::dioCzarnecki; 
      _dioSpectrum.initialize< CzarneckiSpectrum >( _elow, _ehi, _spectrumResolution );
    }
    else if ( _energySpectrum == "flat"            ) { 
      _dioGenId = GenId::dioFlat;          
      _dioSpectrum.initialize< SimpleSpectrum >( _elow, _ehi, _spectrumResolution, SimpleSpectrum::Spectrum::FlatTrunc );
    }
    else if ( _energySpectrum == "pol5"            ) { 
      _dioGenId = GenId::dioE5; 
      _dioSpectrum.initialize< SimpleSpectrum >( _elow, _ehi, _spectrumResolution, SimpleSpectrum::Spectrum::Pol5 );
    }
    else if ( _energySpectrum == "pol58"           ) { 
      _dioGenId = GenId::dioE58;             
      _dioSpectrum.initialize< SimpleSpectrum >( _elow, _ehi, _spectrumResolution, SimpleSpectrum::Spectrum::Pol58 );
    }
    else {
      throw cet::exception("MODEL")
        << "Wrong or not allowed DIO energy spectrum";
    }
    
  }

  DecayInOrbitGun::~DecayInOrbitGun(){
  }

  void DecayInOrbitGun::generate( GenParticleCollection& genParts ){

    // Choose the number of electrons to generate this event.
    static CLHEP::RandPoissonQ randPoissonQ( getEngine(), std::abs(_mean) );

    long n = (_mean < 0 ? static_cast<long>(-_mean): randPoissonQ.fire());

    if ( _doHistograms ) {
      _hMultiplicity->Fill(n);
    }

    //Loop over particles to generate

    for (int i=0; i<n; ++i) {

      //Pick up position and momentum
      CLHEP::Hep3Vector pos(0,0,0);
      double time;
      _fGenerator->generatePositionAndTime(pos, time, _timeFolding);

      //Pick up energy and momentum vector
      static CLHEP::RandGeneral randEnergy       ( getEngine(), _dioSpectrum.getPDF(), _dioSpectrum.getNbins() );
      static RandomUnitSphere   randomUnitSphere ( getEngine(), _czmin, _czmax, _phimin, _phimax );

      const double e = _elow + randEnergy.fire()*(_ehi-_elow);
      const double p = safeSqrt(e*e - _mass*_mass);

      CLHEP::Hep3Vector p3 = randomUnitSphere.fire(p);

      //Set Four-momentum
      CLHEP::HepLorentzVector mom(p3,e);

      // Add the particle to  the list.
      genParts.push_back( GenParticle( PDGCode::e_minus, _dioGenId, pos, mom, time));

      if( _doHistograms ){
        const CLHEP::Hep3Vector detPos(_detSys->toDetector(pos));
        _hEElec     ->Fill( e );
        _hEElecZ    ->Fill( e );
        _hradius    ->Fill( detPos.perp() );
        _hzPosition ->Fill( detPos.z() );
        _hcz        ->Fill( p3.cosTheta() );
        _hphi       ->Fill( p3.phi() );
        _ht         ->Fill( time );
        _hmudelay   ->Fill(_fGenerator->muDelay());
        _hpulsedelay->Fill(_fGenerator->pulseDelay());
      }

    } // End of loop over particles

  } // DecayInOrbitGun::generate

}
