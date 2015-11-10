//
// Generate a (not necessarily stopped) muon
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: StoppedMuonGun.cc,v 1.6 2014/01/27 22:20:17 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/01/27 22:20:17 $
//
// Original author KLG somewhat based on ConversionGun
//

// C++ includes.
#include <iostream>

// Framework includes
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "EventGenerator/inc/StoppedMuonGun.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "StoppingTargetGeom/inc/zBinningForFoils.hh"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

namespace mu2e {

  static const double p0 = 0.0;

  StoppedMuonGun::StoppedMuonGun( art::Run& run, const SimpleConfig& config ):

    // Base class
    GeneratorBase(),

    // Get initializers from the run time configuration.
    // Initialization of tmin and tmax is deferred.
    _p     (config.getDouble("stoppedMuonGun.p", p0 )),
    _czmin (config.getDouble("stoppedMuonGun.czmin", -1.0)),
    _czmax (config.getDouble("stoppedMuonGun.czmax",  1.0)),
    _phimin(config.getDouble("stoppedMuonGun.phimin", 0. )),
    _phimax(config.getDouble("stoppedMuonGun.phimax", CLHEP::twopi )),
    _PStoDSDelay(config.getBool("stoppedMuonGun.PStoDSDelay", true)),
    _pPulseDelay(config.getBool("stoppedMuonGun.pPulseDelay", false)),
    _pPulseShift(config.getDouble("stoppedMuonGun.pPulseShift", 0)),
    _timeFolding(config.getBool("FoilParticleGenerator.foldingTimeOption", true)),
    _foilGen(config.getString("stoppedMuonGun.foilGen", "expoVolWeightFoil")),
    _posGen(config.getString("stoppedMuonGun.posGen", "flatPos")),
    _timeGen(config.getString("stoppedMuonGun.timeGen", "limitedExpoTime")),
    _tmin(0.),
    _tmax(0.),
    _doHistograms(config.getBool("stoppedMuonGun.doHistograms", true )),

    // Random distribution.
    _randomUnitSphere ( getEngine(), _czmin, _czmax, _phimin, _phimax ),
    _STfname(config.getString("FoilParticleGenerator.STfilename","ExampleDataFiles/StoppedMuons/stoppedMuons_02.txt")),
    _nToSkip (config.getInt("stoppedMuonGun.nToSkip",0)),

    // Properly initialized later.
    _mass(0),

    // Histograms
    _hMultiplicity(0),
    _hcz(0),
    _hphi(0),
    _hmomentum(0),
    _hradius(0),
    _hzPos(0),
    _htime(0),
    _hmudelay(0),
    _hpulsedelay(0),
    _hxyPos(0),
    _hrzPos(0){

    // About the ConditionsService:
    // The argument to the constructor is ignored for now.  It will be a
    // data base key.  There is a second argument that I have let take its
    // default value of "current"; it will be used to specify a version number.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    GlobalConstantsHandle<ParticleDataTable> pdt;

    // Initialize data members that could not be initialized correctly in the initiailizer list.

    // Default values for the start and end of the live window.
    _tmin = config.getDouble("stoppedMuonGun.tmin",  0. );
    _tmax = config.getDouble("stoppedMuonGun.tmax",  accPar->deBuncherPeriod );

    //Get particle mass
    const HepPDT::ParticleData& mu_data = pdt->particle(PDGCode::mu_minus).ref();
    _mass = mu_data.mass().value();

    _fGenerator = unique_ptr<FoilParticleGenerator>
      (new FoilParticleGenerator( getEngine(), _tmin, _tmax,
                                  FoilParticleGenerator::findFoilGenByName(_foilGen),
                                  FoilParticleGenerator::findPosGenByName(_posGen),
                                  FoilParticleGenerator::findTimeGenByName(_timeGen),
                                  _PStoDSDelay,
                                  _pPulseDelay,
                                  _pPulseShift,
                                  _STfname,
                                  _nToSkip));
    
    if ( _doHistograms ) bookHistograms();
  }


  StoppedMuonGun::~StoppedMuonGun(){
  }

  void StoppedMuonGun::generate( GenParticleCollection& genParts ){

    // Compute position and momentum
    double time;
    CLHEP::Hep3Vector pos;
    _fGenerator->generatePositionAndTime(pos, time, _timeFolding);

    // Compute momentum 3-vector
    CLHEP::Hep3Vector p3 = _randomUnitSphere.fire(_p);

    // Compute energy
    double e = sqrt( _p*_p + _mass*_mass );

    // Set four-momentum
    CLHEP::HepLorentzVector mom(p3, e);

    // Add the particle to  the list.
    genParts.push_back( GenParticle( PDGCode::mu_minus, GenId::stoppedMuonGun, pos, mom, time));

    if ( !_doHistograms ) return;

    double genRadius = pos.perp();

    _hMultiplicity->Fill(1);
    _hcz->Fill(p3.cosTheta());
    _hphi->Fill(p3.phi());
    _hmomentum->Fill(_p);
    _hradius->Fill( genRadius );
    _hzPos->Fill(pos.z());
    _htime->Fill(time);
    _hmudelay->Fill(_fGenerator->muDelay());
    _hpulsedelay->Fill(_fGenerator->pulseDelay());
    _hxyPos->Fill( pos.x(), pos.y()   );
    _hrzPos->Fill( pos.z(), genRadius );

  }

  void StoppedMuonGun::bookHistograms(){

    // Compute a binning that ensures that the stopping target foils are at bin centers.
    GeomHandle<StoppingTarget> target;
    Binning bins = zBinningForFoils(*target,7);
    Binning bins2 = zBinningForFoils(*target,3);

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir( "StoppedMuonGun" );

    _hMultiplicity = tfdir.make<TH1F>( "hMultiplicity", "Conversion Multiplicity",  10,  0.,  10.  );
    _hcz           = tfdir.make<TH1F>( "hcz",
                                       "Conversion Electron cos(theta) at Production;(MeV)",
                                       100,  -1.,  1.  );
    _hphi          = tfdir.make<TH1F>( "hphi",
                                       "Conversion Electron phi at Production;(MeV)",
                                       100,  -M_PI,  M_PI  );
    _hmomentum     = tfdir.make<TH1F>( "hmomentum",
                                       "Conversion Electron Momentum at Production;(MeV)",
                                       100,  90.,  110.  );
    _hradius       = tfdir.make<TH1F>( "hradius",
                                       "Conversion Electron Radius at Production;(mm)",
                                       60,  0., 120. );
    _hzPos         = tfdir.make<TH1F>( "hzPos",
                                       "Conversion Electron z at Production;(mm)",
                                       bins.nbins(), bins.low(), bins.high() );
    _htime         = tfdir.make<TH1F>( "htime",
                                       "Conversion Electron time at Production;(ns)",
                                       210, -200., 3000. );
    _hmudelay      = tfdir.make<TH1F>( "hmudelay",
                                       "Production delay due to muons arriving at ST;(ns)",
                                       300, 0., 2000. );
    _hpulsedelay   = tfdir.make<TH1F>( "hpdelay",
                                       "Production delay due to the proton pulse;(ns)",
                                       60, -150., 150. );
    _hxyPos        = tfdir.make<TH2F>( "hxyPos",
                                       "Conversion Electron (x,y) at Production;(mm)",
                                       60,  -120., 120., 60, -120., 120. );
    _hrzPos        = tfdir.make<TH2F>( "hrzPos",
                                       "Conversion Electron (z,r) at Production;(mm)",
                                       bins2.nbins(), bins2.low(), bins2.high(), 60, 0., 120. );

  }

} // end namespace mu2e
