//
// Generate an electron with the conversion energy
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: MuonicXRayGun.cc,v 1.1 2014/02/13 17:10:07 rhbob Exp $
// $Author: rhbob $
// $Date: 2014/02/13 17:10:07 $
//
// Original author Rob Kutschke
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
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "EventGenerator/inc/MuonicXRayGun.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "StoppingTargetGeom/inc/zBinningForFoils.hh"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"
#include "TH1F.h"
#include "TH2F.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Units/PhysicalConstants.h"


using namespace std;

namespace mu2e {

  MuonicXRayGun::MuonicXRayGun( art::Run& run, const SimpleConfig& config ):

    // Base class
    GeneratorBase(),

    // Get initializers from the run time configuration.
    // Initialization of tmin and tmax is deferred.
    _czmin (config.getDouble("muonicXRayGun.czmin", -1.0)),
    _czmax (config.getDouble("muonicXRayGun.czmax",  1.0)),
    _phimin(config.getDouble("muonicXRayGun.phimin", 0. )),
    _phimax(config.getDouble("muonicXRayGun.phimax", CLHEP::twopi )),
    _PStoDSDelay(config.getBool("muonicXRayGun.PStoDSDelay", true)),
    _pPulseDelay(config.getBool("muonicXRayGun.pPulseDelay", false)),
    _pPulseShift(config.getDouble("muonicXRayGun.pPulseShift", 0)),
    _timeFolding(config.getBool("FoilParticleGenerator.foldingTimeOption", true)),
    _tmin(0.0),
    _tmax(0.0),
    _foilGen(config.getString("muonicXRayGun.foilGen", "muonFileInputFoil")),
    _posGen(config.getString("muonicXRayGun.posGen", "muonFileInputPos")),
    _timeGen(config.getString("muonicXRayGun.timeGen", "negExp")),
    _doHistograms(config.getBool("muonicXRayGun.doHistograms", true )),

    // Random distribution. 
    _randomUnitSphere ( getEngine(), _czmin, _czmax, _phimin, _phimax ),
    _randFlat( getEngine()),
    _STfname(config.getString("FoilParticleGenerator.STfilename","ExampleDataFiles/StoppedMuons/stoppedMuons_02.txt")),
    _nToSkip (config.getInt("muonicXRayGun.nToSkip",0)),

    // Properly initialized later.
  
    _detSys(),

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
    _tmin = config.getDouble("muonicXRayGun.tmin",  0. );
    _tmax = config.getDouble("muonicXRayGun.tmax",  accPar->deBuncherPeriod );


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

    _detSys = &*GeomHandle<DetectorSystem>();

    if ( _doHistograms ) bookHistograms();
  }


  MuonicXRayGun::~MuonicXRayGun(){
  }

  void MuonicXRayGun::generate( GenParticleCollection& genParts ){

    // Compute position and momentum
    double time;
    CLHEP::Hep3Vector pos;
    _fGenerator->generatePositionAndTime(pos, time, _timeFolding);

    
    //       double _xtemp = -3904.;
    //double _ytemp = 0;
    //double _ztemp = pos.z() + 2000.;
    //pos = CLHEP::Hep3Vector(_xtemp,_ytemp,_ztemp);
    
 // length of foil system is 17*50, double it.  total hack.
    //   std::cout << "pos = " << pos << std::endl;

    int nphotons = 0;
    vector<double> photonEnergy;
    
    // create X-Rays.  Measday tells us values for Al.  Table 3.5

    
    double prob = _randFlat.fire();

    //    cout << "prob = " << prob << std::endl;
    if (prob < 0.625)  // 3d-2p line
      {
	++nphotons;
	photonEnergy.push_back(0.0661);
      }

    prob = _randFlat.fire();
    if (prob < 0.797)  // 2p-1s line
      {
	++nphotons;
	photonEnergy.push_back(0.3468);
      }
    
    for (int ithphoton=0; ithphoton < nphotons; ++ithphoton)
      {
	// Compute momentum 3-vector
	CLHEP::Hep3Vector p3 = _randomUnitSphere.fire(photonEnergy[ithphoton]);

    // Compute energy
	_p = photonEnergy[ithphoton];
	double e = _p; // yes this is stupid, let the optimizer fix it.  keeps code parallel among guns

    // Set four-momentum
	CLHEP::HepLorentzVector mom(p3, e);

    // Add the particle(s) to  the list.
	genParts.push_back( GenParticle( PDGCode::gamma, GenId::muonicXRayGun, pos, mom, time));

      
	if (_doHistograms )
	  {

	    const CLHEP::Hep3Vector detPos(_detSys->toDetector(pos));
	    double genRadius = detPos.perp();
	    _hMultiplicity->Fill(nphotons);
	    _hcz->Fill(p3.cosTheta());
	    _hphi->Fill(p3.phi());
 	    _hmomentum->Fill(_p*1000.);
	    //	    _hmomentum->Fill(_p);
	    _hradius->Fill( genRadius );
	    _hzPos->Fill(detPos.z());
	    _htime->Fill(time);
	    _hmudelay->Fill(_fGenerator->muDelay());
	    _hpulsedelay->Fill(_fGenerator->pulseDelay());
	    _hxyPos->Fill( detPos.x(), detPos.y()   );
	    _hrzPos->Fill( detPos.z(), genRadius );
	  }
      }
  }

  void MuonicXRayGun::bookHistograms(){

    // Compute a binning that ensures that the stopping target foils are at bin centers.
    GeomHandle<StoppingTarget> target;
    Binning bins = zBinningForFoils(*target,7);
    Binning bins2 = zBinningForFoils(*target,3);

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir( "MuonicXRayGun" );

    _hMultiplicity = tfdir.make<TH1F>( "hMultiplicity", "MuonicXRay Multiplicity",  10,  0.,  10.  );
    _hcz           = tfdir.make<TH1F>( "hcz",
                                       "MuonicXRay Photon cos(theta) at Production;(MeV)",
                                       100,  -2.,  2.  );
    _hphi          = tfdir.make<TH1F>( "hphi",
                                       "MuonicXRay Photon phi at Production;(MeV)",
                                       100,  -M_PI,  M_PI  );
    _hmomentum     = tfdir.make<TH1F>( "hmomentum",
                                       "MuonicXRay Photon Momentum at Production;(keV)",
                                       500,  0.,  500.  );
    _hradius       = tfdir.make<TH1F>( "hradius",
                                       "MuonicXRay Photon Radius at Production;(mm)",
                                       60,  0., 120. );
    _hzPos         = tfdir.make<TH1F>( "hzPos",
                                       "MuonicXRay Photon z at Production;(mm)",
                                       bins.nbins(), bins.low(), bins.high() );
    _htime         = tfdir.make<TH1F>( "htime",
                                       "MuonicXRay Photon time at Production;(ns)",
                                       210, -200., 3000. );
    _hmudelay      = tfdir.make<TH1F>( "hmudelay",
                                       "Production delay due to muons arriving at ST;(ns)",
                                       300, 0., 2000. );
    _hpulsedelay   = tfdir.make<TH1F>( "hpdelay",
                                       "Production delay due to the proton pulse;(ns)",
                                       60, 0., 300. );
    _hxyPos        = tfdir.make<TH2F>( "hxyPos",
                                       "MuonicXRay Photon (x,y) at Production;(mm)",
                                       60,  -120., 120., 60, -120., 120. );
    _hrzPos        = tfdir.make<TH2F>( "hrzPos",
                                       "MuonicXRay Photon (z,r) at Production;(mm)",
                                       bins2.nbins(), bins2.low(), bins2.high(), 60, 0., 120. );

  }

} // end namespace mu2e
