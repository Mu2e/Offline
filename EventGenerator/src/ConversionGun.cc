//
// Generate an electron with the conversion energy
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: ConversionGun.cc,v 1.24 2011/05/18 16:11:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 16:11:17 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

// Mu2e includes
#include "EventGenerator/inc/ConversionGun.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TargetGeom/inc/zBinningForFoils.hh"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

namespace mu2e {

  // Need a Conditions entity to hold info about conversions:  // endpoints and lifetimes for different materials etc
  // Grab them from Andrew's minimc package?
  static const double pEndPoint = 104.96;

  ConversionGun::ConversionGun( art::Run& run, const SimpleConfig& config ):

    // Base class
    GeneratorBase(),

    // Get initializers from the run time configuration.
    // Initialization of tmin and tmax is deferred.
    _p     (config.getDouble("conversionGun.p", pEndPoint )),
    _czmin (config.getDouble("conversionGun.czmin", -0.5)),
    _czmax (config.getDouble("conversionGun.czmax",  0.5)),
    _phimin(config.getDouble("conversionGun.phimin", 0. )),
    _phimax(config.getDouble("conversionGun.phimax", CLHEP::twopi )),
    _PStoDSDelay(config.getBool("conversionGun.PStoDSDelay", true)),
    _pPulseDelay(config.getBool("conversionGun.pPulseDelay", true)),
    _tmin(0.),
    _tmax(0.),
    _doHistograms(config.getDouble("conversionGun.doHistograms", true )),

    // Random distribution.
    _randomUnitSphere ( getEngine(), _czmin, _czmax, _phimin, _phimax ),

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
    _hxyPos(0),
    _hrzPos(0){

    // About the ConditionsService:
    // The argument to the constructor is ignored for now.  It will be a
    // data base key.  There is a second argument that I have let take its
    // default value of "current"; it will be used to specify a version number.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    ConditionsHandle<DAQParams>         daqPar("ignored");
    ConditionsHandle<ParticleDataTable> pdt("ignored");

    // Initialize data members that could not be initialized correctly in the initiailizer list.

    // Default values for the start and end of the live window.
    _tmin = config.getDouble("conversionGun.tmin",  daqPar->t0 );
    _tmax = config.getDouble("conversionGun.tmax",  accPar->deBuncherPeriod );

    //Get particle mass
    const HepPDT::ParticleData& e_data = pdt->particle(PDGCode::e_minus).ref();
    _mass = e_data.mass().value();

    _fGenerator = auto_ptr<FoilParticleGenerator>(new FoilParticleGenerator( getEngine(), _tmin, _tmax,
                                                                             FoilParticleGenerator::volWeightFoil,
                                                                             FoilParticleGenerator::flatPos,
                                                                             FoilParticleGenerator::limitedExpoTime,
                                                                             false, //dummy value
                                                                             _PStoDSDelay,
                                                                             _pPulseDelay));
    if ( _doHistograms ) bookHistograms();
  }


  ConversionGun::~ConversionGun(){
  }

  void ConversionGun::generate( ToyGenParticleCollection& genParts ){

    // Compute position and momentum
    double time;
    CLHEP::Hep3Vector pos;
    _fGenerator->generatePositionAndTime(pos, time);

    // Compute momentum 3-vector
    CLHEP::Hep3Vector p3 = _randomUnitSphere.fire(_p);

    // Compute energy
    double e = sqrt( _p*_p + _mass*_mass );

    // Set four-momentum
    CLHEP::HepLorentzVector mom(p3, e);

    // Add the particle to  the list.
    genParts.push_back( ToyGenParticle( PDGCode::e_minus, GenId::conversionGun, pos, mom, time));

    if ( !_doHistograms ) return;

    double genRadius = pos.perp();

    _hMultiplicity->Fill(1);
    _hcz->Fill(p3.cosTheta());
    _hphi->Fill(p3.phi());
    _hmomentum->Fill(_p);
    _hradius->Fill( genRadius );
    _hzPos->Fill(pos.z());
    _htime->Fill(time);
    _hxyPos->Fill( pos.x(), pos.y()   );
    _hrzPos->Fill( pos.z(), genRadius );

  }

  void ConversionGun::bookHistograms(){

    // Compute a binning that ensures that the stopping target foils are at bin centers.
    GeomHandle<Target> target;
    Binning bins = zBinningForFoils(*target,7);
    Binning bins2 = zBinningForFoils(*target,3);

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir( "ConversionGun" );

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
                                       120, 0., 2400. );
    _hxyPos        = tfdir.make<TH2F>( "hxyPos",
                                       "Conversion Electron (x,y) at Production;(mm)",
                                       60,  -120., 120., 60, -120., 120. );
    _hrzPos        = tfdir.make<TH2F>( "hrzPos",
                                       "Conversion Electron (z,r) at Production;(mm)",
                                       bins2.nbins(), bins2.low(), bins2.high(), 60, 0., 120. );

  }

} // end namespace mu2e
