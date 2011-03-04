//
// Generate an electron with the conversion energy
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: ConversionGun.cc,v 1.18 2011/03/04 23:31:34 kutschke Exp $ 
// $Author: kutschke $
// $Date: 2011/03/04 23:31:34 $
//
// Original author Rob Kutschke
// 

// C++ incldues.
#include <iostream>

// Framework includes
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes
#include "EventGenerator/inc/ConversionGun.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;

namespace mu2e {
  
  // Need a Conditions entity to hold info about conversions:  // endpoints and lifetimes for different materials etc
  // Grab them from Andrew's minimc package?
  static const double pEndPoint = 104.96;

  ConversionGun::ConversionGun( edm::Run& run, const SimpleConfig& config ):

    // Base class
    GeneratorBase(),
    

    _doConvs(config.getBool( "conversionGun.do", 1)),
    _p      (config.getDouble("conversionGun.p", pEndPoint )),
    _czmin (config.getDouble("conversionGun.czmin",  0.3)),
    _czmax (config.getDouble("conversionGun.czmax",  0.6)),
    _phimin(config.getDouble("conversionGun.phimin", 0. )),
    _phimax(config.getDouble("conversionGun.phimax", CLHEP::twopi )),
    _randomUnitSphere ( getEngine(), _czmin, _czmax, _phimin, _phimax )
   
  {

    // About the ConditionsService:
    // The argument to the constructor is ignored for now.  It will be a
    // data base key.  There is a second argument that I have let take its
    // default value of "current"; it will be used to specify a version number.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    ConditionsHandle<DAQParams>         daqPar("ignored");
    ConditionsHandle<ParticleDataTable> pdt("ignored");

    //Get particle mass
    const HepPDT::ParticleData& e_data = pdt->particle(PDGCode::e_minus).ref();
    _mass = e_data.mass().value();

    // Default values for the start and end of the live window.
    // Can be overriden by the run-time config; see below.
    _tmin = daqPar->t0;
    _tmax = accPar->deBuncherPeriod;

    _tmin = config.getDouble("conversionGun.tmin",  _tmin );
    _tmax = config.getDouble("conversionGun.tmax",  _tmax );
    
    _fGenerator = auto_ptr<FoilParticleGenerator>(new FoilParticleGenerator( getEngine(), _tmin, _tmax, 
                                                                             FoilParticleGenerator::volWeightFoil, 
                                                                             FoilParticleGenerator::flatPos, 
                                                                             FoilParticleGenerator::limitedExpoTime));
    
  }
  
  ConversionGun::~ConversionGun(){
  }
  
  void ConversionGun::generate( ToyGenParticleCollection& genParts ){
    
    //Pick up position and momentum
    double time;
    CLHEP::Hep3Vector pos(0,0,0);
    _fGenerator->generatePositionAndTime(pos, time);

    //Pick up momentum vector

    CLHEP::Hep3Vector p3 = _randomUnitSphere.fire(_p);

    //Pick up energy
    double e = sqrt( _p*_p + _mass*_mass );    

    //Set four-momentum

    CLHEP::HepLorentzVector mom(p3.x(), p3.y(), p3.z(), e);

    // Add the particle to  the list.
    genParts.push_back( ToyGenParticle( PDGCode::e_minus, GenId::conversionGun, pos, mom, time));
  }

} // end namespace mu2e
