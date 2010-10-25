//
// Generate an electron with the conversion energy
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: ConversionGun.cc,v 1.14 2010/10/25 19:50:21 onoratog Exp $ 
// $Author: onoratog $
// $Date: 2010/10/25 19:50:21 $
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
#include "EventGenerator/inc/FoilParticleGenerator.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;

namespace mu2e {
  
  // Need a Conditions entity to hold info about conversions:
  // endpoints and lifetimes for different materials etc
  // Grab them from Andrew's minimc package?
  static const double pEndPoint = 104.96;

  ConversionGun::ConversionGun( edm::Run& run, const SimpleConfig& config ):

    // Base class
    GeneratorBase(){

    // About the ConditionsService:
    // The argument to the constructor is ignored for now.  It will be a
    // data base key.  There is a second argument that I have let take its
    // default value of "current"; it will be used to specify a version number.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    ConditionsHandle<DAQParams>         daqPar("ignored");

    // Default values for the start and end of the live window.
    // Can be overriden by the run-time config; see below.
    double _tmin = daqPar->t0;
    double _tmax = accPar->deBuncherPeriod;
    
    _doConvs = config.getBool( "conversionGun.do", 1);

    _p      = config.getDouble("conversionGun.p", pEndPoint );
    _czmin  = config.getDouble("conversionGun.czmin",  0.3);
    _czmax  = config.getDouble("conversionGun.czmax",  0.6);
    _phimin = config.getDouble("conversionGun.phimin", 0. );
    _phimax = config.getDouble("conversionGun.phimax", CLHEP::twopi );
    _tmin   = config.getDouble("conversionGun.tmin",  _tmin );
    _tmax   = config.getDouble("conversionGun.tmax",  _tmax );

    cout << "in conversion gun class the values are:\ntmin = " 
         << _tmin << "\ttmax = " << _tmax << endl;

    FoilParticleGenerator generator( PDGCode::e_minus, GenId::conversionGun );

    fGenerator = generator;
    fGenerator._czmin = _czmin;
    fGenerator._czmax = _czmax;
    fGenerator._phimin = _phimin;
    fGenerator._phimax = _phimax;
    fGenerator._tmin = _tmin;
    fGenerator._tmax = _tmax;
    fGenerator._p = _p;

    fGenerator.setRandomEngine( GeneratorBase::getEngine() );

  }
  
  ConversionGun::~ConversionGun(){
  }
  
  void ConversionGun::generate( ToyGenParticleCollection& genParts ){
    
    if (!_doConvs) return;
    fGenerator.generateFromFoil(genParts);

  }

} // end namespace mu2e
