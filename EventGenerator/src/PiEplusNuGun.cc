//
// Generate a positron from pi -> e nu
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: PiEplusNuGun.cc,v 1.1 2009/12/22 17:29:25 rhbob Exp $ 
// $Author: rhbob $
// $Date: 2009/12/22 17:29:25 $
//
// Original author Rob Kutschke heavily modified by R. Bernstein
// 

// C++ incldues.
#include <iostream>

// Framework includes
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes
#include "EventGenerator/inc/PiEplusNuGun.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "TargetGeom/inc/Target.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"

// Other external includes.
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
using CLHEP::RandFlat;
using CLHEP::twopi;

namespace mu2e {

  // Mass of the electron.
  // Once we have the HepPDT package installed, get this number from there.
  static const double m = 0.510999;
  
  // Need a Conditions entity to hold info about conversions:
  // endpoints and lifetimes for different materials etc
  // Grab them from Andrew's minimc package?
  static const double pEplus = 70.0;

  PiEplusNuGun::PiEplusNuGun( edm::Run& run, const SimpleConfig& config ):
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
    
    _doPiEplusNu = config.getBool( "piEplusNuGun.do", 0);
    _p      = config.getDouble("piEplusNuGun.p", pEplus );
    
    _czmin  = config.getDouble("piEplusNuGun.czmin",  0.3);
    _czmax  = config.getDouble("piEplusNuGun.czmax",  0.6);
    _phimin = config.getDouble("piEplusNuGun.phimin", 0. );
    _phimax = config.getDouble("piEplusNuGun.phimax", twopi );
    _tmin   = config.getDouble("piEplusNuGun.tmin",  _tmin );
    _tmax   = config.getDouble("piEplusNuGun.tmax",  _tmax );

    _dcz  = (  _czmax -  _czmin);
    _dphi = ( _phimax - _phimin);
    _dt   = (   _tmax -   _tmin);
    
  }
  
  PiEplusNuGun::~PiEplusNuGun(){
  }
  
  void PiEplusNuGun::generate( ToyGenParticleCollection& genParts ){
    
    if (!_doPiEplusNu) return;

    // Get access to the geometry system.
    GeomHandle<Target> target;
    
    int nFoils = target->nFoils();
    
    // Pick a foil.
    int ifoil = static_cast<int>(nFoils*RandFlat::shoot());
    TargetFoil const& foil = target->foil(ifoil);

    // Foil properties.
    CLHEP::Hep3Vector const& center = foil.center();
    const double r1 = foil.rIn();
    const double dr = foil.rOut() - r1;
    
    // A random point within the foil.
    const double r   = r1 + dr*RandFlat::shoot();
    const double dz  = (-1.+2.*RandFlat::shoot())*foil.halfThickness();
    const double phi = twopi*RandFlat::shoot();
    Hep3Vector pos( center.x()+r*cos(phi), 
		    center.y()+r*sin(phi), 
		    center.z()+dz );
    
    // Random direction.
    // Replace this with RandomUnitSphere from Mu2eUtilities/inc
    const double cz   = _czmin  +  _dcz*RandFlat::shoot();
    const double phi2 = _phimin + _dphi*RandFlat::shoot();
    
    // This should be an exponential decay.
    const double time = _tmin   +   _dt*RandFlat::shoot();

    // Derived quantities.
    const double sz   = safeSqrt(1.- cz*cz);
    double e = sqrt( _p*_p +m*m);
    
    HepLorentzVector mom( _p*sz*cos(phi2), _p*sz*sin(phi2), _p*cz, e);

    // Add the electron to  the list.
    genParts.push_back( ToyGenParticle( PDGCode::e_plus, GenId::piEplusNuGun, pos, mom, time));

  }

}
