//
// Generate a positron from CLHEP::pi -> e nu
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: PiEplusNuGun.cc,v 1.15 2013/05/31 20:04:27 gandr Exp $
// $Author: gandr $
// $Date: 2013/05/31 20:04:27 $
//
// Original author Rob Kutschke heavily modified by R. Bernstein
//

// C++ includes.
#include <iostream>

// Framework includes
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "EventGenerator/inc/PiEplusNuGun.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"

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

  PiEplusNuGun::PiEplusNuGun( art::Run& run, const SimpleConfig& config ):
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

    _p      = config.getDouble("piEplusNuGun.p", pEplus );

    _czmin  = config.getDouble("piEplusNuGun.czmin",  0.3);
    _czmax  = config.getDouble("piEplusNuGun.czmax",  0.6);
    _phimin = config.getDouble("piEplusNuGun.phimin", 0. );
    _phimax = config.getDouble("piEplusNuGun.phimax", CLHEP::twopi );
    _tmin   = config.getDouble("piEplusNuGun.tmin",  _tmin );
    _tmax   = config.getDouble("piEplusNuGun.tmax",  _tmax );

  }

  PiEplusNuGun::~PiEplusNuGun(){
  }

  void PiEplusNuGun::generate( GenParticleCollection& genParts ){

    // getEngine comes from the base class.
    static CLHEP::RandFlat randFlat( getEngine() );
    static RandomUnitSphere randomUnitSphere( getEngine(), _czmin, _czmax, _phimin, _phimax );

    // Get access to the geometry system.
    GeomHandle<StoppingTarget> target;

    int nFoils = target->nFoils();

    // Pick a foil.
    int ifoil = static_cast<int>(nFoils*randFlat.fire());
    TargetFoil const& foil = target->foil(ifoil);

    CLHEP::Hep3Vector const& center = foil.centerInMu2e();
    const double r1 = foil.rIn();
    const double dr = foil.rOut() - r1;

    // A random point within the foil.
    const double r   = r1 + dr*randFlat.fire();
    const double dz  = (-1.+2.*randFlat.fire())*foil.halfThickness();
    const double phi = CLHEP::twopi*randFlat.fire();
    CLHEP::Hep3Vector pos( center.x()+r*cos(phi),
                           center.y()+r*sin(phi),
                           center.z()+dz );

    // This is not the correct distribution but it is good enough for now.
    const double time = _tmin + (_tmax-_tmin)*randFlat.fire();

    // Energy of the positron.
    double e = sqrt( _p*_p +m*m);

    // Generated 4 momentum.
    CLHEP::HepLorentzVector mom( randomUnitSphere.fire(_p), e);

    // Add the electron to  the list.
    genParts.push_back( GenParticle( PDGCode::e_plus, GenId::piEplusNuGun, pos, mom, time));

  }

}
