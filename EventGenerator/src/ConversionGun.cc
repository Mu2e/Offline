//
// Generate an electron with the conversion energy
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: ConversionGun.cc,v 1.11 2010/05/18 21:15:30 kutschke Exp $ 
// $Author: kutschke $
// $Date: 2010/05/18 21:15:30 $
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
#include "Mu2eUtilities/inc/safeSqrt.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
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
  
  // Need a Conditions entity to hold info about conversions:
  // endpoints and lifetimes for different materials etc
  // Grab them from Andrew's minimc package?
  static const double pEndPoint = 104.96;

  ConversionGun::ConversionGun( edm::Run& run, const SimpleConfig& config ):
    GeneratorBase(){

    // About the ConditionsService:
    // The argument to the constructor is ignored for now.  It will be a
    // data base key.  There is a second argument that I have let take its
    // default value of "current"; it will be used to specify a version number.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    ConditionsHandle<DAQParams>         daqPar("ignored");
    ConditionsHandle<ParticleDataTable> pdt("ignored");

    // Get the electron mass from the particle data table (in MeV).
    const HepPDT::ParticleData& e_data = pdt->particle(PDGCode::e_minus);
    _mass = e_data.mass().value();

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

    _dcz  = (  _czmax -  _czmin);
    _dphi = ( _phimax - _phimin);
    _dt   = (   _tmax -   _tmin);
    
  }
  
  ConversionGun::~ConversionGun(){
  }
  
  void ConversionGun::generate( ToyGenParticleCollection& genParts ){
    
    if (!_doConvs) return;

    // Get access to the geometry system.
    GeomHandle<Target> target;
    
    int nFoils = target->nFoils();
    
    // Pick a foil.
    int ifoil = static_cast<int>(nFoils*CLHEP::RandFlat::shoot());
    TargetFoil const& foil = target->foil(ifoil);

    // Foil properties.
    CLHEP::Hep3Vector const& center = foil.center();
    const double r1 = foil.rIn();
    const double dr = foil.rOut() - r1;
    
    // A random point within the foil.
    const double r   = r1 + dr*CLHEP::RandFlat::shoot();
    const double dz  = (-1.+2.*CLHEP::RandFlat::shoot())*foil.halfThickness();
    const double phi = CLHEP::twopi*CLHEP::RandFlat::shoot();
    CLHEP::Hep3Vector pos( center.x()+r*cos(phi), 
                           center.y()+r*sin(phi), 
                           center.z()+dz );
    
    // Random direction.
    // Replace this with RandomUnitSphere from Mu2eUtilities/inc
    const double cz   = _czmin  +  _dcz*CLHEP::RandFlat::shoot();
    const double phi2 = _phimin + _dphi*CLHEP::RandFlat::shoot();
    
    // This should be an exponential decay.
    const double time = _tmin   +   _dt*CLHEP::RandFlat::shoot();

    // Derived quantities.
    const double sz   = safeSqrt(1.- cz*cz);
    double e = sqrt( _p*_p + _mass*_mass );
    
    CLHEP::HepLorentzVector mom( _p*sz*cos(phi2), _p*sz*sin(phi2), _p*cz, e);

    // Add the electron to  the list.
    genParts.push_back( ToyGenParticle( PDGCode::e_minus, GenId::conversionGun, pos, mom, time));

  }

} // end namespace mu2e
