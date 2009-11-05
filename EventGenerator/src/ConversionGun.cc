//
// Generate an electron with the conversion energy
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: ConversionGun.cc,v 1.3 2009/11/05 00:12:54 rhbob Exp $ 
// $Author: rhbob $
// $Date: 2009/11/05 00:12:54 $
//
// Original author Rob Kutschke
// 

#include <iostream>
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "EventGenerator/inc/ConversionGun.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TargetGeom/inc/Target.hh"

#include "CLHEP/Random/RandFlat.h"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
using CLHEP::RandFlat;


namespace mu2e {

  // Need a home for these.  Can we use the root version?
  static const int  pdg_electron = 11;
  static const int  pdg_muon     = 13;
  static const double m = 0.000510999;

  // Also need a home for this - the cycle time of the debuncher.
  static const double tcycle = 1694.;

  ConversionGun::ConversionGun( edm::Run& run, const SimpleConfig& config ):
    GeneratorBase(){

    _doConvs = config.getBool("conversionGun.do",1);
    _p      = config.getDouble("conversionGun.p", 104.9);

    _czmin  = config.getDouble("conversionGun.czmin", 0.3);
    _czmax  = config.getDouble("conversionGun.czmax", 0.6);
    _phimin = config.getDouble("conversionGun.phimin", 0.);
    _phimax = config.getDouble("conversionGun.phimax", 2.*M_PI);
    _tmin   = config.getDouble("conversionGun.tmin",   700. );
    _tmax   = config.getDouble("conversionGun.tmax",  tcycle );

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
    int ifoil = static_cast<int>(nFoils*RandFlat::shoot());
    TargetFoil const& foil = target->foil(ifoil);

    // Foil properties.
    CLHEP::Hep3Vector const& center = foil.center();
    const double r1 = foil.rIn();
    const double dr = foil.rOut() - r1;

    // A random point within the foil.
    const double r   = r1 + dr*RandFlat::shoot();
    const double dz  = (-1.+2.*RandFlat::shoot())*foil.halfThickness();
    const double phi = 2.*M_PI*RandFlat::shoot();
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
    genParts.push_back( ToyGenParticle( pdg_electron, GenId::conversionGun, pos, mom, time));

  }

}
