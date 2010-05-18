//
// Generate an electron with the conversion energy
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: DecayInOrbitGun.cc,v 1.6 2010/05/18 21:15:33 kutschke Exp $ 
// $Author: kutschke $
// $Date: 2010/05/18 21:15:33 $
//
// Original author Rob Kutschke
// 

// C++ incldues.
#include <iostream>

// Framework includes
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes
#include "EventGenerator/inc/DecayInOrbitGun.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "TargetGeom/inc/Target.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"

// General Utilities
#include "GeneralUtilities/inc/pow.hh"

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
  static const double mElectron = 0.510999;
  
  // Need a Conditions entity to hold info about conversions:
  // endpoints and lifetimes for different materials etc
  // Grab them from Andrew's minimc package?

  static const double conversionEnergyAluminum = 104.96;
  DecayInOrbitGun::DecayInOrbitGun( edm::Run& run, const SimpleConfig& config ):
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
    
    _doConvs = config.getBool( "decayinorbitGun.do", 0);
    _p      = config.getDouble("decayinorbitGun.p", conversionEnergyAluminum );
    
    _czmin  = config.getDouble("decayinorbitGun.czmin",  0.3);
    _czmax  = config.getDouble("decayinorbitGun.czmax",  0.6);
    _phimin = config.getDouble("decayinorbitGun.phimin", 0. );
    _phimax = config.getDouble("decayinorbitGun.phimax", CLHEP::twopi );
    _tmin   = config.getDouble("decayinorbitGun.tmin",  _tmin );
    _tmax   = config.getDouble("decayinorbitGun.tmax",  _tmax );

    _mean = config.getDouble("decayinorbitGun.mean",1.);
    _elow = config.getDouble("decayinorbitGun.elow",.100);
    _ehi = config.getDouble("decayinorbitGun.ehi",conversionEnergyAluminum);
    _nbins = config.getInt("decayinorbitGun.nbins",1000);


    _dcz  = (  _czmax -  _czmin);
    _dphi = ( _phimax - _phimin);
    _dt   = (   _tmax -   _tmin);


    double pEndPoint = _ehi; // rename here for naturalness inside, _ehi for uniformity outside


    // set up the generator function
    if (_nbins>0) _bindE = (_ehi - _elow) / _nbins;
    else {
      // I'm sure this isn't the right way to do this...
      throw cms::Exception("RANGE") <<"Nonsense DecayInOrbitGun.nbins requested="<<
        _nbins<<"\n";
    }

    double YFunc[_nbins];
    for (int ib=0; ib<_nbins; ib++) {

      double x = _elow+(ib+0.5) * _bindE;
      if (x > pEndPoint)
        {
          cout << "past endpoint " << endl;
          YFunc[ib] = 0.;
        }
      else
        {
          YFunc[ib] = EnergyDIOFunc(x);
          //           cout << "ib, x, Spectrum = " << ib << " " << x << " " << EnergyDIOFunc(x) << endl;
        }
    }
    _funcGen = auto_ptr<CLHEP::RandGeneral>(new CLHEP::RandGeneral(YFunc,_nbins));

    // Book histograms.
    edm::Service<edm::TFileService> tfs;
    _decayInOrbitMultiplicity = tfs->make<TH1D>( "decayInOrbitMultiplicity", "Decay In Orbit Multiplicity", 20, 0, 20);
    _decayInOrbitEElec = tfs->make<TH1D>( "decayInOrbitEElec", "Decay In Orbit Electron Energy", 10, _elow,0.105);
    _decayInOrbitEElecZ = tfs->make<TH1D>( "decayInOrbitEElecZ", "Decay In Orbit Electron Energy (zoom)", 200, _elow, _ehi);




    
  }
  
  DecayInOrbitGun::~DecayInOrbitGun(){
  }
  
  void DecayInOrbitGun::generate( ToyGenParticleCollection& genParts ){
    if (!_doConvs) return;
    cout << "in DIO gun " << endl;//save for debugging to make sure right genconfig_ij.txt

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
    const double sz   = safeSqrt(1.- cz*cz); // what's the diff between this and boost's sqrtOrThrow?

    // Pick a random energy.  
    const double e = _elow + _funcGen->shoot() * ( _ehi - _elow );
    _decayInOrbitEElec->Fill(e);
    _decayInOrbitEElecZ->Fill(e);

    double electronMomentum = safeSqrt(e*e - mElectron*mElectron);
    
    cout << "electron DIO momentum = " << electronMomentum << endl;
    CLHEP::HepLorentzVector mom( electronMomentum*sz*cos(phi2), electronMomentum*sz*sin(phi2), electronMomentum*cz, e);

    // Add the electron to  the list.
    genParts.push_back( ToyGenParticle( PDGCode::e_minus, GenId::dio1, pos, mom, time));

  }


  const double DecayInOrbitGun::EnergyDIOFunc(const double x)
  {
    return pow<5>(conversionEnergyAluminum - x) ;


  }// DecayInOrbitGun::EnergyDIOFunc




}
