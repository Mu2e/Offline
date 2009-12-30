//
// Generate an ejected proton from the stopping target; use the MECO distributions
// from a random spot within the target system at
// a random time during the accelerator cycle.

//
// $Id: EjectedProtonGun.cc,v 1.3 2009/12/30 19:14:10 rhbob Exp $ 
// $Author: rhbob $
// $Date: 2009/12/30 19:14:10 $
//
// Original author Rob Kutschke, heavily modified by R. Bernstein
// 

// C++ incldues.
#include <iostream>

// Framework includes
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes
#include "EventGenerator/inc/EjectedProtonGun.hh"
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

  // these all need to be in a database eventually
  // Mass of the proton.
  // Once we have the HepPDT package installed, get this number from there.
  static const double mProton = 938.272;

  EjectedProtonGun::EjectedProtonGun( edm::Run& run, const SimpleConfig& config ):
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
    
    _doEjectedProton = config.getBool( "ejectedProtonGun.do", 0);
    _ejectedProtonMomentum      = config.getDouble("ejectedProtonGun.p", .100);
    
    _czmin  = config.getDouble("ejectedProtonGun.czmin",  0.3);
    _czmax  = config.getDouble("ejectedProtonGun.czmax",  0.6);
    _phimin = config.getDouble("ejectedProtonGun.phimin", 0. );
    _phimax = config.getDouble("ejectedProtonGun.phimax", twopi );
    _tmin   = config.getDouble("ejectedProtonGun.tmin",  _tmin );
    _tmax   = config.getDouble("ejectedProtonGun.tmax",  _tmax );

    _mean = config.getDouble("ejectedProtonGun.mean",1.);
    _elow = config.getDouble("ejectedProtonGun.elow",0.);
    _ehi = config.getDouble("ejectedProtonGun.ehi",.100);
    _nbins = config.getInt("ejectedProtonGun.nbins",1000);


    _dcz  = (  _czmax -  _czmin);
    _dphi = ( _phimax - _phimin);
    _dt   = (   _tmax -   _tmin);


    double pEndPoint = _ehi; // rename here for naturalness inside, _ehi for uniformity outside

    cout << "from EjectedProtonGun: _elow, _ehi = " << _elow << " " << _ehi << endl;

    // set up the generator function
    if (_nbins>0) _bindE = (_ehi - _elow) / _nbins;
    else {
       // I'm sure this isn't the right way to do this...
       throw cms::Exception("RANGE") <<"Nonsense EjectedProtonGun.nbins requested="<<
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
	   YFunc[ib] = EnergyEjectedProtonFunc(x);
	   //	   cout << "ib, x, Spectrum from EjectedProtonGun= " << ib << " " << x << " " << EnergyEjectedProtonFunc(x) << endl;
	 }
    }
    _funcGen = auto_ptr<RandGeneral>(new RandGeneral(YFunc,_nbins));

    // Book histograms.
    edm::Service<edm::TFileService> tfs;
    _ejectedProtonMultiplicity = tfs->make<TH1D>( "ejectedProtonMultiplicity", "Ejected Proton Multiplicity", 20, 0, 20);
    _ejectedProtonKE = tfs->make<TH1D>( "ejectedProtonKE", "Ejected Proton Kinetic Energy", 10, _elow,_ehi);
    _ejectedProtonMomentumMeV = tfs->make<TH1D>( "ejectedProtonMomentumMeV", "Ejected Proton Momentum in MeV", 10, _elow,_ehi);
    _ejectedProtonKEZoom = tfs->make<TH1D>( "ejectedProtonEZoom", "Ejected Proton Kinetic Energy (zoom)", 200, _elow, _ehi);
  }
  
  EjectedProtonGun::~EjectedProtonGun(){
  }
  
  void EjectedProtonGun::generate( ToyGenParticleCollection& genParts ){
    if (!_doEjectedProton) return;
    cout << "in EjectedProton gun " << endl;//save for debugging to make sure right genconfig_ij.txt

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
    const double sz   = safeSqrt(1.- cz*cz); // what's the diff between this and boost's sqrtOrThrow?

    // Pick a random energy.  
    _ejectedProtonKineticEnergy = _elow + _funcGen->shoot() * ( _ehi - _elow );// returned in GeV, original code was in MeV
    _ejectedProtonKE->Fill(_ejectedProtonKineticEnergy);
    _ejectedProtonKEZoom->Fill(_ejectedProtonKineticEnergy);
    _ejectedProtonEnergy = _ejectedProtonKineticEnergy + mProton;

    double ejectedProtonMomentum = safeSqrt(_ejectedProtonEnergy*_ejectedProtonEnergy - mProton*mProton);
    _ejectedProtonMomentumMeV->Fill(ejectedProtonMomentum);
    
    HepLorentzVector mom( _ejectedProtonMomentum*sz*cos(phi2), ejectedProtonMomentum*sz*sin(phi2), ejectedProtonMomentum*cz, 
			  _ejectedProtonEnergy);

    // Add the proton to  the list.
    genParts.push_back( ToyGenParticle( PDGCode::p_plus, GenId::ejectedProtonGun, pos, mom, time));

  }


  double EjectedProtonGun::EnergyEjectedProtonFunc(const double protonKineticEnergyInMeV)
  {
    //taken from GMC 
//
//   Ed Hungerford  Houston University May 17 1999 
//   Rashid Djilkibaev New York University (modified) May 18 1999 
//
//   ep - proton kinetic energy (MeV)
//   pp - proton Momentum (MeV/c)
// 
//   Generates a proton spectrum similar to that observed in
//   u capture in Si.  JEPT 33(1971)11 and PRL 20(1967)569

    //these numbers are in MeV!!!!
    static const double emx = 1000.;
    static const double fac = 8.4;
    // 8.4 twice avoids stupid static problem in C++; fac cannot appear in a const
    // expression, and either I can do it here or put it in the class.  need constexpr in the standard!!
    static const double xnorm = 350.*8.4; 
    static const double emn = 1.4; // replacing par1 from GMC
    static const double par2 = 1.3279;
    static const double par3=17844.0;
    static const double par4=.32218;
    static const double par5=100.;
    static const double par6=10.014;
    static const double par7=1050.;
    static const double par8=5.103;

    double spectrumWeight;
      if (protonKineticEnergyInMeV >= 20)
	{
	  spectrumWeight=par5*TMath::Exp(-(protonKineticEnergyInMeV-20.)/par6);
	}

      else if(protonKineticEnergyInMeV >= 8.0 && protonKineticEnergyInMeV <= 20.0)
	{
	  spectrumWeight=par7*exp(-(protonKineticEnergyInMeV-8.)/par8);
	}
      else if (protonKineticEnergyInMeV > emn)
	{
	  double xw=(1.-emn/protonKineticEnergyInMeV);
	  double xu=TMath::Power(xw,par2);
	  double xv=par3*TMath::Exp(-par4*protonKineticEnergyInMeV);
	  spectrumWeight=xv*xu;
	}
      else 
	{
	  spectrumWeight = 0.;
	}


      return spectrumWeight; // GMC code was in MeV


  }// EjectedProtonGun::EnergyEjectedProtonFunc




}
