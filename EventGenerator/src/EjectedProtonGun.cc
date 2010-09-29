//
//
// Simulate the protons that come from the stopping target when muons capture
// on an Al nucleus.  Use the MECO distribution for the kinetic energy of the
// protons.  Production is uniform across the targets and uniform in time;
// this model needs to be improved.
//
// $Id: EjectedProtonGun.cc,v 1.8 2010/09/29 22:55:06 kutschke Exp $ 
// $Author: kutschke $
// $Date: 2010/09/29 22:55:06 $
//
// Original author Rob Kutschke, heavily modified by R. Bernstein
// 
// 1) This code uses a incorrect model of the distribution of proton production over 
//    the targets.  It is uniform in target number and uniform across each target.
//    At a future date this needs to be made more realistic.
// 2) This code uses an incorrect model of the distribution of protons in time.
//    At a future date this needs to be made more realistic.
// 3) About the initialization of _shape.
//    The c'tor of RandGeneral wants, as its second argument, the starting
//    address of an array of doubles that describes the required shape.
//    The method binnedEnergySpectrum returns, by value, a std::vector<double>.
//    We can get the required argument by taking the address of the first element 
//    of the std::vector. There is a subtlety about the return value of
//    binnedEnergySpectrum:  it returns by value to a temporary variable that
//    we cannot see; this variable goes out of scope after the c'tor completes;
//    therefore its lifetime is managed properly.
//

// C++ incldues.
#include <iostream>

// Framework includes
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"

// Mu2e includes
#include "EventGenerator/inc/EjectedProtonGun.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "TargetGeom/inc/Target.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"

// General Utilities
#include "GeneralUtilities/inc/pow.hh"

// Other external includes.
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

//ROOT Includes
#include "TH1D.h"
#include "TMath.h"

using namespace std;

namespace mu2e {

  // these all need to be in a database eventually
  // Mass of the proton.
  // Once we have the HepPDT package installed, get this number from there.
  //static const double mProton = 938.272;

  EjectedProtonGun::EjectedProtonGun( edm::Run& run, const SimpleConfig& config ):

    // Base class.
    GeneratorBase(),
    
    // Configurable parameters
    _mean(config.getDouble("ejectedProtonGun.mean",1.)),
    _elow(config.getDouble("ejectedProtonGun.elow",0.)),
    _ehi(config.getDouble("ejectedProtonGun.ehi",.100)),
    _czmin(config.getDouble("ejectedProtonGun.czmin",  0.3)),
    _czmax(config.getDouble("ejectedProtonGun.czmax",  0.6)),
    _phimin(config.getDouble("ejectedProtonGun.phimin", 0. )),
    _phimax(config.getDouble("ejectedProtonGun.phimax", CLHEP::twopi )),
    _nbins(config.getInt("ejectedProtonGun.nbins",1000)),
    _doHistograms(config.getBool("ejectedProtonGun.doHistograms",true)),

    // Initialize random number distributions; getEngine comes from the base class.
    _randFlat( getEngine() ),
    _randPoissonQ( getEngine(), std::abs(_mean) ),
    _randomUnitSphere( getEngine(), _czmin, _czmax, _phimin, _phimax),

    // See Note 3.
    _shape( GeneratorBase::getEngine(), &(binnedEnergySpectrum()[0]), _nbins),

    // Histogram pointers
    _hMultiplicity(),
    _hKE(),
    _hKEZoom(),
    _hMomentumMeV(),
    _hzPosition(),
    _hcz(),
    _hphi(),
    _htime(){

    // About the ConditionsService:
    // The argument to the constructor is ignored for now.  It will be a
    // data base key.  There is a second argument that I have let take its
    // default value of "current"; it will be used to specify a version number.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    ConditionsHandle<DAQParams>         daqPar("ignored");
    ConditionsHandle<ParticleDataTable> pdt("ignored");
    
    // Default values for the start and end of the live window.
    // Can be overriden by the run-time config; see below.
    _tmin = daqPar->t0;
    _tmax = accPar->deBuncherPeriod;
    
    _tmin = config.getDouble("ejectedProtonGun.tmin",  _tmin );
    _tmax = config.getDouble("ejectedProtonGun.tmax",  _tmax );

    // Get the electron mass from the particle data table (in MeV).
    const HepPDT::ParticleData& e_data = pdt->particle(PDGCode::p_plus);
    _mProton = e_data.mass().value();

    // Book histograms.
    if ( _doHistograms ){
      edm::Service<edm::TFileService> tfs;
      edm::TFileDirectory tfdir  = tfs->mkdir( "EjectedProtonGun" );
      _hMultiplicity = tfdir.make<TH1D>( "hMultiplicity", "Proton Multiplicity",                20,     0,     20  );
      _hKE           = tfdir.make<TH1D>( "hKE",           "Proton Kinetic Energy",              50, _elow,   _ehi  );
      _hMomentumMeV  = tfdir.make<TH1D>( "hMomentumMeV",  "Proton Momentum in MeV",             50, _elow,   _ehi  );
      _hKEZoom       = tfdir.make<TH1D>( "hEZoom",        "Proton Kinetic Energy (zoom)",      200, _elow,   _ehi  );
      _hzPosition    = tfdir.make<TH1D>( "hzPosition",    "Proton z Position (Tracker Coord)", 200, -6600., -5600. );
      _hcz           = tfdir.make<TH1D>( "hcz",           "Proton cos(theta)",                 100,    -1.,     1. );
      _hphi          = tfdir.make<TH1D>( "hphi",          "Proton azimuth",                    100,  -M_PI,  M_PI  );
      _htime         = tfdir.make<TH1D>( "htime",         "Proton time ",                      200,      0,  2000. );
    }
  }
  
  EjectedProtonGun::~EjectedProtonGun(){
  }
  
  void EjectedProtonGun::generate( ToyGenParticleCollection& genParts ){

    // Choose the number of electrons to generate this event.
    long n = _mean < 0 ? static_cast<long>(-_mean): _randPoissonQ.fire();
    if ( _doHistograms ) { 
      _hMultiplicity->Fill(n); 
    }

    // Get information about the target system.
    GeomHandle<Target> target;
    int nFoils = target->nFoils();

    for ( long i=0; i<n; ++i ){

      // Pick a foil.
      int ifoil = static_cast<int>(nFoils*_randFlat.fire());
      TargetFoil const& foil = target->foil(ifoil);
      
      // Foil properties.
      CLHEP::Hep3Vector const& center = foil.center();
      const double r1 = foil.rIn();
      const double dr = foil.rOut() - r1;

      // Pick a random point within the foil.
      const double r   = r1 + dr*_randFlat.fire();
      const double dz  = (-1.+2.*_randFlat.fire())*foil.halfThickness();
      const double phi = CLHEP::twopi*_randFlat.fire();
      CLHEP::Hep3Vector pos( center.x()+r*cos(phi), 
                             center.y()+r*sin(phi), 
                             center.z()+dz );

      // Pick a kinetic energy from the distribution.  Compute energy and momentum.
      double eKine = _elow + _shape.fire() * ( _ehi - _elow );
      double ePr   = eKine + _mProton;
      double pPr   = safeSqrt(ePr*ePr - _mProton*_mProton);

      // Pick a 4 vector uniformly over the requested region of the unit sphere.
      CLHEP::HepLorentzVector mom( _randomUnitSphere.fire(pPr), ePr);

      // This is not the right time structure.
      double time = _tmin + (_tmax - _tmin)*_randFlat.fire();

      // Add the proton to the list of generated particles.
      genParts.push_back( ToyGenParticle( PDGCode::p_plus, GenId::ejectedProtonGun, pos, mom, time));

      cout << "Ejected proton: " 
           << mom  << " "
           << pPr << " "
           << mom.vect().mag() << " "
           << time
           << endl;

      // Fill histograms.
      if ( _doHistograms) {
        _hKE->Fill(eKine);
        _hKEZoom->Fill(eKine);
        _hMomentumMeV->Fill(pPr);
        _hzPosition->Fill(pos.z());
        _hcz->Fill(mom.cosTheta());
        _hphi->Fill(mom.phi());
        _htime->Fill(time);
      }

    } // end loop over generated protons.

} // end generate

  // Input argument is in MeV.
  double EjectedProtonGun::energySpectrum(const double protonKineticEnergyInMeV){

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
    
    
    return spectrumWeight;

  } // EjectedProtonGun::EnergyEjectedProtonFunc

  // Compute a binned representation of the energy spectrum.
  std::vector<double> EjectedProtonGun::binnedEnergySpectrum(){

    // Sanity check.
    if (_nbins <= 0) {
      throw cms::Exception("RANGE") 
        << "Nonsense ejectedProtonGun.nbins requested="
        << _nbins
        << "\n";
    }

    // Bin width.
    double dE = (_ehi - _elow) / _nbins;

    // Vector to hold the binned representation of the energy spectrum.
    std::vector<double> spectrum;
    spectrum.reserve(_nbins);

    for (int ib=0; ib<_nbins; ib++) {
      double e = _elow+(ib+0.5) * dE;
      spectrum.push_back(energySpectrum(e));
    }

    return spectrum;
  } // EjectedProtonGun::binnedEnergySpectrum(){

} // namespace mu2e
