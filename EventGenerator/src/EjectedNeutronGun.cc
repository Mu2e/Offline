//
//
// Simulate the neutrons that come from the stopping target when muons capture
// on an Al nucleus.  Use the MECO distribution for the kinetic energy of the
// neutrons.
//
// $Id: EjectedNeutronGun.cc,v 1.10 2011/05/20 19:26:39 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/20 19:26:39 $
//
// Original author Rob Kutschke (proton gun), adapted to neutron by G. Onorato
//
//

// C++ includes.
#include <iostream>

// Framework includes
#include "art/Framework/Core/Run.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "EventGenerator/inc/EjectedNeutronGun.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eUtilities/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"
#include "TargetGeom/inc/Target.hh"

// Other external includes.
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

//ROOT Includes
#include "TH1D.h"
#include "TMath.h"

using namespace std;

static const double spectrumEndPoint = 0.1;

namespace mu2e {

  EjectedNeutronGun::EjectedNeutronGun( art::Run& run, const SimpleConfig& config ):

    // Base class.
    GeneratorBase(),

    // Configurable parameters
    _mean(config.getDouble("ejectedNeutronGun.mean",1.)),
    _elow(config.getDouble("ejectedNeutronGun.elow",0.)),
    _ehi(config.getDouble("ejectedNeutronGun.ehi",spectrumEndPoint)),
    _czmin(config.getDouble("ejectedNeutronGun.czmin",  -1.)),
    _czmax(config.getDouble("ejectedNeutronGun.czmax",  1.)),
    _phimin(config.getDouble("ejectedNeutronGun.phimin", 0. )),
    _phimax(config.getDouble("ejectedNeutronGun.phimax", CLHEP::twopi )),
    _PStoDSDelay(config.getBool("conversionGun.PStoDSDelay", true)),
    _pPulseDelay(config.getBool("conversionGun.pPulseDelay", true)),
    _nbins(config.getInt("ejectedNeutronGun.nbins",200)),
    _doHistograms(config.getBool("ejectedNeutronGun.doHistograms",true)),

    // Initialize random number distributions; getEngine comes from the base class.
    _randPoissonQ( getEngine(), std::abs(_mean) ),
    _randomUnitSphere ( getEngine(), _czmin, _czmax, _phimin, _phimax ),
    _shape ( getEngine() , &(binnedEnergySpectrum()[0]), _nbins ),

    // Histogram pointers
    _hMultiplicity(),
    _hKE(),
    _hKEZoom(),
    _hMomentumMeV(),
    _hzPosition(),
    _hcz(),
    _hphi(),
    _htime(){


    if (_nbins!=((_ehi-_elow)/0.0005)) {
      throw cet::exception("RANGE")
        << "Number f bins for the energy spectrum must be consistent with the data table binning: "
        << "0.5 KeV" ;
    }

    // About the ConditionsService:
    // The argument to the constructor is ignored for now.  It will be a
    // data base key.  There is a second argument that I have let take its
    // default value of "current"; it will be used to specify a version number.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    ConditionsHandle<DAQParams>         daqPar("ignored");
    ConditionsHandle<ParticleDataTable> pdt("ignored");

    //Set particle mass
    const HepPDT::ParticleData& p_data = pdt->particle(PDGCode::n0).ref();
    _mass = p_data.mass().value();


    // Default values for the start and end of the live window.
    // Can be overriden by the run-time config; see below.
    _tmin = daqPar->t0;
    _tmax = accPar->deBuncherPeriod;

    _tmin = config.getDouble("ejectedNeutronGun.tmin",  _tmin );
    _tmax = config.getDouble("ejectedNeutronGun.tmax",  _tmax );


    // Book histograms.
    if ( _doHistograms ){
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir  = tfs->mkdir( "EjectedNeutronGun" );
      _hMultiplicity = tfdir.make<TH1D>( "hMultiplicity", "Neutron Multiplicity",                20,     0,     20  );
      _hKE           = tfdir.make<TH1D>( "hKE",           "Neutron Kinetic Energy",              50, _elow,   _ehi  );
      _hMomentumMeV  = tfdir.make<TH1D>( "hMomentumMeV",  "Neutron Momentum in MeV",             50, _elow,   _ehi  );
      _hKEZoom       = tfdir.make<TH1D>( "hEZoom",        "Neutron Kinetic Energy (zoom)",      200, _elow,   _ehi  );
      _hzPosition    = tfdir.make<TH1D>( "hzPosition",    "Neutron z Position (Tracker Coord)", 200, -6600., -5600. );
      _hcz           = tfdir.make<TH1D>( "hcz",           "Neutron cos(theta)",                 100,    -1.,     1. );
      _hphi          = tfdir.make<TH1D>( "hphi",          "Neutron azimuth",                    100,  -M_PI,  M_PI  );
      _htime         = tfdir.make<TH1D>( "htime",         "Neutron time ",                      210,   -200.,  2000. );
    }

    _fGenerator = auto_ptr<FoilParticleGenerator>(new FoilParticleGenerator( getEngine(), _tmin, _tmax,
                                                                             FoilParticleGenerator::volWeightFoil,
                                                                             FoilParticleGenerator::flatPos,
                                                                             FoilParticleGenerator::limitedExpoTime,
                                                                             false, //dummy value
                                                                             _PStoDSDelay,
                                                                             _pPulseDelay));
  }

  EjectedNeutronGun::~EjectedNeutronGun(){
  }

  void EjectedNeutronGun::generate( GenParticleCollection& genParts ){

    // Choose the number of neutrons to generate this event.
    long n = _mean < 0 ? static_cast<long>(-_mean): _randPoissonQ.fire();
    if ( _doHistograms ) {
      _hMultiplicity->Fill(n);
    }

    //Loop over particles to generate

    for (int i=0; i<n; ++i) {

      //Pick up position and momentum
      CLHEP::Hep3Vector pos(0,0,0);
      double time;
      _fGenerator->generatePositionAndTime(pos, time);

      //Pick up energy
      double eKine = _elow + _shape.fire() * ( _ehi - _elow );
      double e   = eKine + _mass;

      //Pick up momentum vector

      _p = safeSqrt(e*e - _mass*_mass);
      CLHEP::Hep3Vector p3 = _randomUnitSphere.fire(_p);

      //Set Four-momentum
      CLHEP::HepLorentzVector mom(0,0,0,0);
      mom.setPx( p3.x() );
      mom.setPy( p3.y() );
      mom.setPz( p3.z() );
      mom.setE( e );

      // Add the particle to  the list.
      genParts.push_back( GenParticle(PDGCode::n0, GenId::ejectedNeutronGun, pos, mom, time));

      // Fill histograms.
      if ( _doHistograms) {
        _hKE->Fill( eKine );
        _hKEZoom->Fill( eKine );
        _hMomentumMeV->Fill( _p );
        _hzPosition->Fill( pos.z() );
        _hcz->Fill( mom.cosTheta() );
        _hphi->Fill( mom.phi() );
        _htime->Fill( time );

      }
    } // end of loop on particles

  } // end generate


  // Compute a binned representation of the energy spectrum of the neutron.
  //Energy in MeV
  std::vector<double> EjectedNeutronGun::binnedEnergySpectrum(){

    vector<double> neutronSpectrum;;
    ConfigFileLookupPolicy spectrumFileName;
    string NeutronFileFIP = spectrumFileName("ConditionsService/data/neutronSpectrum.txt");
    fstream infile(NeutronFileFIP.c_str(), ios::in);
    if (infile.is_open()) {
      double en, val;
      int i=0;
      while (i!=_nbins) {
        if (infile.eof()) break;
        infile >> en >> val;
        if (en >= _elow && en <= _ehi) {
          neutronSpectrum.push_back(val);
          i++;
        }
      }
    } else {
      cout << "No file associated for the ejected neutron spectrum" << endl;
      for (int i=0; i < _nbins; ++i) {
        neutronSpectrum.push_back(1);
      }
    }


    return neutronSpectrum;
  } // EjectedNeutronGun::binnedEnergySpectrum




} // namespace mu2e
