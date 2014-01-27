//
//
// Simulate the neutrons that come from the stopping target when muons capture
// on an Al nucleus.  Use the MECO distribution for the kinetic energy of the
// neutrons.
//
// $Id: EjectedNeutronGun.cc,v 1.35 2014/01/27 22:20:17 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/01/27 22:20:17 $
//
// Original author Rob Kutschke (proton gun), adapted to neutron by G. Onorato
//
//

// C++ includes.
#include <iostream>

// Framework includes
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "EventGenerator/inc/EjectedNeutronGun.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "StoppingTargetGeom/inc/zBinningForFoils.hh"

// Other external includes.
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

//ROOT Includes
#include "TH1D.h"
#include "TMath.h"

#include "cetlib/pow.h"

using namespace std;

static const double spectrumEndPoint = 100.;

namespace mu2e {

  EjectedNeutronGun::EjectedNeutronGun( art::Run& run, const SimpleConfig& config ):

    // Base class.
    GeneratorBase(),

    // Configurable parameters
    _mean(config.getDouble("ejectedNeutronGun.mean",-1.0)),
    _kelow(config.getDouble("ejectedNeutronGun.elow",0.)),
    _kehi(config.getDouble("ejectedNeutronGun.ehi",spectrumEndPoint)),
    _czmin(config.getDouble("ejectedNeutronGun.czmin",  -1.)),
    _czmax(config.getDouble("ejectedNeutronGun.czmax",  1.)),
    _phimin(config.getDouble("ejectedNeutronGun.phimin", 0. )),
    _phimax(config.getDouble("ejectedNeutronGun.phimax", CLHEP::twopi )),
    _PStoDSDelay(config.getBool("ejectedNeutronGun.PStoDSDelay", false)),
    _pPulseDelay(config.getBool("ejectedNeutronGun.pPulseDelay", false)),
    _pPulseShift(config.getDouble("ejectedNeutronGun.pPulseShift", 0)),
    _timeFolding(config.getBool("FoilParticleGenerator.foldingTimeOption", true)),
    _foilGen(config.getString("ejectedNeutronGun.foilGen", "muonFileInputFoil")),
    _posGen(config.getString("ejectedNeutronGun.posGen", "muonFileInputPos")),
    _timeGen(config.getString("ejectedNeutronGun.timeGen", "negExp")),
    _nbins(config.getInt("ejectedNeutronGun.nbins",200)),
    _tmin(config.getDouble("ejectedNeutronGun.tmin",  0. )),
    _tmax(0.),
    _doHistograms(config.getBool("ejectedNeutronGun.doHistograms",true)),
    _spectrumModel(checkSpectrumModel(config.getInt("ejectedNeutronGun.spectrumNumber",0))),
    _fGenerator(nullptr),
    _filetoread (config.getString("ejectedNeutronGun.spectrumFile","ConditionsService/data/neutronSpectrum.txt")),
    // Initialize random number distributions; getEngine comes from the base class.
    _randPoissonQ( getEngine(), std::abs(_mean) ),
    _randomUnitSphere ( getEngine(), _czmin, _czmax, _phimin, _phimax ),
    _shape ( getEngine() , &(binnedEnergySpectrum()[0]), _nbins ),
    _STfname(config.getString("FoilParticleGenerator.STfilename","ExampleDataFiles/StoppedMuons/stoppedMuons_02.txt")),
    _nToSkip (config.getInt("ejectedNeutronGun.nToSkip",0)),

    _detSys(),

     // Histogram pointers
    _hMultiplicity(),
    _hKE(),
    _hKEZoom(),
    _hMomentumMeV(),
    _hzPosition(),
    _hcz(),
    _hphi(),
    _htime(),
    _hmudelay(),
    _hpulsedelay()  {


    // About the ConditionsService:
    // The argument to the constructor is ignored for now.  It will be a
    // data base key.  There is a second argument that I have let take its
    // default value of "current"; it will be used to specify a version number.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    GlobalConstantsHandle<ParticleDataTable> pdt;

    //Set particle mass
    const HepPDT::ParticleData& p_data = pdt->particle(PDGCode::n0).ref();
    _mass = p_data.mass().value();

    // High side of the live window.
    _tmax = config.getDouble("ejectedNeutronGun.tmax",  accPar->deBuncherPeriod );

    // Limits of the momentum histotgram
    double eLow  = _kelow + _mass;
    double eHigh = _kehi  + _mass;
    double pLow  = (eLow  > _mass) ? sqrt(cet::diff_of_squares(eLow, _mass)) : 0.;
    double pHigh = (eHigh > _mass) ? sqrt(cet::diff_of_squares(eHigh,_mass)) : _mass;

    _detSys = &*GeomHandle<DetectorSystem>();

    // Book histograms.
    if ( _doHistograms ){
      // Compute a binning that ensures that the stopping target foils are at bin centers.
      GeomHandle<StoppingTarget> target;
      Binning bins = zBinningForFoils(*target,7);
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir  = tfs->mkdir( "EjectedNeutronGun" );
      _hMultiplicity = tfdir.make<TH1D>( "hMultiplicity", "Neutron Multiplicity",                20,        0,     20   );
      _hKE           = tfdir.make<TH1D>( "hKE",           "Neutron Kinetic Energy",             100,   _kelow,   _kehi  );
      _hMomentumMeV  = tfdir.make<TH1D>( "hMomentumMeV",  "Neutron Momentum in MeV",             50,     pLow,   pHigh  );
      _hKEZoom       = tfdir.make<TH1D>( "hEZoom",        "Neutron Kinetic Energy (zoom)",      200,   _kelow,   _kehi  );
      _hzPosition    = tfdir.make<TH1D>( "hzPosition",    "Neutron z Position (Tracker Coord)", bins.nbins(), bins.low(), bins.high());
      _hcz           = tfdir.make<TH1D>( "hcz",           "Neutron cos(theta)",                 100,      -1.,       1. );
      _hphi          = tfdir.make<TH1D>( "hphi",          "Neutron azimuth",                    100,    -M_PI,    M_PI  );
      _htime         = tfdir.make<TH1D>( "htime",         "Neutron time ",                      210,    -200.,    3000. );
      _hmudelay      = tfdir.make<TH1D>( "hmudelay",      "Production delay due to muons arriving at ST;(ns)", 300, 0., 2000. );
      _hpulsedelay   = tfdir.make<TH1D>( "hpdelay",       "Production delay due to the proton pulse;(ns)",      60, -150.,  150. );
    }

    _fGenerator = unique_ptr<FoilParticleGenerator>
      (new FoilParticleGenerator( getEngine(), _tmin, _tmax,
                                  FoilParticleGenerator::findFoilGenByName(_foilGen),
                                  FoilParticleGenerator::findPosGenByName(_posGen),
                                  FoilParticleGenerator::findTimeGenByName(_timeGen),
                                  _PStoDSDelay,
                                  _pPulseDelay,
                                  _pPulseShift,
                                  _STfname,
                                  _nToSkip));

  }

  EjectedNeutronGun::~EjectedNeutronGun(){
  }

  // Check bounds of spectrum number and convert to enum
  EjectedNeutronGun::SpectrumType EjectedNeutronGun::checkSpectrumModel(const int i) {

    if (i >= last_enum || i < 0) {
      throw cet::exception("RANGE")
       << "Invalid spectrum model given" ;
    }

    return SpectrumType(i);

  }

  void EjectedNeutronGun::generate( GenParticleCollection& genParts ){

    // Choose the number of neutrons to generate this event.

    long n = _mean < 0 ? static_cast<long>(-_mean):_randPoissonQ.fire();

    if ( _doHistograms ) {
      _hMultiplicity->Fill(n);
    }

    //Loop over particles to generate

    for (int i=0; i<n; ++i) {

      //Pick up position and momentum
      CLHEP::Hep3Vector pos(0,0,0);
      double time;
      _fGenerator->generatePositionAndTime(pos, time, _timeFolding);

      //Pick up energy
      double eKine = _kelow + _shape.fire() * ( _kehi - _kelow );
      double e     = eKine + _mass;


      //Pick up momentum vector
      double p = safeSqrt(e*e - _mass*_mass);
      CLHEP::Hep3Vector p3 = _randomUnitSphere.fire(p);

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
        const CLHEP::Hep3Vector detPos(_detSys->toDetector(pos));
        _hKE->Fill( eKine );
        _hKEZoom->Fill( eKine );
        _hMomentumMeV->Fill( p );
        _hzPosition->Fill( detPos.z() );
        _hcz->Fill( mom.cosTheta() );
        _hphi->Fill( mom.phi() );
        _htime->Fill( time );
        _hmudelay   ->Fill(_fGenerator->muDelay());
        _hpulsedelay->Fill(_fGenerator->pulseDelay());
      }
    } // end of loop on particles

  } // end generate


  // Continuous function of new energy spectrum from doc-db 1619
  double EjectedNeutronGun::energySpectrum(const double x) {
    // Parameters from doc-db 1619

    if (x < 10.0) {
      static const double A  = 0.0583;
      static const double B =   3.8278;
      static const double b1 =   2.1068;
      static const double b2   =   1.114;
      static const double c1s    =   pow(6.1167,2);
      static const double c2s    =  pow(2.4447,2);

      return A*(exp(-pow((x-b1),2)/c1s) + B*exp(-pow((x-b2),2)/c2s));
    }

    else {
      static const double Ah  = 0.151;
      static const double Bh =   0.0795;
      static const double c1h    =   3.8292;
      static const double c2h    =  20.0653;

      return Ah*(exp(-x/c1h) + Bh*exp(-x/c2h));
    }

  } // EjectedNeutronGun::energySpectrum


  // Compute a binned representation of the energy spectrum of the neutron.
  // Energy in MeV
  std::vector<double> EjectedNeutronGun::binnedEnergySpectrum(){


    // Vector to hold the binned representation of the energy spectrum.
    vector<double> neutronSpectrum;
    // Sanity check.
    if (_nbins <= 0) {
      throw cet::exception("RANGE")
        << "Nonsense EjectedNeutronGun.nbins requested="
        << _nbins
        << "\n";
    }

    if (_spectrumModel == docdb1619)
    {
      // Bin width.
      double dE = (_kehi - _kelow) / _nbins;
      neutronSpectrum.reserve(_nbins);

      for (int ib=0; ib<_nbins; ib++) {
        double x = (_kelow+(ib+0.5) * dE)*1000.0; //Function takes energy in MeV
        neutronSpectrum.push_back(energySpectrum(x));
      }

      return neutronSpectrum;

    }
    else
    {
      ConfigFileLookupPolicy spectrumFileName;
      string NeutronFileFIP =
        spectrumFileName(_filetoread);
      fstream infile(NeutronFileFIP.c_str(), ios::in);
      double supposedBinning = (_kehi - _kelow) / _nbins;
      if (infile.is_open()) {
        double en, val;
        bool read_on = true ;
        while (read_on) {
          if (infile.eof()) {
            en += supposedBinning;
            val = 0;
          } else {
            infile >> en >> val;
          }
          if (en >= _kelow && en <= _kehi) {
            neutronSpectrum.push_back(val);
          }
          if (en > _kehi) read_on = false;
        }
      }
      else {
        throw cet::exception("RANGE")
        << "No file associated for the ejected neutron spectrum" << endl;
      }

      return neutronSpectrum;

    }
  } // EjectedNeutronGun::binnedEnergySpectrum



} // namespace mu2e
